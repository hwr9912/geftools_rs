use ahash::RandomState;
use anyhow::*;
use clap::Parser;
use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use hashbrown::HashMap;
use hdf5::types::VarLenUnicode;
use hdf5::File as H5File;
use itertools::Itertools;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

type Key = u64; // (bx<<32)|by

#[derive(Parser)]
struct Args {
    /// 输入 GEM 或 GEM.GZ
    #[arg(short, long)]
    input: String,
    /// 输出 bGEF (HDF5)
    #[arg(short, long)]
    output: String,
    /// 逗号分隔的 bin 列表
    #[arg(short, long, value_delimiter = ',', default_value = "1,20,50,100")]
    bins: Vec<u32>,
    /// 可选区域：minx,maxx,miny,maxy
    #[arg(long)]
    region: Option<String>,
    /// 分块大小（行）
    #[arg(long, default_value_t = 2_000_000)]
    chunksize: usize,
    /// 顶层属性：resolution
    #[arg(long, default_value_t = 1)]
    resolution: i32,
    /// 顶层属性：omics
    #[arg(long, default_value = "Transcriptomics")]
    omics: String,
}

fn open_text(path: &str) -> Result<Box<dyn Read>> {
    let f = File::open(path).with_context(|| format!("open {}", path))?;
    if path.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}

#[derive(Debug, Clone, Copy)]
struct Header {
    offset_x: i32,
    offset_y: i32,
    has_exon: bool,
    header_line_index: usize,
}

fn parse_header(path: &str) -> Result<Header> {
    let rdr = open_text(path)?;
    let mut br = BufReader::new(rdr);
    let mut line = String::new();
    let (mut ox, mut oy) = (0, 0);
    let mut header_line: Option<(usize, String)> = None;
    let mut i = 0usize;
    loop {
        line.clear();
        let n = br.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if let Some(rest) = line.strip_prefix("#OffsetX=") {
            ox = rest.trim().parse()?;
        } else if let Some(rest) = line.strip_prefix("#OffsetY=") {
            oy = rest.trim().parse()?;
        } else if line.starts_with("geneID") {
            header_line = Some((i, line.clone()));
            break;
        }
        i += 1;
    }
    let (idx, header) = header_line.ok_or_else(|| anyhow!("未找到表头 geneID"))?;
    let cols = header.matches('\t').count() + 1;
    Ok(Header {
        offset_x: ox,
        offset_y: oy,
        has_exon: cols == 5,
        header_line_index: idx,
    })
}

#[derive(Clone, Copy, Default)]
struct Agg {
    mid: u32,
    exon: u32,
}

#[derive(Clone, Copy)]
struct Range {
    min_x: u32,
    max_x: u32,
    min_y: u32,
    max_y: u32,
}

fn key_xy(x: u32, y: u32, bin: u32) -> Key {
    let bx = x / bin;
    let by = y / bin;
    ((bx as u64) << 32) | (by as u64)
}

fn parse_region(s: &str) -> Result<[u32; 4]> {
    let v: Vec<u32> = s.split(',').map(|t| t.trim().parse()).try_collect()?;
    if v.len() != 4 {
        bail!("region 需为 minx,maxx,miny,maxy");
    }
    Ok([v[0], v[1], v[2], v[3]])
}

/// 第一遍：仅求整体范围（可选；这里只是打印参考）
fn first_pass_range(path: &str, skip_header: usize) -> Result<Range> {
    let rdr = open_text(path)?;
    let mut br = BufReader::new(rdr);
    let mut line = String::new();
    let mut i = 0usize;
    let (mut min_x, mut min_y) = (u32::MAX, u32::MAX);
    let (mut max_x, mut max_y) = (0u32, 0u32);

    while br.read_line(&mut line)? != 0 {
        if i <= skip_header {
            i += 1;
            line.clear();
            continue;
        }
        if line.as_bytes().get(0) == Some(&b'#') || line.trim().is_empty() {
            line.clear();
            i += 1;
            continue;
        }
        let mut iter = line.split('\t');
        let _gid = iter.next();
        let x: u32 = iter.next().ok_or_else(|| anyhow!("缺 x"))?.trim().parse()?;
        let y: u32 = iter.next().ok_or_else(|| anyhow!("缺 y"))?.trim().parse()?;
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
        line.clear();
        i += 1;
    }
    if max_x == 0 && max_y == 0 {
        bail!("未读取到有效坐标；检查文件或 region");
    }
    Ok(Range {
        min_x,
        max_x,
        min_y,
        max_y,
    })
}

/// 第二遍：读取并聚合到各 bin 的稀疏网格
fn gem_to_sparse_bins(
    path: &str,
    bins: &[u32],
    region: Option<[u32; 4]>,
    has_exon: bool,
    chunksize: usize,
) -> Result<(Range, Vec<HashMap<Key, Agg, RandomState>>)> {
    // 1) 为每个 bin 初始化一个 HashMap；使用 AHash（ahash）作为更快的哈希器
    //    maps[i] 对应 bins[i]，用于存放该 bin 下的稀疏网格
    let mut maps: Vec<HashMap<Key, Agg, RandomState>> = bins
        .iter()
        .map(|_| HashMap::with_hasher(RandomState::new()))
        .collect();

    // 2) 构造一个 TSV 读取器（无表头，制表符分隔）
    //    open_text 支持 .gz 与常规文本
    let rdr = open_text(path)?;
    let mut csv = ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#')) // 忽略以 # 开头的行
        .delimiter(b'\t')
        .from_reader(rdr);

    // 3) 初始化全局范围；max_* 为 0、min_* 为 u32::MAX，以便后续取 min/max
    let mut r = Range {
        min_x: u32::MAX,
        max_x: 0,
        min_y: u32::MAX,
        max_y: 0,
    };
    let mut nrow = 0usize;

    // 4) 主循环：逐条记录读取
    //    `csv.records()` 产生 Result<StringRecord, Error>；这里直接 `?` 抛错
    for rec in csv.records() {
        let rec = rec?;
        // 4.1) 跳过表头：
        if rec.get(0) == Some("geneID") {
            continue;
        }
        // 4.2) 解析所需列：x, y, MIDCount, [ExonCount]
        //      这里假定数据列齐全且合法（不考虑错误分支）
        let x: u32 = rec.get(1).ok_or_else(|| anyhow!("缺 x"))?.parse()?;
        let y: u32 = rec.get(2).ok_or_else(|| anyhow!("缺 y"))?.parse()?;
        let mid: u32 = rec.get(3).ok_or_else(|| anyhow!("缺 MIDCount"))?.parse()?;
        let exon: u32 = if has_exon {
            rec.get(4).unwrap_or("0").parse()?
        } else {
            0
        };

        // 4.3) 可选的矩形裁剪：落在区域外的点直接丢弃
        if let Some([minx, maxx, miny, maxy]) = region {
            if x < minx || x > maxx || y < miny || y > maxy {
                continue;
            }
        }

        // 4.4) 更新全局范围（用于回填二维矩阵时确定边界与大小）
        r.min_x = r.min_x.min(x);
        r.max_x = r.max_x.max(x);
        r.min_y = r.min_y.min(y);
        r.max_y = r.max_y.max(y);

        // 4.5) 对每个 bin 更新对应网格格子的聚合值
        //      key_xy 会把 (x,y) 映射到 (bx,by) 并编码为 64 位 Key
        for (i, &bin) in bins.iter().enumerate() {
            let k = key_xy(x, y, bin);
            let e = maps[i].entry(k).or_default();
            e.mid = e.mid.saturating_add(mid);
            if has_exon {
                e.exon = e.exon.saturating_add(exon);
            }
        }
        // 4.6) 可选：进度日志
        nrow += 1;
        if nrow % (chunksize * 2) == 0 {
            eprintln!(".. processed {} rows", nrow);
        }
    }
    if r.max_x == 0 && r.max_y == 0 {
        bail!("未读取到有效坐标；检查文件或 region");
    }
    Ok((r, maps))
}

/// 写 /dnb/bin{bin}（把稀疏回填成 2D，再一次写入）
fn write_dnb_h5(
    path: &str,
    bin: u32,
    range: &Range,
    map: &HashMap<Key, Agg, RandomState>,
) -> Result<()> {
    let file = if std::path::Path::new(path).exists() {
        H5File::open_rw(path)?
    } else {
        H5File::create(path)?
    };

    let min_x = (range.min_x / bin) * bin;
    let min_y = (range.min_y / bin) * bin;
    let len_x = (range.max_x / bin) - (range.min_x / bin) + 1;
    let len_y = (range.max_y / bin) - (range.min_y / bin) + 1;

    // 回填为一维缓冲（行主序：iy * len_x + ix）
    let mut buf = vec![0u32; (len_x as usize) * (len_y as usize)];
    for (&k, v) in map.iter() {
        let bx = (k >> 32) as u32;
        let by = (k & 0xFFFF_FFFF) as u32;
        let ix = (bx - (min_x / bin)) as usize;
        let iy = (by - (min_y / bin)) as usize;
        buf[iy * (len_x as usize) + ix] = v.mid as u32;
    }

    let g = file.create_group("dnb").or_else(|_| file.group("dnb"))?;
    // 关键：先定义二维 shape，再写一维 buf
    let dset = g
        .new_dataset::<u32>()
        .shape([len_y as usize, len_x as usize])
        .chunk([len_y.min(1024) as usize, len_x.min(1024) as usize])
        .deflate(4)
        .create(format!("bin{bin}").as_str())?;

    // Writes a 1-dimensional array view into a dataset/attribute in memory order.
    // 用.write()会报:shape mismatch when writing: memory = [298286586], destination = [17043, 17502]
    dset.write_raw(&buf)?; // 一维扁平缓冲会按顺序填满 2D 数据集

    dset.new_attr::<i32>()
        .create("min_x")?
        .write_scalar(&(min_x as i32))?;
    dset.new_attr::<i32>()
        .create("min_y")?
        .write_scalar(&(min_y as i32))?;
    dset.new_attr::<i32>()
        .create("len_x")?
        .write_scalar(&(len_x as i32))?;
    dset.new_attr::<i32>()
        .create("len_y")?
        .write_scalar(&(len_y as i32))?;
    dset.new_attr::<i32>()
        .create("bin")?
        .write_scalar(&(bin as i32))?;
    Ok(())
}

fn main() -> Result<()> {
    // 命令行参数
    let args = Args::parse();
    // 读取gem文件头
    let hdr = parse_header(&args.input)?;
    // 输出配准参数和外显子数据
    eprintln!(
        "header: offset=({},{}), has_exon={}",
        hdr.offset_x, hdr.offset_y, hdr.has_exon
    );

    // 可选：粗略范围（只是日志）
    let r0 = first_pass_range(&args.input, hdr.header_line_index)?;
    // 输出横纵坐标大致范围
    eprintln!(
        "coarse range: x=[{},{}], y=[{},{}]",
        r0.min_x, r0.max_x, r0.min_y, r0.max_y
    );

    let region = args.region.as_ref().map(|s| parse_region(s)).transpose()?;
    let (range, maps) = gem_to_sparse_bins(
        &args.input,
        &args.bins,
        region,
        hdr.has_exon,
        args.chunksize,
    )?;

    // 顶层属性（一次写）
    {
        let f = if std::path::Path::new(&args.output).exists() {
            H5File::open_rw(&args.output)?
        } else {
            H5File::create(&args.output)?
        };

        f.new_attr::<i32>()
            .create("resolution")?
            .write_scalar(&args.resolution)?;
        let vstr = unsafe { VarLenUnicode::from_str_unchecked(args.omics.as_str()) };
        f.new_attr::<VarLenUnicode>()
            .create("omics")?
            .write_scalar(&vstr)?;
        f.new_attr::<i32>()
            .create("offset_x")?
            .write_scalar(&hdr.offset_x)?;
        f.new_attr::<i32>()
            .create("offset_y")?
            .write_scalar(&hdr.offset_y)?;
    }

    for (i, &bin) in args.bins.iter().enumerate() {
        eprintln!("writing /dnb/bin{bin} ...");
        write_dnb_h5(&args.output, bin, &range, &maps[i])?;
    }
    eprintln!("done: {}", &args.output);
    Ok(())
}
