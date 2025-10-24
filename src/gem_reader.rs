use anyhow::*;
use flate2::read::GzDecoder;
use hdf5::types::VarLenUnicode;
use hdf5::H5Type;
use std::cmp::{max, min};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

/// 处理gz文件
fn open_text(path: &str) -> Result<Box<dyn Read>> {
    let f = File::open(path).with_context(|| format!("open {}", path))?;
    if path.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}

#[derive(Debug, Clone)]
pub struct Header {
    pub bin_type: String,         // #BinType
    pub bin_size: u32,            // #BinSize
    pub omics: String,            // #Omics
    pub stereo_seq_chip: String,  // #Stereo-seqChip
    pub offset_x: i32,            // #OffsetX
    pub offset_y: i32,            // #OffsetY
    pub has_exon: bool,           // gene 表头列数==5 则有 ExonCount
    pub header_line_index: usize, // gene 表头所在的行号（0-based）
}

/// 读取gem文件头输出文件头注释信息
/// 输入：
///     文件地址
/// 返回：
///     bin_type,
///     bin_size,
///     omics,
///     stereo_seq_chip: sn,
///     offset_x: ox,
///     offset_y: oy,
///     has_exon: cols == 5,
///     header_line_index: idx,
pub fn parse_header(path: &str) -> Result<Header> {
    // 读取函数实例加载
    let rdr = open_text(path)?;
    // 读取函数实例加载进缓冲区
    let mut br = BufReader::new(rdr);
    // 新建暂存字符串
    let mut line = String::new();
    // 关键参数
    let mut bin_type = String::new();
    let mut bin_size: u32 = 1;
    let mut omics = String::from("Transcriptomics");
    let mut sn = String::new();
    let (mut ox, mut oy) = (0, 0);
    let mut header_line: Option<(usize, String)> = None;
    // 计数器
    let mut i = 0usize;
    loop {
        line.clear();
        // 从br读取一行，追加到line，然后输出总字节数为n，如果是空行就break
        if br.read_line(&mut line)? == 0 {
            break;
        }
        if let Some(rest) = line.strip_prefix("#BinType=") {
            bin_type = rest.trim().to_owned();
        } else if let Some(rest) = line.strip_prefix("#BinSize=") {
            bin_size = rest.trim().parse()?;
        } else if let Some(rest) = line.strip_prefix("#Omics=") {
            omics = rest.trim().to_owned();
        } else if let Some(rest) = line.strip_prefix("#Stereo-seqChip=") {
            sn = rest.trim().to_owned();
        } else if let Some(_rest) = line.strip_prefix("#OffsetX=") {
            ox = _rest.trim().parse()?;
        } else if let Some(_rest) = line.strip_prefix("#OffsetY=") {
            oy = _rest.trim().parse()?;
        } else if line.starts_with("geneID") {
            header_line = Some((i, line.clone()));
            break;
        }
        i += 1;
    }
    let (idx, header) = header_line.ok_or_else(|| anyhow!("未找到表头 geneID"))?;
    let cols = header.matches('\t').count() + 1;
    Ok(Header {
        bin_type,
        bin_size,
        omics,
        stereo_seq_chip: sn,
        offset_x: ox,
        offset_y: oy,
        has_exon: cols == 5,
        header_line_index: idx,
    })
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
pub struct Expression {
    pub x: i32,
    pub y: i32,
    pub count: u32, // MIDCount
}

// gene 复合类型（注意 repr(C) + VarLenUnicode）
#[repr(C)]
#[derive(H5Type, Clone, Debug)]
pub struct GeneRec {
    pub gene_id: VarLenUnicode, // 比如 "ENSMUSG..."；现阶段可与 gene_name 相同
    pub gene_name: VarLenUnicode, // 比如 "Pdgfrb"；若无映射，先同 gene_id
    pub offset: u32,
    pub count: u32,
}

/// 遍历gem所有表达量行
pub fn get_expression(
    path: &str,
    header_line_index: usize,
    has_exon: bool,
) -> Result<(
    BTreeMap<String, HashMap<(i32, i32), (u32, u32)>>,
    i32,
    i32,
    i32,
    i32,
    u32,
    u32,
)> {
    // 读取函数实例加载
    let rdr = open_text(path)?;
    // 读取函数实例加载进缓冲区
    let mut br = BufReader::new(rdr);
    // 计数器
    let mut i = 0usize;
    // 缓存变量
    //
    let mut gene_bins: BTreeMap<String, HashMap<(i32, i32), (u32, u32)>> = BTreeMap::new();
    let mut min_x = i32::MAX;
    let mut min_y = i32::MAX;
    let mut max_x = i32::MIN;
    let mut max_y = i32::MIN;
    let mut max_exp = u32::MIN;
    let mut max_exon = u32::MIN;

    // 逐行读取正文
    let mut line = String::new();
    loop {
        line.clear();
        // 从br读取一行，追加到line，然后输出总字节数为n，如果是空行就break
        if br.read_line(&mut line)? == 0 {
            break;
        }
        i += 1;
        if i <= header_line_index + 1 {
            continue;
        } else if i % 10000000 == 0 {
            println!("Processed {:>10} lines...", i);
        }
        // 注意去掉行尾 \r\n
        let mut it = line.trim_end_matches(&['\r', '\n'][..]).split('\t');
        // 依次读取：geneID  x  y  MIDCount [exon]
        let gene = match it.next() {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let x: i32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing x"))?
            .parse()?;
        let y: i32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing y"))?
            .parse()?;
        let mid: u32 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing MIDCount"))?
            .parse()?;
        // 更新最小和最大坐标范围
        min_x = min(min_x, x);
        min_y = min(min_y, y);
        max_x = max(max_x, x);
        max_y = max(max_y, y);
        // gene_bins.entry(k)：进入最外层 BTreeMap 的“入口”，只查一次键 k
        // .or_default()：如果这个基因不存在，就插入默认值（HashMap::default()）；存在就直接返回那个值的可变引用
        let inner = gene_bins.entry(gene.to_string()).or_default();
        // 索引特定bin上的gene表达量
        let vals = inner.entry((x, y)).or_insert((0, 0));
        vals.0 += mid; // 合并本行 MID 到该格
        max_exp = max(max_exp, vals.0); // 更新合并后的最大值
        if has_exon {
            // 有exon的情况下新建
            let exon: u32 = it
                .next()
                .ok_or_else(|| anyhow::anyhow!("missing Exon"))?
                .parse()?;
            vals.1 += exon; // 合并本行 exon 到该格
            max_exon = max(max_exon, vals.1); // 更新合并后的最大值
        }
    }
    println!("Processed all lines: {:>10} ", i);
    Ok((gene_bins, min_x, max_x, min_y, max_y, max_exp, max_exon))
}

/// 把坐标形式的hashmap展开成二维稠密矩阵
///
/// Params:
///
///     spot_map: HashMap<(i32, i32), T>, 坐标必须是(i32, i32);
///     len_x, len_y: 稠密矩阵将为[len_x * len_y];
///
/// Returns:
///
///     Result<Vec<Vec<T>>>, 展开后的矩阵
pub fn map2mat<T>(spot_map: &HashMap<(i32, i32), T>, len_x: usize, len_y: usize) -> Result<Vec<T>>
where
    T: Default + Clone,
{
    // 初始化一个填满默认值的二维矩阵
    let mut mat: Vec<Vec<T>> = vec![vec![T::default(); len_x]; len_y];

    // 填充矩阵
    for (&(x, y), v) in spot_map.iter() {
        if x >= 0 && y >= 0 {
            let (xi, yi) = (x as usize, y as usize);
            if xi < len_x && yi < len_y {
                mat[yi][xi] = v.clone();
            }
        }
    }

    // 输出为
    let vector = mat.into_iter().flat_map(|v| v).collect();

    Ok(vector)
}
