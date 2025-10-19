use anyhow::Result;
use clap::Parser;
use hdf5::{File as H5File, H5Type, types::VarLenUnicode};
use ndarray::Array2;

mod bgef_writer;
mod gem_reader;

use bgef_writer::GEFTOOL_RS_VERSION;
use gem_reader::parse_header;

#[derive(Parser)]
struct Args {
    /// 输入 GEM 或 GEM.GZ
    #[arg(short, long, default_value = "test10000.gem.gz")]
    input: String,
    /// 输出 bGEF (HDF5)
    #[arg(short, long, default_value = "dummy.bgef")]
    output: String,
    /// 逗号分隔的 bin 列表
    #[arg(short, long, value_delimiter = ',', default_value = "1,20,50,100")]
    bins: Vec<u32>,
    // /// 可选区域：minx,maxx,miny,maxy
    // #[arg(long)]
    // region: Option<String>,
    // /// 分块大小（行）
    // #[arg(long, default_value_t = 2_000_000)]
    // chunksize: usize,
    // /// 顶层属性：resolution
    // #[arg(long, default_value_t = 1)]
    // resolution: i32,
    // /// 顶层属性：omics
    // #[arg(long, default_value = "Transcriptomics")]
    // omics: String,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct Expr {
    x: i32,
    y: i32,
    count: u32, // MIDCount
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct GeneRec {
    gene_id: u32,   // 简化：用数字代替字符串
    gene_name: u32, // 简化：用数字代替字符串
    offset: i32,    // 在 expression 数组中的起始下标
    count: i32,     // 该基因对应多少条 expression 记录
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct StatGene {
    gene_id: u32,
    gene_name: u32,
    mid_count: u32,
    log10_mid: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct Pixel {
    mid: u32,
    exon: u32,
}

fn main() -> Result<()> {
    // 命令行参数
    let args = Args::parse();
    // 读取 gem 文件头
    let hdr = parse_header(&args.input)?;
    // 输出配准参数和外显子数据
    eprintln!(
        "Header info:\n  BinType={}  BinSize={}\n  Omics={}  Chip={}\n  Offset=({}, {})  HasExon={}  HeaderLineIndex={}",
        hdr.bin_type,
        hdr.bin_size,
        hdr.omics,
        hdr.stereo_seq_chip,
        hdr.offset_x,
        hdr.offset_y,
        hdr.has_exon,
        hdr.header_line_index
    );

    // // 可选：粗略范围（只是日志）
    // let r0 = first_pass_range(&args.input, hdr.header_line_index)?;
    // // 输出横纵坐标大致范围
    // eprintln!(
    //     "coarse range: x=[{},{}], y=[{},{}]",
    //     r0.min_x, r0.max_x, r0.min_y, r0.max_y
    // );

    // let region = args.region.as_ref().map(|s| parse_region(s)).transpose()?;
    // let (range, maps) = gem_to_sparse_bins(
    //     &args.input,
    //     &args.bins,
    //     region,
    //     hdr.has_exon,
    //     args.chunksize,
    // )?;

    // 顶层属性（一次写）
    let f = H5File::create(&args.output)?;

    let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.bin_type.as_str()) };
    f.new_attr::<VarLenUnicode>() // 一维数组属性
        .create("bin_type")?
        .write_scalar(&vstr)?;
    f.new_attr::<f32>().create("gef_area")?.write_scalar(&0)?;
    f.new_attr::<[u32; 3]>()
        .create("geftool_ver")?
        .write_scalar(&[1, 1, 20])?;
    let bytes = hdr.omics.as_bytes();
    f.new_attr::<u8>()
        .shape([bytes.len()])
        .create("omics")?
        .write(bytes)?;
    f.new_attr::<u32>()
        .create("version")?
        .write_scalar(&GEFTOOL_RS_VERSION)?;
    let bytes = hdr.stereo_seq_chip.as_bytes();
    f.new_attr::<u8>()
        .shape([bytes.len()])
        .create("sn")?
        .write(bytes)?;

    // ---------- geneExp/bin1 ----------
    let g_geneexp = f.create_group("geneExp")?;
    let g_bin1 = g_geneexp.create_group("bin1")?;

    // 随便构造几条 expression，假设 5 条
    let expressions: Vec<Expr> = vec![
        Expr {
            x: 10,
            y: 20,
            count: 1,
        },
        Expr {
            x: 11,
            y: 21,
            count: 2,
        },
        Expr {
            x: 12,
            y: 22,
            count: 1,
        },
        Expr {
            x: 30,
            y: 40,
            count: 3,
        },
        Expr {
            x: 31,
            y: 41,
            count: 1,
        },
    ];
    let ds_expr = g_bin1
        .new_dataset_builder()
        .with_data(&expressions)
        .create("expression")?;

    // exon 标记（与 expression 一一对应）
    let exon_flags: Vec<u8> = vec![1, 0, 1, 1, 0];
    let _ds_exon = g_bin1
        .new_dataset_builder()
        .with_data(&exon_flags)
        .create("exon")?;

    // gene 索引（简化：3 个基因，映射到 expression 的分段）
    // gene0 -> 0..2，gene1 -> 2..4，gene2 -> 4..5
    let genes: Vec<GeneRec> = vec![
        GeneRec {
            gene_id: 1,
            gene_name: 1001,
            offset: 0,
            count: 2,
        },
        GeneRec {
            gene_id: 2,
            gene_name: 1002,
            offset: 2,
            count: 2,
        },
        GeneRec {
            gene_id: 3,
            gene_name: 1003,
            offset: 4,
            count: 1,
        },
    ];
    let _ds_gene = g_bin1
        .new_dataset_builder()
        .with_data(&genes)
        .create("gene")?;

    // ---------- stat/gene ----------
    let g_stat = f.create_group("stat")?;
    let stat_genes: Vec<StatGene> = vec![
        StatGene {
            gene_id: 1,
            gene_name: 1001,
            mid_count: 3,
            log10_mid: (3.0f32).log10(),
        },
        StatGene {
            gene_id: 2,
            gene_name: 1002,
            mid_count: 4,
            log10_mid: (4.0f32).log10(),
        },
        StatGene {
            gene_id: 3,
            gene_name: 1003,
            mid_count: 1,
            log10_mid: (1.0f32).log10(),
        },
    ];
    let _ds_stat_gene = g_stat
        .new_dataset_builder()
        .with_data(&stat_genes)
        .create("gene")?;

    // ---------- wholeExp/bin1 ----------
    let g_whole = f.create_group("wholeExp")?;

    // 做一个 8x6 的网格，随便填
    let (h, w) = (8, 6);
    let whole: Array2<Pixel> = Array2::from_shape_fn((h, w), |(y, x)| Pixel {
        mid: (x as u32 + y as u32) * 10,
        exon: ((x as u32 + y as u32) % 3),
    });
    let _ds_whole = g_whole
        .new_dataset_builder()
        .with_data(&whole)
        .create("bin1")?; // 数据集名和组名同名在 HDF5 是允许的；若你不喜欢可改 "data"

    // ---------- wholeExpExon/bin1 ----------
    let g_exon = f.create_group("wholeExpExon")?;
    let whole_exon: Array2<u32> =
        Array2::from_shape_fn((h, w), |(y, x)| if (x + y) % 2 == 0 { 1 } else { 0 });
    let _ds_exon_grid = g_exon
        .new_dataset_builder()
        .with_data(&whole_exon)
        .create("bin1")?;

    // 确保 expression 被实际使用，避免警告
    let _ = ds_expr.size();

    println!("wrote {}!", &args.output);
    Ok(())
}
