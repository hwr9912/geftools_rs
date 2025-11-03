use std::{collections::HashMap, error::Error};

use anyhow::Result;
use clap::Parser;
use hdf5::H5Type;

mod bgef_writer;
mod gem_reader;
mod log;

use crate::{
    bgef_writer::{str2fa64, BgefWriter, Expression, GeneRec, SpotGene},
    gem_reader::{get_expression, parse_header},
    log::log_msg,
};

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
    /// 顶层属性：resolution
    #[arg(long, default_value_t = 500)]
    resolution: u16,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct StatGene {
    gene_id: u32,
    gene_name: u32,
    mid_count: u32,
    log10_mid: f32,
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1. 命令行参数
    let args = Args::parse();
    // 2. 读取 gem 文件头
    let hdr = parse_header(&args.input)?;
    log_msg(&format!(
    "Header info:\n  BinType={}  BinSize={}\n  Omics={}  Chip={}\n  Offset=({}, {})  HasExon={}  HeaderLineIndex={}",
    hdr.bin_type,
    hdr.bin_size,
    hdr.omics,
    hdr.stereo_seq_chip,
    hdr.offset_x,
    hdr.offset_y,
    hdr.has_exon,
    hdr.header_line_index,
    ));

    // 读取和处理 geneExp 数据
    let (gene_bins, min_x, max_x, min_y, max_y, mut max_exp, mut max_exon) =
        get_expression(&args.input, hdr.header_line_index, hdr.has_exon)
            .expect("get_expression 输入有问题");

    // 计算总条数
    let total: usize = gene_bins.values().map(|coord_map| coord_map.len()).sum();

    // 预分配 (为 geneExp/bin1 准备)
    let mut expressions: Vec<Expression> = Vec::with_capacity(total);
    let mut exons: Vec<u32> = if hdr.has_exon {
        Vec::with_capacity(total)
    } else {
        Vec::new()
    };
    let mut genes_meta: Vec<GeneRec> = Vec::with_capacity(gene_bins.len());

    // 预分配 (为 wholeExp 和 wholeExpExon 准备)
    let mut spot_mid_map: HashMap<(i32, i32), SpotGene> = HashMap::new();
    let mut spot_exon_map: HashMap<(i32, i32), u32> = HashMap::new();
    let mut max_mid_per_bin_n = u32::MIN;

    // 串接顺序：外层按基因（BTreeMap 已排序），内层按 (x,y) 排序
    // 表达数据计数器
    let mut offset_u32: u32 = 0;
    for (gene_key, coord_map) in &gene_bins {
        // 1) 收集并排序该基因的所有 (x,y) 记录，保证稳定性
        let mut recs: Vec<_> = coord_map.into_iter().collect(); // Vec<((x,y), (mid,exon))>
        recs.sort_by_key(|&((x, y), _)| (x, y));

        // 2) 记录起始 offset: gene_key基因的数据从哪里开始
        let start = offset_u32;

        // 3) 逐条推入 expression / exon，并维护 max 值
        for (&(x, y), &(mid, exon_cnt)) in recs {
            // --- (原 Loop 1 的逻辑) ---
            expressions.push(Expression { x, y, count: mid });
            if hdr.has_exon {
                exons.push(exon_cnt);
            }
            if mid > max_exp {
                max_exp = mid;
            }
            if hdr.has_exon && exon_cnt > max_exon {
                max_exon = exon_cnt;
            }
            offset_u32 = offset_u32.saturating_add(1);

            // --- (原 Loop 2 的逻辑 - 合并于此) ---
            let entry = spot_mid_map.entry((x, y)).or_insert(SpotGene::default());
            entry.MIDcount += mid;
            entry.genecount += 1;
            if entry.MIDcount > max_mid_per_bin_n {
                max_mid_per_bin_n = entry.MIDcount;
            }
            // exon 累加 (注意变量名是 exon_cnt)
            *spot_exon_map.entry((x, y)).or_insert(0) =
                spot_exon_map.get(&(x, y)).copied().unwrap_or(0).saturating_add(exon_cnt);
        }

        // 4) 构建 gene 行
        let gene_id = str2fa64(gene_key); // 这里需要补充转化为geneid的代码
        let gene_name = str2fa64(gene_key);

        genes_meta.push(GeneRec {
            geneID: gene_id,
            geneName: gene_name,
            offset: start,
            count: offset_u32 - start,
        });
    }

    log_msg(&format!(
        "/geneExp/bin1 info:\n  minX={}  maxX={}\n  minY={}  maxY={}\n  maxExp={} maxExon={} resolution={}",
        min_x, max_x,
        min_y, max_y,
        max_exp, max_exon, args.resolution,
        ));

    // 4. 处理 wholeExp 和 wholeExpExon 数据
    let len_x = (max_x - min_x + 1) as i32;
    let len_y = (max_y - min_y + 1) as i32;

    // 计算统计数据
    let max_gene = spot_mid_map.values().map(|&x| x.genecount).max().unwrap_or(0);
    let number = spot_mid_map.len() as u64;
    let max_exon_spot = spot_exon_map.values().max().copied().unwrap_or(0);

    log_msg(&format!(
        "/wholeExp/bin1 info:\n  number={}\n  minX={}  lenX={}\n  minY={}  lenY={}\n  maxMID={}  maxGene={} resolution={}",
        number,
        min_x, len_x,
        min_y, len_y,
        max_mid_per_bin_n, max_gene, args.resolution,
        ));
    log_msg(&format!(
        "/wholeExpExon/bin1 info:\n  maxExon={}",
        max_exon_spot
    ));

    // 5. 实例化 BgefWriter 并传入所有数据
    let writer = BgefWriter::new(
        args.output.clone(),
        expressions,
        genes_meta,
        exons,
        spot_mid_map,
        spot_exon_map,
        min_x,
        min_y,
        max_x,
        max_y,
        max_exp,
        max_exon,
        args.resolution,
        hdr.has_exon,
    );

    // 6. 执行写入
    writer.write_all(&hdr)?;

    println!("wrote {}!", &args.output);
    Ok(())
}
