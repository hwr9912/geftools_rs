use std::collections::HashMap;

use anyhow::Result;
use clap::Parser;
use hdf5::{types::VarLenUnicode, File as H5File, H5Type};
use ndarray::Array2;

mod bgef_writer;
mod gem_reader;
mod log;

use bgef_writer::{Expression, GeneRec, GEFTOOL_RS_VERSION};
use gem_reader::{map2mat, parse_header};
use log::log_msg;

use crate::{bgef_writer::SpotGene, gem_reader::get_expression};

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
    /// 顶层属性：resolution
    #[arg(long, default_value_t = 500)]
    resolution: u16,
    // /// 顶层属性：omics
    // #[arg(long, default_value = "Transcriptomics")]
    // omics: String,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
struct StatGene {
    gene_id: u32,
    gene_name: u32,
    mid_count: u32,
    log10_mid: f32,
}

fn main() -> Result<()> {
    // 命令行参数
    let args = Args::parse();
    // 读取 gem 文件头
    let hdr = parse_header(&args.input)?;
    // 输出配准参数和外显子数据
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
    let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.omics.as_str()) };
    f.new_attr::<VarLenUnicode>() // 一维数组属性
        .create("omics")?
        .write_scalar(&vstr)?;
    f.new_attr::<u32>()
        .create("version")?
        .write_scalar(&GEFTOOL_RS_VERSION)?;
    let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.stereo_seq_chip.as_str()) };
    f.new_attr::<VarLenUnicode>() // 一维数组属性
        .create("sn")?
        .write_scalar(&vstr)?;

    // ---------- geneExp/bin1 ----------
    // 关键：输出按基因和坐标排布的表达量
    let (gene_bins, min_x, max_x, min_y, max_y, mut max_exp, mut max_exon) =
        get_expression(&args.input, hdr.header_line_index, hdr.has_exon)
            .expect("get_expression 输入有问题");

    // 计算总条数（按记录条数，而非坐标数）
    let total: usize = gene_bins.values().map(|coord_map| coord_map.len()).sum();

    // 预分配
    let mut expressions: Vec<Expression> = Vec::with_capacity(total);
    let mut exons: Vec<u32> = if hdr.has_exon {
        Vec::with_capacity(total)
    } else {
        Vec::new()
    };

    let mut genes_meta: Vec<GeneRec> = Vec::with_capacity(gene_bins.len());

    // 串接顺序：外层按基因（BTreeMap 已排序），内层按 (x,y) 排序
    let mut offset_u32: u32 = 0;
    // 后面还有用gene_bins，这里只能用引用
    for (gene_key, coord_map) in &gene_bins {
        // 1) 收集并排序该基因的所有 (x,y) 记录，保证稳定性
        let mut recs: Vec<_> = coord_map.into_iter().collect(); // Vec<((x,y), (mid,exon))>
        recs.sort_by_key(|&((x, y), _)| (x, y));

        // 2) 记录起始 offset
        let start = offset_u32;

        // 3) 逐条推入 expression / exon，并维护 max 值，为了赋值，必须预先解引用
        for (&(x, y), &(mid, exon_cnt)) in recs {
            expressions.push(Expression { x, y, count: mid });
            if hdr.has_exon {
                // 按先gene后坐标位置顺序推入外显子计数
                exons.push(exon_cnt);
            }
            // 逐个比较得到所有binN的所有gene中表达最大值
            if mid > max_exp {
                max_exp = mid;
            }
            // 同上，这里换成了外显子表达最大值
            if hdr.has_exon && exon_cnt > max_exon {
                max_exon = exon_cnt;
            }
            offset_u32 = offset_u32.saturating_add(1);
        }

        // 4) 构建 gene 行（安全 from_str；如需映射，替换 gene_name_here 即可）
        // 这里先让 gene_id = gene_key，gene_name 同 gene_id；后续你有映射表时把 gene_name 换成符号
        let gene_id_v = unsafe { VarLenUnicode::from_str_unchecked(gene_key.as_str()) };
        let gene_name_v = gene_id_v.clone(); // 现阶段同值；有符号映射时改这里

        let cnt = offset_u32 - start;
        genes_meta.push(GeneRec {
            gene_id: gene_id_v,
            gene_name: gene_name_v,
            offset: start,
            count: cnt,
        });
    }

    // 5) 写 /geneExp/binN/expression
    let gene_exp = f.create_group("geneExp")?;
    let gene_exp_bin1 = gene_exp.create_group("bin1")?;
    let ds_expr = gene_exp_bin1
        .new_dataset_builder()
        .with_data(&expressions)
        .create("expression")?;

    // 写属性
    ds_expr
        .new_attr::<i32>()
        .create("minX")?
        .write_scalar(&min_x)?;
    ds_expr
        .new_attr::<i32>()
        .create("minY")?
        .write_scalar(&min_y)?;
    ds_expr
        .new_attr::<i32>()
        .create("maxX")?
        .write_scalar(&max_x)?;
    ds_expr
        .new_attr::<i32>()
        .create("maxY")?
        .write_scalar(&max_y)?;
    ds_expr
        .new_attr::<u32>()
        .create("maxExp")?
        .write_scalar(&max_exp)?;
    ds_expr
        .new_attr::<u32>()
        .create("resolution")?
        .write_scalar(&args.resolution)?;

    // 6) 写 /geneExp/binN/exon（可选）
    if hdr.has_exon {
        debug_assert_eq!(exons.len(), expressions.len());
        let ds_exon = gene_exp_bin1
            .new_dataset_builder()
            .with_data(&exons)
            .create("exon")?;
        ds_exon
            .new_attr::<u32>()
            .create("maxExon")?
            .write_scalar(&max_exon)?;
    }

    // 7) 写 /geneExp/binN/gene（复合类型一把写）
    let _ds_gene = gene_exp_bin1
        .new_dataset_builder()
        .with_data(&genes_meta)
        .create("gene")?;

    log_msg(&format!(
        "/geneExp/bin1 info:\n  minX={}  maxX={}\n  minY={}  maxY={}\n  maxExp={} maxExon={} resolution={}",
        min_x, max_x,
        min_y, max_y,
        max_exp, max_exon, args.resolution,
        ));

    // ---------- wholeExp/bin1 & wholeExpExon/bin1 ----------
    {
        // wholeExp/bin1 属性 lenX, lenY, maxGene
        let len_x = (max_x - min_x + 1) as i32;
        let len_y = (max_y - min_y + 1) as i32;

        // 总非0坐标点
        let mut unique_coords = std::collections::HashSet::new();
        let mut spot_mid: HashMap<(i32, i32), SpotGene> = HashMap::new();
        let mut spot_exon: HashMap<(i32, i32), u32> = HashMap::new();
        let mut max_mid_per_bin_n = u32::MIN;

        // 使用
        for (_gene, coord_map) in &gene_bins {
            for (&(x, y), &(mid, exon)) in coord_map {
                // 总非0坐标点计数
                unique_coords.insert((x, y));
                // MID 累加 + 基因计数（注意拿可变引用来改）
                let entry = spot_mid.entry((x, y)).or_insert(SpotGene::default());
                entry.mid_count += mid;
                entry.gene_count += 1;
                // 逐个比较得到所有binN的所有gene中表达最大值
                if entry.mid_count > max_mid_per_bin_n {
                    max_mid_per_bin_n = entry.mid_count;
                }
                // wholeExpExon/bin1: exon 累加
                *spot_exon.entry((x, y)).or_insert(0) = spot_exon
                    .get(&(x, y))
                    .copied()
                    .unwrap_or(0)
                    .saturating_add(exon);
            }
        }

        // ###### 写 wholeExp/bin1 属性到 HDF5 ######
        // 每个坐标下基因种类数
        let max_gene = spot_mid.values().map(|&x| x.gene_count).max().unwrap_or(0);
        let number = unique_coords.len() as u64;

        // 创建组
        let whole_exp = f.create_group("wholeExp")?;
        let mat: Vec<SpotGene> = map2mat(&spot_mid, len_x as usize, len_y as usize).unwrap(); // 先解包 Result，mat 是 Vec<T>
        let whole_exp_bin1 = whole_exp
            .new_dataset::<SpotGene>()
            .shape([len_x as usize, len_y as usize])
            .create("bin1")?;
        whole_exp_bin1
            .write(&Array2::from_shape_vec([len_x as usize, len_y as usize], mat).unwrap())?;
        // 写属性
        whole_exp_bin1
            .new_attr::<u64>()
            .create("number")?
            .write_scalar(&number)?;
        whole_exp_bin1
            .new_attr::<i32>()
            .create("minX")?
            .write_scalar(&min_x)?;
        whole_exp_bin1
            .new_attr::<i32>()
            .create("lenX")?
            .write_scalar(&len_x)?;
        whole_exp_bin1
            .new_attr::<i32>()
            .create("minY")?
            .write_scalar(&min_y)?;
        whole_exp_bin1
            .new_attr::<i32>()
            .create("lenY")?
            .write_scalar(&len_y)?;
        whole_exp_bin1
            .new_attr::<u32>()
            .create("maxMID")?
            .write_scalar(&max_mid_per_bin_n)?;
        whole_exp_bin1
            .new_attr::<u32>()
            .create("maxGene")?
            .write_scalar(&max_gene)?;
        whole_exp_bin1
            .new_attr::<u32>()
            .create("resolution")?
            .write_scalar(&args.resolution)?;

        log_msg(&format!(
        "/wholeExp/bin1 info:\n  number={}\n  minX={}  lenX={}\n  minY={}  lenY={}\n  maxMID={}  maxGene={} resolution={}",
        number,
        min_x, len_x,
        min_y, len_y,
        max_mid_per_bin_n, max_gene, args.resolution,
        ));

        // ###### 写 wholeExpExon/bin1 属性到 HDF5 ######
        let whole_exon = f.create_group("wholeExpExon")?;
        let whole_exp_bin1 = whole_exon
            .new_dataset::<u32>()
            .shape([len_x as usize, len_y as usize])
            .create("bin1")?;
        // 计算和输出spot上的最大exon表达量
        let max_exon = spot_exon.values().map(|&x| x).max().unwrap();
        whole_exp_bin1
            .new_attr::<u32>()
            .create("maxExon")?
            .write_scalar(&max_exon)?;

        log_msg(&format!("/wholeExpExon/bin1 info:\n  maxExon={}", max_exon,));
    }

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

    // // 确保 expression 被实际使用，避免警告
    // let _ = ds_expr.size();

    println!("wrote {}!", &args.output);
    Ok(())
}
