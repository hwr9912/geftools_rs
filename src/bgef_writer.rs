use std::collections::{BTreeMap, HashMap};
use std::error::Error;

use hdf5::{types::VarLenUnicode, File as H5File, H5Type};

pub const GEFTOOL_RS_VERSION: u32 = 4;

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
#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type, Default)]
pub struct SpotGene {
    pub mid_count: u32,
    pub gene_count: u16,
}

struct BgefWriter {
    output: String,
    gene_bins: BTreeMap<String, HashMap<(i32, i32), (u32, u32)>>,
    pub h5file: H5File,
    bin: Vec<u16>,
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    max_exp: u32,
    max_exon: u32,
    resolution: u16,
}

impl BgefWriter {
    pub fn new(
        output: String,
        gene_bins: BTreeMap<String, HashMap<(i32, i32), (u32, u32)>>,
        bin: Vec<u16>,
        min_x: i32,
        min_y: i32,
        max_x: i32,
        max_y: i32,
        resolution: u16,
        has_exon: bool,
    ) -> Result<Self, Box<dyn Error>> {
        // 创建 HDF5 文件
        let h5file = H5File::create(&output)?;
        let mut max_exp = u32::MIN;
        let mut max_exon = u32::MIN;

        // 计算总条数（按记录条数，而非坐标数）
        let total: usize = gene_bins.values().map(|coord_map| coord_map.len()).sum();

        // 预分配
        let mut expressions: Vec<Expression> = Vec::with_capacity(total);
        let mut exons: Vec<u32> = if has_exon {
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
                if has_exon {
                    // 按先gene后坐标位置顺序推入外显子计数
                    exons.push(exon_cnt);
                }
                // 逐个比较得到所有binN的所有gene中表达最大值
                if mid > max_exp {
                    max_exp = mid;
                }
                // 同上，这里换成了外显子表达最大值
                if has_exon && exon_cnt > max_exon {
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
        // 返回初始化后的结构体
        Ok(Self {
            output,
            gene_bins,
            h5file,
            bin,
            min_x,
            min_y,
            max_x,
            max_y,
            max_exp,
            max_exon,
            resolution,
        })
    }

    fn write_bgef_gene_exp(self) -> Result<(), Box<dyn std::error::Error>> {
        let gene_exp = self.h5file.create_group("geneExp")?;
        let gene_exp_bin1 = gene_exp.create_group("bin1")?;
        // let ds_expr = gene_exp_bin1
        //     .new_dataset_builder()
        //     .with_data(&expressions)
        //     .create("expression")?;

        // // 写属性
        // ds_expr
        //     .new_attr::<i32>()
        //     .create("minX")?
        //     .write_scalar(&min_x)?;
        // ds_expr
        //     .new_attr::<i32>()
        //     .create("minY")?
        //     .write_scalar(&min_y)?;
        // ds_expr
        //     .new_attr::<i32>()
        //     .create("maxX")?
        //     .write_scalar(&max_x)?;
        // ds_expr
        //     .new_attr::<i32>()
        //     .create("maxY")?
        //     .write_scalar(&max_y)?;
        // ds_expr
        //     .new_attr::<u32>()
        //     .create("maxExp")?
        //     .write_scalar(&max_exp)?;
        // ds_expr
        //     .new_attr::<u32>()
        //     .create("resolution")?
        //     .write_scalar(&args.resolution)?;
        Ok(())
    }
}
