use std::collections::HashMap;
use std::error::Error;

use hdf5::types::VarLenUnicode;
use hdf5::{File as H5File, H5Type}; // 导入 Location trait
use ndarray::Array2;

// 假设 gem_reader 模块提供了 pub fn map2mat 和 pub struct Header
use crate::gem_reader::{map2mat, Header};

pub const GEFTOOL_RS_VERSION: u32 = 4;

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type)]
pub struct Expression {
    pub x: i32,
    pub y: i32,
    pub count: u32, // MIDCount
}

#[repr(C)]
#[derive(H5Type, Clone, Debug)]
#[allow(non_snake_case)]
pub struct GeneRec {
    pub geneID: VarLenUnicode,
    pub geneName: VarLenUnicode,
    pub offset: u32,
    pub count: u32,
}
#[repr(C)]
#[allow(non_snake_case)]
#[derive(Clone, Copy, Debug, H5Type, Default)]
pub struct SpotGene {
    pub MIDcount: u32,
    pub genecount: u16,
}

// BgefWriter 现在持有所有预处理过的数据，准备写入
pub struct BgefWriter {
    output: String,
    expressions: Vec<Expression>,
    genes_meta: Vec<GeneRec>,
    exons: Vec<u32>,
    spot_mid_map: HashMap<(i32, i32), SpotGene>,
    spot_exon_map: HashMap<(i32, i32), u32>,
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    max_exp: u32,
    max_exon: u32,
    resolution: u16,
    has_exon: bool,
}

impl BgefWriter {
    // new 函数现在只接收数据，不进行处理
    #[allow(clippy::too_many_arguments)] // 参数较多是合理的，因为数据都在 main 中处理
    pub fn new(
        output: String,
        expressions: Vec<Expression>,
        genes_meta: Vec<GeneRec>,
        exons: Vec<u32>,
        spot_mid_map: HashMap<(i32, i32), SpotGene>,
        spot_exon_map: HashMap<(i32, i32), u32>,
        min_x: i32,
        min_y: i32,
        max_x: i32,
        max_y: i32,
        max_exp: u32,
        max_exon: u32,
        resolution: u16,
        has_exon: bool,
    ) -> Self {
        Self {
            output,
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
            resolution,
            has_exon,
        }
    }

    /// 将所有数据写入 HDF5 文件
    pub fn write_all(self, hdr: &Header) -> Result<(), Box<dyn Error>> {
        // 1. 创建 HDF5 文件
        let f = H5File::create(&self.output)?;

        // 2. 写入根属性 (使用安全且简洁的方式)
        // bin类型
        let vstr = hdr.bin_type.parse::<VarLenUnicode>()?;
        f.new_attr::<VarLenUnicode>().create("bin_type")?.write_scalar(&vstr)?;
        // 组织区域+工具版本（硬编码）
        f.new_attr::<f32>().create("gef_area")?.write_scalar(&0.0)?;
        f.new_attr::<[u32; 3]>().create("geftool_ver")?.write_scalar(&[1, 1, 20])?;
        // 组学类型
        let vstr = hdr.omics.parse::<VarLenUnicode>()?;
        f.new_attr::<VarLenUnicode>().create("omics")?.write_scalar(&vstr)?;
        // 版本号
        f.new_attr::<u32>().create("version")?.write_scalar(&GEFTOOL_RS_VERSION)?;
        // 芯片代号
        let vstr = hdr.stereo_seq_chip.parse::<VarLenUnicode>()?;
        f.new_attr::<VarLenUnicode>().create("sn")?.write_scalar(&vstr)?;

        // 3. 写入 /geneExp/bin1
        let gene_exp = f.create_group("geneExp")?;
        let gene_exp_bin1 = gene_exp.create_group("bin1")?;

        let ds_expr = gene_exp_bin1
            .new_dataset_builder()
            .with_data(&self.expressions)
            .create("expression")?;

        // 写入 expression 属性
        ds_expr.new_attr::<i32>().create("minX")?.write_scalar(&self.min_x)?;
        ds_expr.new_attr::<i32>().create("minY")?.write_scalar(&self.min_y)?;
        ds_expr.new_attr::<i32>().create("maxX")?.write_scalar(&self.max_x)?;
        ds_expr.new_attr::<i32>().create("maxY")?.write_scalar(&self.max_y)?;
        ds_expr.new_attr::<u32>().create("maxExp")?.write_scalar(&self.max_exp)?;
        ds_expr.new_attr::<u32>().create("resolution")?.write_scalar(&self.resolution)?;

        // 写入 /geneExp/bin1/exon (可选)
        if self.has_exon {
            debug_assert_eq!(self.exons.len(), self.expressions.len());
            let ds_exon =
                gene_exp_bin1.new_dataset_builder().with_data(&self.exons).create("exon")?;
            ds_exon.new_attr::<u32>().create("maxExon")?.write_scalar(&self.max_exon)?;
        }

        // 写入 /geneExp/bin1/gene
        let _ds_gene =
            gene_exp_bin1.new_dataset_builder().with_data(&self.genes_meta).create("gene")?;

        // 4. 写入 /wholeExp/bin1
        let len_x = (self.max_x - self.min_x + 1) as i32;
        let len_y = (self.max_y - self.min_y + 1) as i32;

        let max_mid_per_bin_n = self.spot_mid_map.values().map(|&x| x.MIDcount).max().unwrap_or(0);
        let max_gene = self.spot_mid_map.values().map(|&x| x.genecount).max().unwrap_or(0);
        let number = self.spot_mid_map.len() as u64;

        let whole_exp = f.create_group("wholeExp")?;
        let mat: Vec<SpotGene> = map2mat(&self.spot_mid_map, len_x as usize, len_y as usize)?;
        let whole_exp_bin1 = whole_exp
            .new_dataset::<SpotGene>()
            .shape([len_x as usize, len_y as usize])
            .create("bin1")?;
        whole_exp_bin1
            .write(&Array2::from_shape_vec([len_x as usize, len_y as usize], mat).unwrap())?;

        // 写入 wholeExp 属性
        // 稠密矩阵中非零点的数量
        whole_exp_bin1.new_attr::<u64>().create("number")?.write_scalar(&number)?;
        // 非零点x和y坐标最小值
        whole_exp_bin1.new_attr::<i32>().create("minX")?.write_scalar(&self.min_x)?;
        whole_exp_bin1.new_attr::<i32>().create("minY")?.write_scalar(&self.min_y)?;
        // 非零点x和y坐标极差
        whole_exp_bin1
            .new_attr::<i32>()
            .create("lenX")?
            .write_scalar(&(self.max_x - self.min_x + 1i32))?;
        whole_exp_bin1
            .new_attr::<i32>()
            .create("lenY")?
            .write_scalar(&(self.max_y - self.min_y + 1i32))?;
        // spot中最大的 MID 计数
        whole_exp_bin1.new_attr::<u32>().create("maxMID")?.write_scalar(&max_mid_per_bin_n)?;
        // spot中最大的基因类型计数
        whole_exp_bin1.new_attr::<u32>().create("maxGene")?.write_scalar(&max_gene)?;
        whole_exp_bin1.new_attr::<u32>().create("resolution")?.write_scalar(&self.resolution)?;

        // 5. 写入 /wholeExpExon/bin1
        let whole_exon = f.create_group("wholeExpExon")?;
        let whole_exp_bin1_exon = whole_exon
            .new_dataset::<u32>()
            .shape([len_x as usize, len_y as usize])
            .create("bin1")?;

        let exon_mat: Vec<u32> = map2mat(&self.spot_exon_map, len_x as usize, len_y as usize)?;
        whole_exp_bin1_exon
            .write(&Array2::from_shape_vec([len_x as usize, len_y as usize], exon_mat).unwrap())?;
        // 当分箱大小为 N 时，斑点中的最大外显子表达计数
        let max_exon_spot = self.spot_exon_map.values().max().copied().unwrap_or(0);
        whole_exp_bin1.new_attr::<u32>().create("maxExon")?.write_scalar(&max_exon_spot)?;

        // (stat/gene 逻辑被注释掉了，这里也忽略)

        Ok(())
    }
}
