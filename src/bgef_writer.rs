use std::collections::{BTreeMap, HashMap};

use hdf5::{types::VarLenUnicode, File as H5File, H5Type};

pub const GEFTOOL_RS_VERSION: u32 = 4;

#[repr(C)]
#[derive(Clone, Copy, Debug, H5Type, Default)]
pub struct SpotGene {
    pub mid_count: u32,
    pub gene_count: u16,
}

struct BgefWriter {
    output: String,
    gene_bins: BTreeMap<String, HashMap<(i32, i32), (u32, u32)>>,
    h5file: H5File,
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    max_exp: u32,
    max_exon: u32,
}

impl BgefWriter {
    fn BgefWriter(
        mut self,
        output: String,
        bin: Vec<u16>,
        resolution: u16,
    ) -> Result<(), Box<dyn std::error::Error>> {
        self.output = output;
        let f = H5File::create(&self.output)?;
        // let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.bin_type.as_str()) };
        // f.new_attr::<VarLenUnicode>() // 一维数组属性
        //     .create("bin_type")?
        //     .write_scalar(&vstr)?;
        // f.new_attr::<f32>().create("gef_area")?.write_scalar(&0)?;
        // f.new_attr::<[u32; 3]>()
        //     .create("geftool_ver")?
        //     .write_scalar(&[1, 1, 20])?;
        // let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.omics.as_str()) };
        // f.new_attr::<VarLenUnicode>() // 一维数组属性
        //     .create("omics")?
        //     .write_scalar(&vstr)?;
        // f.new_attr::<u32>()
        //     .create("version")?
        //     .write_scalar(&GEFTOOL_RS_VERSION)?;
        // let vstr = unsafe { VarLenUnicode::from_str_unchecked(hdr.stereo_seq_chip.as_str()) };
        // f.new_attr::<VarLenUnicode>() // 一维数组属性
        //     .create("sn")?
        //     .write_scalar(&vstr)?;
        Ok(())
    }

    fn write_bgef_whole_exp() {}
}
