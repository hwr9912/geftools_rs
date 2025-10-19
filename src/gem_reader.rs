use anyhow::*;
use flate2::read::GzDecoder;
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
