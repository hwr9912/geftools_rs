use anyhow::*;
use flate2::read::GzDecoder;
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
        let x: i32 = it.next().ok_or_else(|| anyhow::anyhow!("missing x"))?.parse()?;
        let y: i32 = it.next().ok_or_else(|| anyhow::anyhow!("missing y"))?.parse()?;
        let mid: u32 = it.next().ok_or_else(|| anyhow::anyhow!("missing MIDCount"))?.parse()?;
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
            let exon: u32 = it.next().ok_or_else(|| anyhow::anyhow!("missing Exon"))?.parse()?;
            vals.1 += exon; // 合并本行 exon 到该格
            max_exon = max(max_exon, vals.1); // 更新合并后的最大值
        }
    }
    println!("Processed all lines: {:>10} ", i);
    Ok((gene_bins, min_x, max_x, min_y, max_y, max_exp, max_exon))
}

/// 将 HashMap 坐标数据截取在指定范围 (min_x, max_x) × (min_y, max_y) 内，
/// 并转换为局部矩阵坐标（从 0 开始）。
///
/// # Params
///
/// - `spot_map`: `HashMap<(i32, i32), T>`  
///   原始坐标与值的映射，键为 (x, y)，必须为 `i32` 类型；  
///
/// - `min_x`, `min_y`: `i32`  
///   截取矩形区域的左上角坐标；  
///
/// - `max_x`, `max_y`: `i32`  
///   截取矩形区域的右下角坐标；  
///
/// - `len_x`, `len_y`: `usize`  
///   输出矩阵的宽度与高度（可由 `max_x - min_x + 1`、`max_y - min_y + 1` 计算得到）；  
///
/// # Returns
///
/// - `Result<Vec<Vec<T>>>`  
///   返回截取并重置坐标后的稠密矩阵，
///   其中矩阵尺寸为 `[len_y][len_x]`，空缺位置以 `T::default()` 填充。
pub fn map2mat<T>(
    spot_map: &HashMap<(i32, i32), T>,
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
) -> Result<(Vec<T>, usize, usize)>
where
    T: Default + Clone,
{
    // 计算宽高（注意这里 +1 是为了包含边界）
    let width = (max_x - min_x + 1) as usize;
    let height = (max_y - min_y + 1) as usize;

    // 初始化默认值矩阵
    let mut mat: Vec<Vec<T>> = vec![vec![T::default(); width]; height];

    // 遍历所有点，筛选落在指定范围内的
    for (&(x, y), value) in spot_map.iter() {
        if x >= min_x && x <= max_x && y >= min_y && y <= max_y {
            // 转换为矩阵坐标（从 0 开始）
            let xi = (x - min_x) as usize;
            let yi = (y - min_y) as usize;
            mat[yi][xi] = value.clone();
        }
    }

    // 将二维矩阵展开为一维向量
    let vector = mat.into_iter().flat_map(|row| row).collect();

    Ok((vector, width, height))
}

