use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Serialize, Deserialize)]
struct GeneMap {
    /// 小写基因名 -> ENSMUSG...
    by_symbol: HashMap<String, String>,
}

/// 读取 Ensembl 小鼠 GTF（或 .gtf.gz），抽取 gene 级注释，更新/重建本地 JSON 表。
/// 返回写入的条目数。
pub fn update_table_from_gtf(gtf_path: &str, db_json_path: &str) -> Result<usize> {
    // 允许覆盖写入：不存在则新建，存在则重建（简单粗暴，避免脏并发）
    let mut by_symbol: HashMap<String, String> = HashMap::new();

    // GTF 第9列 attributes 的简单解析： gene_id "ENSMUSG..."; gene_name "P2ry12";
    let re_gene_id = Regex::new(r#"gene_id\s+"([^"]+)""#).unwrap();
    let re_gene_name = Regex::new(r#"gene_name\s+"([^"]+)""#).unwrap();

    let reader: Box<dyn Read> = if gtf_path.ends_with(".gz") {
        Box::new(GzDecoder::new(
            File::open(gtf_path).with_context(|| gtf_path.to_string())?,
        ))
    } else {
        Box::new(File::open(gtf_path).with_context(|| gtf_path.to_string())?)
    };

    let buf = BufReader::new(reader);
    for line in buf.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // GTF: chr, source, feature, start, end, score, strand, frame, attributes
        // 只看 feature == "gene"
        let mut it = line.split('\t');
        let _chr = match it.next() {
            Some(v) => v,
            None => continue,
        };
        let _src = match it.next() {
            Some(v) => v,
            None => continue,
        };
        let feature = match it.next() {
            Some(v) => v,
            None => continue,
        };
        if feature != "gene" {
            continue;
        }
        // 跳过 start,end,score,strand,frame
        for _ in 0..5 {
            it.next();
        }
        let attr = match it.next() {
            Some(v) => v,
            None => continue,
        };

        let gid = re_gene_id
            .captures(attr)
            .and_then(|c| c.get(1).map(|m| m.as_str().to_string()));
        let gname = re_gene_name
            .captures(attr)
            .and_then(|c| c.get(1).map(|m| m.as_str().to_string()));

        let (gid, gname) = match (gid, gname) {
            (Some(gid), Some(gname)) => (gid, gname),
            _ => continue,
        };

        // 只收 gene_id 形如 ENSMUSG...
        if !gid.starts_with("ENSMUSG") {
            continue;
        }

        // 基因名大小写不敏感；若重复，保持首次写入（避免把主名被别名覆盖）
        let key = gname.to_lowercase();
        by_symbol.entry(key).or_insert(gid);
    }

    let map = GeneMap { by_symbol };
    let json = serde_json::to_string_pretty(&map)?;
    fs::write(db_json_path, json).with_context(|| db_json_path.to_string())?;
    Ok(map.by_symbol.len())
}

/// 查表：给定 symbol（如 "P2ry12"），返回 Some("ENSMUSG...") 或 None。
pub fn query_gene_id(db_json_path: &str, symbol: &str) -> Result<Option<String>> {
    if !Path::new(db_json_path).exists() {
        // 没有表直接返回 None（或你也可以改成 Err 提示先更新）
        return Ok(None);
    }
    let f = File::open(db_json_path).with_context(|| db_json_path.to_string())?;
    let mut s = String::new();
    BufReader::new(f).read_to_string(&mut s)?;
    let map: GeneMap = serde_json::from_str(&s)?;
    Ok(map.by_symbol.get(&symbol.to_lowercase()).cloned())
}
