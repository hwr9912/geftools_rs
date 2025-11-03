# 华大 cpp 库 geftools 的 rust 实现

geftools 对于来源于 stereo-seq 某些 gem 数据，可能会出现一些溢出错误

```bash
(st) PS D:\Rust\geftools_rs> & E:/Anaconda_envs/envs/st/python.exe d:/Rust/geftools_rs/gem2gef_bgi_bug.py 
geneID  x       y       MIDCount        ExonCount
 5
create bgef file: out/Y00855N1_bgi.bgef
can not find sn info from input file out/Y00855N1.tissue.gem.gz
bin 1 matrix: min_x=0 len_x=-2147483648 min_y=0 len_y=-2147483648 matrix_len=0
```

本代码仓库使用 rust 重写了 geftools 的部分内容。开发过程中主要参考了：

> 1. [geftools 原始代码仓库](https://github.com/STOmics/geftools)
> 2. [stereopy 使用文档](https://stereopy.readthedocs.io/en/latest/index.html)
> 3. [gefpy: geftools 的 python 封装](https://github.com/STOmics/gefpy)
> 4. [SAW 8.0 用户手册](https://www.stomics.tech/service/saw_8_1/docs/gao-ji-she-zhi/expression-matrix-format.html)
> 5. [官方 GEF 文件结构思维导图](https://www.processon.com/view/link/610cc49c7d9c087bbd1ab7ab#map)

目前已完成的功能：

- gem 转换为 bgef

即将开发的功能：

- 并行化
- cgef 转换

使用说明

```bash
cargo build --release
```

```
(base) william_han@192 geftools_rs % target/release/gem2gef -h
Usage: gem2gef [OPTIONS]

Options:
  -i, --input <INPUT>            输入 GEM 或 GEM.GZ [default: test10000.gem.gz]
  -o, --output <OUTPUT>          输出 bGEF (HDF5) [default: dummy.bgef]
  -b, --bins <BINS>              逗号分隔的 bin 列表 [default: 1,20,50,100]
      --resolution <RESOLUTION>  顶层属性：resolution [default: 500]
  -h, --help                     Print help
```

