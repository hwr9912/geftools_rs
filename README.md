# 华大cpp库geftools的rust实现

geftools对于来源于stereo-seq某些gem数据，可能会出现一些溢出错误，因此本代码仓库使用rust重写了geftools的部分内容。开发过程中主要参考了：

> 1. [geftools原始代码仓库](https://github.com/STOmics/geftools)
> 2. [stereopy使用文档](https://stereopy.readthedocs.io/en/latest/index.html)
> 3. [gefpy: geftools的python封装](https://github.com/STOmics/gefpy)
> 4. [SAW 8.0 用户手册](https://www.stomics.tech/service/saw_8_1/docs/gao-ji-she-zhi/expression-matrix-format.html)
> 5. [官方GEF文件结构思维导图](https://www.processon.com/view/link/610cc49c7d9c087bbd1ab7ab#map)

目前已完成的功能：

- gem转换为bgef

即将开发的功能：

- 并行化
- cgef转换