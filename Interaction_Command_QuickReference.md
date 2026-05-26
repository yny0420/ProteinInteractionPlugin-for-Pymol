# Protein Interaction Analyzer 互作分析命令速查

Copyright (c) yangyu. All rights reserved.

本文档总结插件中用于寻找和分析蛋白互作的主要命令。

## 1. `ppi_analyze`

最主要的分析命令，用来分析两条链或多条链之间的蛋白互作。

```pymol
ppi_analyze 1brs, A+B
```

含义：

- `1brs`：PyMOL 中的 object 名称
- `A+B`：分析 A 链和 B 链之间的互作

会输出和显示：

- 界面残基
- 整体互作面积
- 每个界面残基的局部面积贡献
- 氢键
- 盐桥/离子相互作用
- 疏水接触
- close contact / clash
- 互作强度评分
- PTM 位点，常见修饰残基会自动纳入主分析

分析多条链的所有两两互作：

```pymol
ppi_analyze 1brs, A+B+C+D
```

会自动分析：

```text
A-B, A-C, A-D, B-C, B-D, C-D
```

适合用来判断一个多聚体中哪两个链的互作界面更强。

## 2. `ppi_batch`

批量分析多个结构或多个 object。

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D"
```

含义：

- `001:A+B`：分析 object `001` 中 A-B 链互作
- `002:A+B`：分析 object `002` 中 A-B 链互作
- `003:C+D`：分析 object `003` 中 C-D 链互作

适合比较：

- 多个突变体
- 多个 AlphaFold/ColabFold 模型
- 多个 docking 结果
- 多个构象
- 多个同源复合物

## 3. `ppi_region`

分析一个自定义区域和蛋白其它部位的互作。

先创建 PyMOL selection：

```pymol
select my_region, 1brs and chain A and resi 35-45
```

再分析：

```pymol
ppi_region my_region
```

插件会盲找 `my_region` 之外所有与该区域接触的蛋白残基。

也可以限定搜索范围：

```pymol
ppi_region "query=my_region; target=1brs and chain D"
```

适合：

- 查看某个 loop 和其它区域的互作
- 查看某个突变位点周围互作
- 查看某个结构域、口袋或 motif 周围接触
- 从局部区域出发寻找潜在互作热点

输出内容包括：

- 氢键
- 盐桥/离子相互作用
- 疏水接触
- 一般接触
- close contact / clash
- 互作残基和距离

## 4. `ppi_ptm`

专门从 PTM 位点出发，寻找 PTM 和蛋白其它区域的互作。

```pymol
ppi_ptm 7mp9
```

作用：自动识别 object `7mp9` 中的常见 PTM 位点，并搜索它们周围的互作残基。

只看磷酸化位点：

```pymol
ppi_ptm "object=7mp9; ptm=phospho"
```

限定某条链：

```pymol
ppi_ptm "object=7mp9; chain=A"
```

限定某个位点：

```pymol
ppi_ptm "object=7mp9; chain=A; resi=205"
```

如果你的 PTM 残基名称比较特殊，可以直接指定残基名：

```pymol
ppi_ptm "object=your_object; ptm=SEP+TPO+PTR"
```

也可以先自己创建 selection：

```pymol
select my_ptm, your_object and chain A and resi 205
ppi_ptm "query=my_ptm"
```

目前内置识别的常见 PTM 包括：

- 磷酸化：`SEP`、`TPO`、`PTR`
- 乙酰化：`ALY`
- 甲基化：`MLY`、`M3L`、`MLZ`、`M2L`
- 其它常见修饰：`HYP`、`CSO`、`CME`、`KCX`

注意：`ppi_ptm` 默认会排除同一条链上 PTM 前后 1 个残基，避免把肽键相邻残基误判为互作。如果需要保留相邻残基：

```pymol
ppi_ptm "object=7mp9; exclude_neighbors=0"
```

## 5. `ppi_ub`

分析泛素链和其它蛋白链之间的互作，并尝试识别可能的泛素化连接位点。

```pymol
ppi_ub 2o6v
```

作用：自动识别结构中的 ubiquitin-like chains，分析它们和其它链之间的界面。

会输出：

- 检测到的泛素链
- 泛素链和其它链的界面面积
- 氢键、盐桥、疏水接触、close contact
- 互作强度评分
- 可能的 Gly76-C 到 Lys-NZ 异肽键
- 可能的 linear ubiquitin N-terminal linkage

手动指定泛素链：

```pymol
ppi_ub "object=your_object; ub_chains=A+B"
```

指定泛素链和 target chains：

```pymol
ppi_ub "object=your_object; ub_chains=A+B; target_chains=C+D"
```

导出汇总：

```pymol
ppi_ub "object=2o6v; summary_csv=/tmp/ubiquitin_summary.csv"
```

导出每个泛素链-目标链界面的详细表格：

```pymol
ppi_ub "object=2o6v; export_dir=/tmp/ubiquitin_details"
```

调整泛素化连接识别距离：

```pymol
ppi_ub "object=2o6v; linkage_cutoff=2.8"
```

注意：如果 PDB 说明中有泛素化连接，但该连接没有被建模到坐标中，`ppi_ub` 仍然可以分析界面，但不会报告 covalent linkage。

## 6. 导出结果

### 6.1 单对链互作详细导出

```pymol
ppi_analyze 1brs, A+B, export_path=/tmp/1brs_AB.csv
```

导出内容包括：

- 总界面面积
- 总 buried area
- 评分
- 氢键、盐桥、疏水接触、close contact
- 每个互作的残基、原子、距离
- 局部界面面积

### 6.2 多链两两分析汇总导出

```pymol
ppi_analyze 1brs, A+B+C+D, summary_csv=/tmp/1brs_pairs.csv
```

### 6.3 批量分析汇总导出

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D", summary_csv=/tmp/batch_summary.csv
```

### 6.4 批量导出详细结果

```pymol
ppi_batch "001:A+B;002:A+B", export_dir=/tmp/ppi_details
```

### 6.5 区域互作导出

```pymol
ppi_region "query=my_region; target=1brs and chain D; export_path=/tmp/region_contacts.csv"
```

### 6.6 PTM 互作导出

```pymol
ppi_ptm "object=7mp9; ptm=phospho; export_path=/tmp/ptm_contacts.csv"
```

### 6.7 泛素链互作导出

```pymol
ppi_ub "object=2o6v; summary_csv=/tmp/ubiquitin_summary.csv"
```

## 7. 辅助命令

### 7.1 在线下载 PDB

```pymol
ppi_fetch 1brs
```

等价于从 PDB 下载 `1brs` 并加载到 PyMOL。

### 7.2 加载本地结构

```pymol
ppi_load /path/to/model.pdb
```

支持常见结构文件：

- `.pdb`
- `.cif`
- `.mmcif`

### 7.3 查看评分解释

```pymol
ppi_score_help
```

会解释插件中的互作强度评分如何构成。

### 7.4 打开图形界面

```pymol
ppi_panel
```

如果 PyMOL 的 Tk/Tcl 或 X11 环境不可用，图形界面可能打不开。此时建议直接使用本文档中的命令行方式。

## 8. 推荐使用顺序

如果只是分析一个复合物中两条链：

```pymol
ppi_fetch 1brs
ppi_analyze 1brs, A+B
```

如果想比较多条链之间哪个界面更强：

```pymol
ppi_analyze 1brs, A+B+C+D
```

如果想比较多个模型或突变体：

```pymol
ppi_batch "model1:A+B;model2:A+B;model3:A+B"
```

如果想看某一段区域周围的互作：

```pymol
select my_region, your_object and chain A and resi 100-130
ppi_region my_region
```

如果想专门看 PTM 位点周围互作：

```pymol
ppi_ptm "object=your_object; ptm=phospho"
```

如果想专门看泛素链互作：

```pymol
ppi_ub "object=your_object"
```
