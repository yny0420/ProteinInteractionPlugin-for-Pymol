# Protein Interaction Analyzer 用户指南

本指南面向所有 `ProteinInteractionPlugin.py` 脚本的用户。  
不需要懂 Python，也不需要安装额外 Python 包；只需要会在 PyMOL 中打开结构并输入命令。

## 1. 这个插件能做什么

Protein Interaction Analyzer 是一个 PyMOL 插件，用于分析蛋白链之间的相互作用界面。

主要功能：

- 识别蛋白链之间的界面残基
- 主界面分析会把常见 PTM 残基作为蛋白残基一起纳入计算
- 计算整体界面面积
- 计算每个界面残基的局部面积贡献
- 识别氢键
- 识别盐桥/离子相互作用
- 识别疏水接触
- 标记过近接触或潜在 clash
- 给出一个 0-100 的启发式互作强度评分
- 支持一个结构内多条链的两两分析
- 支持多个 PDB/object 的批量比较
- 支持从自定义区域或 PTM 位点出发搜索周围互作
- 支持泛素链界面分析和可能的泛素化连接位点识别
- 支持导出 CSV 表格

## 2. 文件清单


```text
ProteinInteractionPlugin.py
```
```text
USER_GUIDE.md
```

```

## 3. 环境要求

需要：

- PyMOL 2.x
- 可以正常打开 PDB/mmCIF 文件
- 如果需要在线下载 PDB，PyMOL 需要联网

不需要：

- 不需要安装 `numpy`
- 不需要安装 `pandas`
- 不需要安装 `scipy`
- 不需要配置 Conda 环境

## 4. 安装与加载

### 4.1 推荐方式：用 `run` 命令加载

把 `ProteinInteractionPlugin.py` 放到任意文件夹，例如：

macOS / Linux：

```text
/Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

Windows：

```text
C:/Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

然后在 PyMOL 命令行运行：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin.py
```

macOS / Linux 示例：

```pymol
run /Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

Windows 示例：

```pymol
run C:/Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

如果路径中有空格，请加引号：

```pymol
run "/Users/yourname/Desktop/my plugins/interfacefinder/ProteinInteractionPlugin.py"
```

每次重新打开 PyMOL 后，都需要重新运行一次 `run` 命令。

### 4.2 验证是否加载成功

加载后运行：

```pymol
ppi_score_help
```

如果 PyMOL 控制台输出评分说明，表示插件加载成功。

也可以运行：

```pymol
ppi_fetch 1brs
ppi_analyze 1brs, A+D
```

如果出现表格输出，并且结构界面被着色，说明插件工作正常。

## 5. 加载结构

### 5.1 从 PDB 下载结构

```pymol
ppi_fetch 1brs
```

也可以使用 PyMOL 自带命令：

```pymol
fetch 1brs, async=0
```

### 5.2 打开本地结构文件

```pymol
ppi_load /path/to/model.pdb
```

支持常见结构文件：

- `.pdb`
- `.cif`
- `.mmcif`

也可以直接用 PyMOL 的 `load` 命令打开结构。

## 6. 分析两条链

假设已经加载了对象 `1brs`，并想分析 A 链和 B 链：

```pymol
ppi_analyze 1brs, A+B
```

含义：

```text
1brs  = PyMOL 中的对象名
A+B   = 分析 A 链和 B 链
```

如果对象名是数字开头，例如：

```text
001
```

可以运行：

```pymol
ppi_analyze 001, A+B
```

如果你的 PyMOL 对逗号参数解析异常，可以使用引号：

```pymol
ppi_analyze "001, A+B"
```

## 7. 分析一个结构中的多条链

如果想分析 A、B、C、D 四条链的所有两两互作：

```pymol
ppi_analyze 1brs, A+B+C+D
```

插件会自动分析：

```text
A-B
A-C
A-D
B-C
B-D
C-D
```

适合用来判断一个多聚体中哪两个链的界面最强。

## 8. 批量分析多个结构

如果 PyMOL 中已经加载了多个对象：

```text
001
002
003
```

可以批量分析：

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D"
```

含义：

```text
001:A+B  分析对象 001 的 A-B 链
002:A+B  分析对象 002 的 A-B 链
003:C+D  分析对象 003 的 C-D 链
```

批量分析适合比较：

- 多个预测模型
- 多个突变体
- 多个构象
- 多个 docking 结果
- 多个同源复合物

## 9. 导出 CSV

### 9.1 单对分析导出

```pymol
ppi_analyze 1brs, A+B, export_path=/tmp/1brs_AB.csv
```

### 9.2 多链两两分析导出

```pymol
ppi_analyze 1brs, A+B+C+D, summary_csv=/tmp/1brs_pairs.csv
```

### 9.3 批量分析导出

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D", summary_csv=/tmp/batch_summary.csv
```

### 9.4 批量导出详细结果

```pymol
ppi_batch "001:A+B;002:A+B", export_dir=/tmp/ppi_details
```

## 9.5 分析自定义区域周围的互作位点

如果你想选中某个 loop、口袋、motif、突变位点或一段残基，然后查看其它任何与该区域互作的位点，可以使用 `ppi_region`。

推荐先创建一个 PyMOL selection：

```pymol
select my_region, 1brs and chain A and resi 35-45
```

然后运行：

```pymol
ppi_region my_region
```

插件会自动搜索 `my_region` 之外所有与它接触的蛋白残基，并显示：

- 氢键
- 盐桥/离子相互作用
- 疏水接触
- 一般接触
- close contact / clash

也可以显式指定搜索范围：

```pymol
ppi_region "query=my_region; target=1brs and chain D"
```

导出 CSV：

```pymol
ppi_region "query=my_region; target=1brs and chain D; export_path=/tmp/region_contacts.csv"
```

常见用途：

- 查看某个突变位点周围所有互作残基
- 查看某个结构域、loop 或 motif 与其它链的接触
- 查看小范围界面附近是否有氢键、盐桥或疏水接触
- 从一个功能区域出发寻找潜在互作热点

## 9.6 分析 PTM 位点和蛋白其它区域的互作

插件会在主界面分析中自动把常见 PTM 残基作为蛋白残基处理，包括面积计算、界面残基显示、氢键/盐桥/接触分析和批量分析。也就是说，`ppi_analyze` 和 `ppi_batch` 从一开始就会纳入这些修饰位点。

如果想专门从 PTM 位点出发搜索周围互作，可以额外使用 `ppi_ptm`。当前内置识别的常见修饰包括：

- 磷酸化：`SEP`、`TPO`、`PTR`
- 乙酰化：`ALY`
- 甲基化：`MLY`、`M3L`、`MLZ`、`M2L`
- 羟基化/氧化等常见修饰：`HYP`、`CSO`、`CME`、`KCX`

盲找某个 object 中所有 PTM 位点与其它蛋白区域的互作：

```pymol
ppi_ptm 7mp9
```

只看磷酸化位点：

```pymol
ppi_ptm "object=7mp9; ptm=phospho"
```

限定某条链或某个位点：

```pymol
ppi_ptm "object=7mp9; chain=A"
ppi_ptm "object=7mp9; chain=A; resi=205"
```

如果你的 PTM 残基名称比较特殊，也可以直接指定残基名：

```pymol
ppi_ptm "object=your_object; ptm=SEP+TPO+PTR"
```

或者先自己建立 selection，再把这个 selection 当作 PTM 位点来分析：

```pymol
select my_ptm, your_object and chain A and resi 205
ppi_ptm "query=my_ptm"
```

导出 CSV：

```pymol
ppi_ptm "object=7mp9; ptm=phospho; export_path=/tmp/ptm_contacts.csv"
```

输出和显示内容包括 PTM 位点、周围互作残基、氢键、盐桥/离子相互作用、疏水接触、一般接触、close contact/clash 以及对应距离。

对于同一条链上的 PTM，`ppi_ptm` 默认会排除 PTM 前后 1 个残基，避免把肽键相邻残基误判为互作。如果确实想保留相邻残基，可以设置：

```pymol
ppi_ptm "object=7mp9; exclude_neighbors=0"
```

## 9.7 分析泛素化/泛素链互作

泛素化在结构中通常不是一个小的三字母 PTM 残基，而是一条或多条 ubiquitin protein chain。因此插件增加了 `ppi_ub`，用于自动识别 ubiquitin-like 链，分析这些泛素链和其它链之间的界面，并寻找可能的泛素化连接位点。

自动识别结构中的泛素链并分析：

```pymol
ppi_ub 2o6v
```

插件会输出：

- 检测到的 ubiquitin-like chains
- 每条泛素链和其它链的界面面积
- 氢键、盐桥、疏水接触、close contact
- 互作强度评分
- 可能的 Gly76-C 到 Lys-NZ 异肽连接
- 可能的 linear ubiquitin N-terminal linkage

如果想手动指定哪些链是泛素链：

```pymol
ppi_ub "object=your_object; ub_chains=A+B"
```

如果只想分析泛素链和指定 target chains：

```pymol
ppi_ub "object=your_object; ub_chains=A+B; target_chains=C+D"
```

导出汇总 CSV：

```pymol
ppi_ub "object=2o6v; summary_csv=/tmp/ubiquitin_summary.csv"
```

导出每个泛素链-目标链界面的详细 CSV：

```pymol
ppi_ub "object=2o6v; export_dir=/tmp/ubiquitin_details"
```

调整异肽键距离阈值：

```pymol
ppi_ub "object=2o6v; linkage_cutoff=2.8"
```

说明：

- `ppi_ub` 会把 ubiquitin chain 当作蛋白链来做界面分析。
- 如果结构中真实建模了 Gly76-C 到 Lys-NZ 的连接，插件会在 linkage 表中列出。
- 有些 PDB 文字说明中存在泛素化连接，但连接柔性太大或没有建模到电子密度中，此时插件可以分析界面，但不会报告 covalent linkage。

## 10. PyMOL 中生成的结果

以命令为例：

```pymol
ppi_analyze 001, A+B
```

默认结果名称会以 `PPI...` 开头，例如：

```text
PPI001ABchainA
PPI001ABchainB
PPI001ABiface
PPI001ABhbonds
PPI001ABsalt
PPI001ABhydrophobic
PPI001ABclashes
```

说明：

| 名称 | 类型 | 含义 |
|---|---|---|
| `...chainA` | selection | A 链选择，不会复制成新 object |
| `...chainB` | selection | B 链选择，不会复制成新 object |
| `...iface` | selection | 界面残基 |
| `...hbonds` | distance object | 氢键 |
| `...salt` | distance object | 盐桥/离子相互作用 |
| `...hydrophobic` | distance object | 疏水接触 |
| `...clashes` | distance object | 过近接触或潜在 clash |

颜色：

| 类型 | 颜色 |
|---|---|
| 链 A | lightblue |
| 链 B | wheat |
| 链 A 界面残基 | marine |
| 链 B 界面残基 | orange |
| 氢键 | cyan |
| 盐桥 | magenta |
| 疏水接触 | yellow |
| close contact / clash | red |

## 11. 输出表格怎么看

### 11.1 整体面积

输出中会看到：

```text
Whole interaction area
Total buried area
```

解释：

```text
Total buried area = 两条链界面残基 dASA 总和
Whole interaction area = Total buried area / 2
```

`Whole interaction area` 可以用来粗略比较界面大小。

### 11.2 局部面积

输出中会看到：

```text
Local interface area hot spots
```

其中：

```text
LocalArea = residue dASA / 2
```

LocalArea 较大的残基，通常是界面中贡献较大的局部残基，可作为潜在 hot spot 候选。

### 11.3 相互作用类型

表格会列出：

- H-bond
- Salt bridge
- Hydrophobic
- Close contact

每一行包含：

- 相互作用类型
- 两侧残基
- 两侧原子
- 距离
- 局部面积估计

## 12. 评分怎么理解

插件会给出一个 0-100 的启发式评分：

```text
Score
```

评分组成：

| 指标 | 含义 |
|---|---|
| AreaSc | 界面面积贡献，最高 40 |
| PolarSc | 氢键 + 盐桥贡献 |
| Hydrophobic | 疏水接触贡献 |
| ClashPenalty | 过近接触扣分 |

可以运行：

```pymol
ppi_score_help
```

查看评分解释。

### 12.1 如何判断哪个界面更强

建议综合看：

1. `Score`
2. `Whole interaction area`
3. `PolarSc`
4. 氢键数量
5. 盐桥数量
6. 疏水接触数量
7. clash 是否过多
8. LocalArea hot spots 是否集中

一般来说，更强、更可信的界面通常表现为：

- 面积较大
- 氢键/盐桥较多
- 疏水接触合理
- clash 较少
- 有多个较大的 LocalArea 残基

注意：评分是结构几何启发式指标，不等同于真实结合自由能。

## 13. 常用参数

### dASA 阈值

```pymol
ppi_analyze 1brs, A+B, cutoff=1.0
```

### 氢键距离阈值

```pymol
ppi_analyze 1brs, A+B, hbond_cutoff=3.6
```

### 盐桥距离阈值

```pymol
ppi_analyze 1brs, A+B, salt_cutoff=4.0
```

### 疏水接触距离阈值

```pymol
ppi_analyze 1brs, A+B, hydrophobic_cutoff=4.2
```

### clash 距离阈值

```pymol
ppi_analyze 1brs, A+B, clash_cutoff=2.2
```

## 14. 常见问题

### 14.1 `ppi_panel` 打不开

如果看到：

```text
Tk/Tcl is unavailable
libX11.6.dylib not found
```

说明 PyMOL 的图形面板依赖不可用。

解决方式：

- 继续使用命令行模式
- 或安装/修复 XQuartz / libX11

命令行模式不受影响。

### 14.2 结果为 0

可能原因：

- 链名写错
- 两条链不接触
- cutoff 设置过高
- 结构缺失残基或原子
- 分析对象不是你想分析的那个 object

可以检查链：

```pymol
get_chains 1brs
```

也可以手动显示链：

```pymol
show cartoon, 1brs and chain A
show cartoon, 1brs and chain B
```

### 14.3 批量命令解析出错

推荐格式：

```pymol
ppi_batch "001:A+B;002:A+B"
```

如果对象名或路径复杂，尽量避免空格和特殊符号。

### 14.4 `Atom` object has no attribute `model`

如果看到类似报错：

```text
AttributeError: 'Atom' object has no attribute 'model'
```

这是 PyMOL 版本差异导致的。某些 PyMOL 版本中 `cmd.get_model()` 返回的 Atom 对象没有 `model` 字段。

解决方法：使用兼容版插件：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin_v2.py
```

然后按原来的命令运行：

```pymol
ppi_analyze 001, A+B
ppi_batch "001:A+B;002:A+B"
```

更多说明见：

```text
ProteinInteractionPlugin_v2_CompatibilityNote.md
```

### 14.5 重新运行后对象太多

可以清理旧结果：

```pymol
delete PPI*
```

然后重新运行分析。

## 15. 推荐工作流

### 单个结构

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin.py
ppi_fetch 1brs
ppi_analyze 1brs, A+D
```

然后检查：

1. PyMOL 中的界面残基
2. 氢键/盐桥/疏水接触
3. Whole interaction area
4. Score
5. LocalArea hot spots

### 多个模型比较

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin.py
ppi_batch "001:A+B;002:A+B;003:A+B", summary_csv=/tmp/ppi_compare.csv
```

然后比较：

1. `Score`
2. `Area(A2)`
3. `PolarSc`
4. `Clash`
5. `top_local_area_residues`

## 16. 方法限制

这个插件适合快速结构分析和模型比较，但不是严格的结合能计算工具。

它不能替代：

- 实验亲和力测定
- Rosetta interface score
- FoldX
- MM/GBSA
- 分子动力学自由能计算

适合用途：

- 快速查看蛋白界面
- 比较同一批结构或突变体
- 找潜在界面热点残基
- 判断界面是否明显变强或变弱
- 生成相互作用表格用于汇报或后续分析
