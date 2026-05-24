# Protein Interaction Analyzer 使用指南

本文档适用于 PyMOL 插件：

```text
/Users/yangyu/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

如果把插件分享给别人，对方不需要使用上面的固定路径；只需要把 `ProteinInteractionPlugin.py` 放在自己电脑上的任意文件夹，然后在 PyMOL 中用该文件的实际路径加载即可。

插件用于快速分析蛋白链之间的界面和相互作用，包括：

- dASA 界面残基
- 整体界面面积
- 局部残基面积贡献
- 氢键
- 盐桥/离子相互作用
- 疏水接触
- 过近接触或潜在 clash
- 0-100 的启发式互作强度评分
- 单结构多链两两分析
- 多个 PDB/object 批量分析

## 1. 分享和安装前准备

建议分享整个文件夹：

```text
interfacefinder/
├── ProteinInteractionPlugin.py
├── ProteinInteractionPlugin_UserGuide.md
└── README_ProteinInteractionPlugin.md
```

其中必须文件是：

```text
ProteinInteractionPlugin.py
```

其余 Markdown 文件是说明文档，便于使用者查阅。

### 1.1 运行环境

使用者需要安装：

- PyMOL 2.x
- 能正常打开 PDB/mmCIF 结构文件
- 如果需要从 PDB 在线下载结构，需要 PyMOL 能联网

插件不依赖额外 Python 包，不需要安装 `numpy`、`scipy`、`pandas` 等第三方库。

### 1.2 可选 GUI 依赖

插件提供命令行功能和一个简单图形面板：

```pymol
ppi_panel
```

如果 PyMOL 的 Tk/Tcl 或 XQuartz 环境不可用，面板可能打不开。这个不会影响核心分析功能，继续使用命令行即可：

```pymol
ppi_analyze 1brs, A+B
ppi_batch "001:A+B;002:A+B"
```

macOS 用户如果看到类似：

```text
Tk/Tcl is unavailable
libX11.6.dylib not found
```

说明是 PyMOL 图形面板依赖问题，不是插件分析功能的问题。可以安装或修复 XQuartz / libX11，或者直接使用命令行模式。

### 1.3 分享给别人的最小清单

把下面两个文件发给对方即可：

```text
ProteinInteractionPlugin.py
ProteinInteractionPlugin_UserGuide.md
```

建议让对方把它们放在一个不含中文和特殊符号的简单路径中，例如：

```text
Desktop/interfacefinder/
```

避免放在云盘同步目录、压缩包内部、只读目录或系统应用内部目录。

推荐目录示例：

macOS / Linux：

```text
/Users/yourname/Desktop/interfacefinder/
```

Windows：

```text
C:/Users/yourname/Desktop/interfacefinder/
```

## 2. 安装方式

### 2.1 推荐方式：用 `run` 命令临时加载

这是最稳妥、最容易排错的方式。

假设插件放在：

```text
/Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

在 PyMOL 命令行运行：

```pymol
run /Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

Windows 示例：

```pymol
run C:/Users/yourname/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

路径中如果有空格，建议使用引号：

```pymol
run "/Users/yourname/Desktop/my plugins/interfacefinder/ProteinInteractionPlugin.py"
```

每次重新打开 PyMOL 后，都需要重新运行一次 `run` 命令。

### 2.2 放入 PyMOL 启动目录

如果希望每次打开 PyMOL 自动加载，可以把 `ProteinInteractionPlugin.py` 放到 PyMOL 的启动脚本目录。

常见位置类似：

```text
~/.pymol/startup/
```

或 PyMOL 应用内部的 startup/plugin 目录。

不同系统和 PyMOL 发行版路径可能不同。推荐先用 `run` 命令确认插件能正常工作，再考虑自动加载。

### 2.3 使用 PyMOL Plugin Manager

部分 PyMOL 版本可以通过菜单安装：

```text
Plugin > Plugin Manager > Install New Plugin
```

选择：

```text
ProteinInteractionPlugin.py
```

如果 Plugin Manager 安装后出现路径、权限或 Tk/Tcl 问题，建议退回 `run` 命令方式。

### 2.4 更新插件

如果收到新版 `ProteinInteractionPlugin.py`：

1. 关闭 PyMOL，或至少不要在 PyMOL 中运行旧插件命令
2. 用新版文件覆盖旧的 `ProteinInteractionPlugin.py`
3. 重新打开 PyMOL
4. 重新运行：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin.py
```

如果 PyMOL 中已经有旧结果对象或 selection，可以先清理：

```pymol
delete PPI*
```

然后重新分析。

### 2.5 卸载插件

如果只是用 `run` 命令临时加载，关闭 PyMOL 即可，不需要卸载。

如果放进了启动目录或通过 Plugin Manager 安装：

1. 删除对应目录里的 `ProteinInteractionPlugin.py`
2. 重启 PyMOL
3. 如果仍能调用 `ppi_analyze`，说明还有另一个副本在 PyMOL 启动路径中，需要继续查找并删除

## 3. 首次加载和验证

在 PyMOL 命令行运行：

```pymol
run /Users/yangyu/Desktop/interfacefinder/ProteinInteractionPlugin.py
```

如果你修改过插件文件，需要重新运行上面这条命令，让 PyMOL 重新加载新版脚本。

加载后，测试命令是否存在：

```pymol
ppi_score_help
```

如果控制台打印评分说明，说明插件已经加载成功。

也可以做一个最小测试：

```pymol
fetch 1brs, async=0
ppi_analyze 1brs, A+D
```

成功时会看到：

- PyMOL 控制台输出表格
- 对象/selection 列表中出现 `PPI...` 开头的结果
- 结构界面被着色
- 氢键、盐桥、疏水接触等以 distance object 显示

## 4. 加载结构

### 4.1 从 PDB 下载

```pymol
ppi_fetch 1brs
```

这等价于从 PDB 获取 `1brs` 并加载到 PyMOL。

也可以直接用 PyMOL 自带命令：

```pymol
fetch 1brs, async=0
```

### 4.2 加载本地结构文件

```pymol
ppi_load /path/to/model.pdb
```

支持常见的 `.pdb`、`.cif`、`.mmcif` 文件。

## 5. 分析两条链

最常用命令：

```pymol
ppi_analyze 1brs, A+B
```

含义：

- `1brs`：PyMOL 中的对象名
- `A+B`：分析 A 链和 B 链的界面

如果对象名是数字开头，例如 `001`，也可以这样写：

```pymol
ppi_analyze 001, A+B
```

如果 PyMOL 对参数解析不稳定，可以使用引号：

```pymol
ppi_analyze "001, A+B"
```

## 6. 分析同一结构中的多条链

如果想分析 A、B、C、D 四条链之间所有两两互作：

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

输出为一个 pairwise summary 表格。

## 7. 批量分析多个 PDB/object

如果已经在 PyMOL 中加载了多个对象，例如：

```text
001
002
003
```

可以一键分析每个对象指定的两条链：

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D"
```

含义：

- `001:A+B`：分析对象 `001` 的 A-B 链
- `002:A+B`：分析对象 `002` 的 A-B 链
- `003:C+D`：分析对象 `003` 的 C-D 链

批量分析适合比较多个模型、突变体、构象或不同复合物的界面强弱。

## 8. 导出 CSV 表格

### 6.1 单对链导出

```pymol
ppi_analyze 1brs, A+B, export_path=/tmp/1brs_AB.csv
```

### 6.2 多链两两汇总导出

```pymol
ppi_analyze 1brs, A+B+C+D, summary_csv=/tmp/1brs_pair_summary.csv
```

### 6.3 批量分析汇总导出

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D", summary_csv=/tmp/batch_summary.csv
```

### 6.4 批量导出每一对的详细结果

```pymol
ppi_batch "001:A+B;002:A+B", export_dir=/tmp/ppi_details
```

每个 object/pair 会输出一个独立 CSV 文件。

## 9. PyMOL 中会生成什么

以命令为例：

```pymol
ppi_analyze 001, A+B
```

默认 prefix 是 `PPI`，所以会生成类似：

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

| 名称 | 类型 | 作用 |
|---|---|---|
| `...chainA` | selection | A 链选择，不会复制成新 object |
| `...chainB` | selection | B 链选择，不会复制成新 object |
| `...iface` | selection | dASA 识别的界面残基 |
| `...hbonds` | distance object | 氢键，显示距离 |
| `...salt` | distance object | 盐桥/离子相互作用，显示距离 |
| `...hydrophobic` | distance object | 疏水接触，显示距离 |
| `...clashes` | distance object | 过近接触或潜在 clash，显示距离 |

颜色：

| 相互作用 | 颜色 |
|---|---|
| 氢键 | cyan |
| 盐桥/离子相互作用 | magenta |
| 疏水接触 | yellow |
| close contact / clash | red |
| 链 A | lightblue |
| 链 B | wheat |
| 链 A 界面残基 | marine |
| 链 B 界面残基 | orange |

## 10. 输出表格怎么看

### 8.1 单对分析输出

运行：

```pymol
ppi_analyze 001, A+B
```

会输出详细表格，包括：

```text
Interface residues
Whole interaction area
Total buried area
Score
Score components
H-bonds / Salt bridges / Hydrophobic contacts / Close contacts
Interaction table
Local interface area hot spots
Top interface residues by dASA
```

### 8.2 批量分析输出

运行：

```pymol
ppi_batch "001:A+B;002:A+B;003:C+D"
```

会输出汇总表，例如：

```text
Object  Pair  Residues  Area(A2)  HBond  Salt  Hydrophobic  Clash  Score  AreaSc  PolarSc
001     A-B   ...
002     A-B   ...
003     C-D   ...
```

后面还会追加：

```text
Local interface area hot spots
```

用于显示每个 interface 中局部面积贡献最高的残基。

## 11. 面积指标解释

### 9.1 Total buried area

总埋藏面积，来自两条链界面残基 dASA 的总和。

简单理解：

```text
Total buried area = chain A dASA + chain B dASA
```

### 9.2 Whole interaction area

整体界面面积，插件中估算为：

```text
Whole interaction area = Total buried area / 2
```

这是常见的复合物界面面积近似方式。

### 9.3 LocalArea

局部残基面积贡献，插件中估算为：

```text
LocalArea = residue dASA / 2
```

`LocalArea` 较大的残基通常是该界面的主要贡献残基，可以作为潜在 hot spot 候选。

注意：这是基于结构几何的近似指标，不等同于实验突变能量。

## 12. 评分怎么看

插件给出一个 0-100 的启发式分数：

```text
Score
```

它综合考虑：

| 项目 | 作用 |
|---|---|
| AreaSc | 界面面积贡献，最高 40 |
| H-bond score | 氢键贡献 |
| Salt score | 盐桥贡献 |
| Hydrophobic score | 疏水接触贡献 |
| Clash penalty | 过近接触扣分 |

控制台中也可以运行：

```pymol
ppi_score_help
```

查看评分解释。

### 10.1 如何判断哪个 interface 更强

建议按以下顺序看：

1. 先看 `Score`
2. 再看 `Area(A2)` 或 `Whole interaction area`
3. 再看 `PolarSc`
4. 再看氢键、盐桥数量
5. 再看疏水接触数量
6. 最后检查 `Clash` 是否过多

一般来说，更可靠的强界面通常具有：

- 较大的界面面积
- 较多且合理的氢键/盐桥
- 一定数量的疏水接触
- 较低的 clash 数
- 多个高 LocalArea 残基集中在界面核心

### 10.2 不建议只看一个指标

例如：

- 面积大但 clash 很多，可能是模型不合理
- 氢键多但面积很小，可能只是局部接触
- 疏水接触多但没有极性配对，可能需要结合具体结构判断

所以推荐用 `Score + Area + PolarSc + LocalArea hot spots` 综合判断。

## 13. 常用参数

### 11.1 dASA cutoff

默认：

```text
cutoff=1.0
```

残基 dASA 大于该阈值才算 interface residue。

示例：

```pymol
ppi_analyze 001, A+B, cutoff=2.0
```

### 11.2 氢键距离

默认：

```text
hbond_cutoff=3.6
```

示例：

```pymol
ppi_analyze 001, A+B, hbond_cutoff=3.4
```

### 11.3 盐桥距离

默认：

```text
salt_cutoff=4.0
```

### 11.4 疏水接触距离

默认：

```text
hydrophobic_cutoff=4.2
```

### 11.5 close contact / clash 距离

默认：

```text
clash_cutoff=2.2
```

## 14. 常见问题

### 12.1 GUI 面板打不开

如果运行：

```pymol
ppi_panel
```

出现类似：

```text
Tk/Tcl is unavailable
libX11.6.dylib not found
```

这是 PyMOL 的 Tk/X11 图形依赖问题，不是分析功能本身的问题。

可以继续用命令行：

```pymol
ppi_analyze 001, A+B
ppi_batch "001:A+B;002:A+B"
```

如果想修复 GUI，需要安装或修复 XQuartz / libX11 环境。

### 12.2 对象名以数字开头

例如对象名是：

```text
001
```

插件内部会尽量用安全方式处理。推荐命令：

```pymol
ppi_analyze 001, A+B
```

如果 PyMOL 参数解析异常，可以用：

```pymol
ppi_analyze "001, A+B"
```

### 12.3 批量命令引号问题

推荐：

```pymol
ppi_batch "001:A+B;002:A+B"
```

如果仍然解析异常，可以不用引号：

```pymol
ppi_batch 001:A+B;002:A+B
```

但某些 PyMOL 版本会把分号当命令分隔符，所以通常还是推荐带引号。

### 12.4 结果为 0

如果某一对链结果为 0，可能原因包括：

- 链名写错
- 两条链本身没有接触
- cutoff 过高
- 结构中缺失相关残基或原子
- 分析对象不是你以为的那个 object

可以先检查链：

```pymol
get_chains 001
```

或者在 PyMOL 中单独显示链：

```pymol
show cartoon, 001 and chain A
show cartoon, 001 and chain B
```

## 15. 推荐工作流

### 单个复合物

```pymol
run /Users/yangyu/Desktop/interfacefinder/ProteinInteractionPlugin.py
ppi_fetch 1brs
ppi_analyze 1brs, A+B
```

然后看：

1. PyMOL 视图中的界面残基和相互作用虚线
2. 控制台表格中的面积和相互作用数量
3. `Local interface area hot spots`
4. `Score components`

### 多个模型比较

```pymol
run /Users/yangyu/Desktop/interfacefinder/ProteinInteractionPlugin.py
ppi_batch "001:A+B;002:A+B;003:A+B", summary_csv=/tmp/ppi_compare.csv
```

然后比较：

1. `Score`
2. `Area(A2)`
3. `PolarSc`
4. `Clash`
5. `top_local_area_residues`

## 16. 重要限制

这个插件是结构几何分析工具，不是严格的能量计算软件。

它不能替代：

- MM/GBSA
- Rosetta interface score
- FoldX
- 分子动力学自由能计算
- 实验亲和力测定

适合用途：

- 快速查看界面
- 比较一组相关模型
- 找潜在界面热点残基
- 判断突变前后界面是否明显变差
- 生成交互作用表格用于汇报或进一步分析
