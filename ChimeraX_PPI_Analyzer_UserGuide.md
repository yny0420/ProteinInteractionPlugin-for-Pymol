# ChimeraX Protein Interaction Analyzer 使用指南

脚本路径：

```text
/Users/yangyu/Desktop/interfacefinder/ChimeraX_PPI_Analyzer.py
```

## 加载

在 ChimeraX 命令行运行：

```chimerax
open /Users/yangyu/Desktop/interfacefinder/ChimeraX_PPI_Analyzer.py
```

加载后会注册三个命令：

```text
cxppi
cxppi_batch
cxppi_help
```

## 单个 interface 分析

先打开结构，例如：

```chimerax
open 1brs
```

然后分析 #1 模型中的 A-B 链：

```chimerax
cxppi "model=#1 chains=A+B"
```

也可以设置参数：

```chimerax
cxppi "model=#1 chains=A+B cutoff=1.0 hbond_cutoff=3.6 salt_cutoff=4.0"
```

## 批量分析

如果已经加载了多个模型：

```chimerax
cxppi_batch "#1:A+B;#2:A+B;#3:C+D"
```

## 输出内容

日志中会输出：

- Whole interaction area
- Total buried area
- Score
- Score components
- H-bonds
- Salt bridges
- Hydrophobic contacts
- Close contacts
- 每个相互作用的残基、原子、距离、局部面积
- Local interface area hot spots

## 可视化

脚本会：

- 给链 A/B 上不同颜色
- 给界面残基上更深颜色
- 用 pseudobonds 显示：
  - 氢键
  - 盐桥
  - 疏水接触
  - 过近接触/clash

## 面积和评分

```text
Whole interaction area = total buried area / 2
LocalArea = residue dASA / 2
```

Score 是启发式指标：

```text
Score = area + H-bond + salt + hydrophobic - clash penalty
```

它适合比较同一批模型或突变体，不等同于真实结合自由能。

