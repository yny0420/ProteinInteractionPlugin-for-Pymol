# ProteinInteractionPlugin_v3 更新说明

## 版本目的

`ProteinInteractionPlugin_v3.py` 基于 `ProteinInteractionPlugin_v2.py` 修改。

本版本专门解决不同 PyMOL 版本、主题或显示设置下，界面残基 label 不明显或不显示的问题。

## 修改内容

仅修改界面残基 label 的显示设置：

- 强制显示界面残基 CA 的 label
- label 使用单字母残基缩写，例如 `A:K45`
- label 颜色固定为白色
- label 大小固定为 18
- label 描边颜色固定为黑色
- label 位置略微上移，减少被 sticks 或 cartoon 遮挡

## 使用方法

在 PyMOL 中加载：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin_v3.py
```

然后正常运行：

```pymol
ppi_analyze 001, A+B
```

或：

```pymol
ppi_batch "001:A+B;002:A+B"
```

## 如果 label 仍然不可见

先确认右侧对象列表中对应结果的 `L` 按钮没有关闭。

也可以手动运行：

```pymol
show labels, PPI*iface and name CA
set label_color, white
set label_size, 18
set label_outline_color, black
```

## 兼容性

v3 保留 v2 的兼容性修复：

- 修复部分 PyMOL 版本中 `Atom` object 没有 `model` 字段导致的报错
- 保留原有分析命令、评分、局部面积和批量分析功能

