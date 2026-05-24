# ProteinInteractionPlugin_v2 兼容性说明

## 适用情况

如果运行原版插件时出现下面的错误：

```text
AttributeError: 'Atom' object has no attribute 'model'
```

例如：

```text
File "ProteinInteractionPlugin.py", line ..., in _atom_sel
    return "(model %s and index %d)" % (atom.model, atom.index)
AttributeError: 'Atom' object has no attribute 'model'
```

请使用兼容版：

```text
ProteinInteractionPlugin_v2.py
```

## 报错原因

不同 PyMOL 版本中，`cmd.get_model()` 返回的 Atom 对象字段不完全一致。

有些版本的 Atom 对象包含：

```text
atom.model
```

有些版本没有这个字段。原版插件在创建 distance object 时使用了 `atom.model` 来构造选择语句，所以在这些 PyMOL 版本中会报错。

## v2 修复内容

`ProteinInteractionPlugin_v2.py` 不再假设 Atom 对象一定有 `model` 字段。

它会：

1. 从原始 selection 中解析 object/model 名
2. 给读取到的 Atom 对象临时记录来源 model
3. 如果仍然没有 model 信息，则退回到 index 选择

这样可以兼容更多 PyMOL 版本。

## 如何使用 v2

在 PyMOL 中加载：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin_v2.py
```

例如 Windows：

```pymol
run E:/liangzhu_Project/interface_analysis/ProteinInteractionPlugin_v2.py
```

然后按原来的方式运行：

```pymol
ppi_analyze 001, A+B
```

或：

```pymol
ppi_batch "001:A+B;002:A+B"
```

## 是否需要删除原版

不需要。

为了避免混淆，建议每次打开 PyMOL 后只加载一个版本：

```pymol
run /path/to/interfacefinder/ProteinInteractionPlugin_v2.py
```

如果已经加载过原版，再加载 v2，PyMOL 会用 v2 的命令覆盖同名命令。

## 建议

如果不确定使用的 PyMOL 版本，建议优先加载 v2。

