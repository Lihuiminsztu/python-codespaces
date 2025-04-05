---
layout: post
title: "CO分子在冰Ih表面的高级定位控制"
date: 2025-04-05
categories: [python, chemistry, simulation]
tags: [computational-chemistry, molecular-dynamics, ase]
---

## 背景介绍

在冰表面吸附模拟中，分子的初始位置和方向对结果有显著影响。今天，我改进了`absorption.py`模块，
增加了更丰富的分子定位和定向选项。

## 关键改进

1. 增加了垂直/水平/倾斜三种分子定向方式
2. 添加了控制哪个原子面向表面的选项（C或O朝向表面）
3. 添加了分子旋转角度控制
4. 支持分子方向反转

## 代码实现示例

```python
def place_molecule_on_slab(slab, molecule, position, orientation='vertical', 
                           angle=0.0, facing_atom=0, invert=False):
    '''Place a molecule at a specific position on the slab with controlled orientation
    
    Parameters:
    ----------
    slab: Atoms
        The surface slab
    molecule: Atoms
        The molecule to place on the slab
    position: array-like
        The position for the molecule's first atom
    orientation: str
        Orientation of the molecule: 'vertical' (default), 'horizontal', or 'tilted'
    angle: float
        Rotation angle in degrees (for 'horizontal' or 'tilted' orientations)
        For 'horizontal': rotation around vertical axis
        For 'tilted': tilt angle from vertical position
    facing_atom: int
        Index of the atom in the molecule that should face the surface (default: 0)
    invert: bool
        If True, inverts the molecule orientation (useful for non-diatomic molecules)
    
    Returns:
    ----------
    Combined: Atoms object
        The combined system with molecule on slab in specified orientation
    '''
    # ...代码细节...
```

## 使用示例

以下命令可以在冰表面模拟C-向下的CO分子:

```bash
python test_with_logging.py --miller 0,0,1 --site ontop --orient vertical --face 0 --dist 2.5
```

而以下命令则使O原子朝向表面:

```bash
python test_with_logging.py --miller 0,0,1 --site ontop --orient vertical --face 1 --dist 2.5
```

## 下一步计划

我计划进一步扩展此功能，支持更复杂分子的定向控制以及多分子系统的构建。
