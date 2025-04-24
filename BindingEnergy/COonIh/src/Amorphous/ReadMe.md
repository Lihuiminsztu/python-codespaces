# Amorphous Ih Studies

## Overview
This files open a new direction of Building ASW(Amorphous Solid Water) work. As of now, i only updated a script which used to change the Ih crystal come from BuildIh function in the other directory "Crystalline".

## Main Features

I am trying to interim building object from Crystalline Ih to Amorphous Ih. I used the tetrahedral method to change the orientation of H atoms connected to bilayer oxygen atoms. In this python script, codes can come ture three main features:
- 1. Moving positions of Hydrogen atoms to special tetrahedral vertices.
- 2. Orient top layer hydrogen facing directions
- 3. Analyze atomic distribution in the system

## Dependecies 
- Atomic Simulation Environment(ASE)
- Numpy 

## Notes

**if u try clone this code, please change the default saved filepath.**

## 中文版本

# 非晶态 Ih 研究

## 概述
本目录开辟了构建非晶态固态水（ASW）模型的新方向。目前，我仅更新了一个脚本，用于修改"Crystalline"目录中 BuildIh 函数生成的 Ih 晶体结构。

## 主要功能

我正在努力将晶态 Ih 转变为非晶态 Ih。当前目标是使用四面体法改变连接到双层氧原子的氢原子的取向。该 Python 脚本目前实现了三个主要功能：
- 1. 在氧原子或氢原子位置添加特殊的四面体顶点
- 2. 调整顶层氢原子的朝向方向
- 3. 分析系统中的原子分布

## 依赖项
- 原子模拟环境 (ASE)
- NumPy

## 注意事项

**如果你克隆此代码，请更改默认的保存文件路径。**