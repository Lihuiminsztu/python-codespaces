# Computational Chemistry Projects

This repository contains Python code and computational tips for Astrophysics and Astrochemistry Projects. Calculates structural and chemical parameters of interstellar ices with quantum mechanical methods.

## Project Structure
As for now, there is just one file that can be used for calculations.

- `BindingEnergy/`: Main simulation code for calculating binding energies
  - `COonIh/`: Specific modules for CO adsorption on ice Ih
    - `src/`: Source code including absorption.py and test_with_logging.py etc...
    - `ReadMe.md`: Detailed documentation
  
This project employs GFN1/GFN2-xTB semi-empirical quantum mechanical methods within the TBLite framework to efficiently calculate binding energies of interstellar molecules on ice surfaces, optimizing convergence through modified SCF parameters.

## Future work
- Future integration with machine learning for enhanced predictions.Sametime, i will import DFT method to calculate more models.

## Blog

I've documented my learning journey and the development of this project in a blog:

[**View My Computational Chemistry Learning Blog**](https://lihuiminsztu.github.io/python-codespaces/)

The blog includes detailed explanations of the computational methods, code evolution, and scientific insights gained during this project.

## 中文版本

# 计算化学项目

本仓库包含天体物理学和天体化学项目的 Python 代码和计算技巧。使用量子力学方法计算星际冰的结构和化学参数。

## 项目结构
目前，只有一个文件可用于计算。

- `BindingEnergy/`：计算结合能的主要模拟代码
  - `COonIh/`：CO 在 Ih 冰表面吸附的特定模块
    - `src/`：源代码，包括 absorption.py 和 test_with_logging.py 等
    - `ReadMe.md`：详细文档
  
本项目采用 GFN1/GFN2-xTB 半经验量子力学方法，结合 TBLite 框架，高效计算星际分子在冰表面的结合能，并通过修改 SCF 参数优化收敛性。

## 未来工作
- 未来将整合机器学习以增强预测能力。同时，我将引入 DFT 方法来计算更多模型。

## 博客

我在博客中记录了我的学习历程和项目开发过程：

[**查看我的计算化学学习博客**](https://lihuiminsztu.github.io/python-codespaces/)

该博客包含了计算方法的详细解释、代码演变以及在此项目进行中获得的科学见解。
