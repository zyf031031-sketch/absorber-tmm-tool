# Multilayer Anisotropic Absorber Calculator v2

这是多层吸波结构理论计算工具第二版。

## 新增功能

- 支持斜入射；
- 支持 TE / TM 极化；
- 支持单轴各向异性材料；
- 支持 `epsilon_x = epsilon_y != epsilon_z`；
- 支持 `mu_x = mu_y != mu_z`；
- 支持一次计算多个角度，例如 `0, 30, 60`；
- 内置论文五层优化结构材料库模式；
- 输出吸收率、反射率、回波损耗和输入阻抗；
- 可导出 CSV 结果。

## 模型范围

当前模型为：

- 多层均匀平板；
- 金属背板；
- 平面波入射；
- 单轴各向异性材料，光轴沿厚度 z 方向；
- 层顺序：第 1 行靠近空气入射侧，最后 1 行靠近金属背板。

金属背板下：

```text
A = 1 - R = 1 - |S11|^2
```

## 安装

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

安装依赖：

```bash
pip install -r requirements.txt
```

运行：

```bash
python -m streamlit run app.py
```

## 各向异性 CSV 格式

```csv
layer,thickness_mm,eps_xy_real,eps_xy_loss,eps_z_real,eps_z_loss,mu_xy_real,mu_xy_loss,mu_z_real,mu_z_loss
1,1.0,4.0,0.3,3.5,0.2,1.0,0.0,1.0,0.0
2,1.0,8.0,1.0,6.0,0.8,1.0,0.0,1.0,0.0
```

## 损耗符号

程序内部使用：

```text
epsilon = epsilon_real - j * epsilon_loss
mu      = mu_real      - j * mu_loss
```

所以损耗项填正数。
