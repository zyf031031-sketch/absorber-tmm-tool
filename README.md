# Multilayer Absorber TMM Tool

这是一个用于计算多层吸波结构理论吸收率的小工具。

## 第一版模型范围

- 多层均匀平板材料；
- 垂直入射；
- 金属背板；
- 最后透射率近似为 0，因此吸收率 `A = 1 - R`；
- 使用传输线 / 传输矩阵思想递推输入阻抗。

## 安装

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

macOS / Linux:

```bash
source .venv/bin/activate
```

安装依赖：

```bash
pip install -r requirements.txt
```

## 运行

```bash
streamlit run app.py
```

浏览器会自动打开工具页面。

## CSV 输入格式

可以上传 CSV，列名如下：

```csv
layer,thickness_mm,eps_real,eps_loss,mu_real,mu_loss
1,1.0,4.0,0.3,1.0,0.0
2,1.0,8.0,1.2,1.0,0.0
```

`layer` 列可以省略。

## 参数说明

程序内部使用：

```text
epsilon_r = eps_real - j * eps_loss
mu_r      = mu_real - j * mu_loss
```

所以损耗项填正数即可。




