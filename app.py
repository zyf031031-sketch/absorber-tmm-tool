import io
import re
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

C0 = 299_792_458.0
Z0 = 376.730313668

st.set_page_config(page_title="Multilayer Anisotropic Absorber Calculator", layout="wide")

def parse_angles(text):
    out = []
    for x in re.split(r"[,，;；\s]+", text.strip()):
        if x:
            out.append(float(x))
    return out or [0.0]

def cplx(real, loss):
    return np.asarray(real, dtype=complex) - 1j * np.asarray(loss, dtype=complex)

def sqrt_branch(z):
    q = np.sqrt(z.astype(complex) if isinstance(z, np.ndarray) else complex(z))
    if isinstance(q, np.ndarray):
        q = np.where(np.imag(q) > 0, -q, q)
        q = np.where((np.abs(np.imag(q)) < 1e-14) & (np.real(q) < 0), -q, q)
    else:
        if q.imag > 0:
            q = -q
        if abs(q.imag) < 1e-14 and q.real < 0:
            q = -q
    return q

def paper_material(num, f_ghz):
    f = np.maximum(np.asarray(f_ghz, dtype=float), 1e-12)
    if int(num) == 8:      # Mg-Zn ferrite
        er = np.full_like(f, 12.0)
        ei = np.zeros_like(f)
        mr = 60.0 / (f ** 0.008)
        mi = 15.0 / (f ** 0.015)
    elif int(num) == 10:   # Fe-Si-Al
        er = np.full_like(f, 2.0)
        ei = np.full_like(f, 0.1)
        mr = 80.0 / (f ** 0.003)
        mi = 5.0 / (f ** 0.005)
    elif int(num) == 13:   # Polyaniline
        er = 5.0 / (f ** 0.861)
        ei = 8.0 / (f ** 0.569)
        mr = np.ones_like(f)
        mi = np.zeros_like(f)
    else:
        raise ValueError("Paper library currently supports material 8, 10 and 13 only.")
    eps = cplx(er, ei)
    mu = cplx(mr, mi)
    return eps, eps, mu, mu

def kz_zc(freq_ghz, angle_deg, pol, eps_xy, eps_z, mu_xy, mu_z):
    f = np.asarray(freq_ghz, dtype=float)
    k0 = 2 * np.pi * f * 1e9 / C0
    s = np.sin(np.deg2rad(angle_deg))
    pol = pol.upper()
    if pol == "TE":
        q2 = mu_xy * eps_xy - (mu_xy / mu_z) * s**2
        q = sqrt_branch(q2)
        zc = Z0 * mu_xy / q
    else:
        q2 = eps_xy * mu_xy - (eps_xy / eps_z) * s**2
        q = sqrt_branch(q2)
        zc = Z0 * q / eps_xy
    return k0 * q, zc

def z_air(angle_deg, pol):
    c = np.cos(np.deg2rad(angle_deg))
    if abs(c) < 1e-12:
        raise ValueError("Angle is too close to 90 degrees.")
    return Z0 / c if pol.upper() == "TE" else Z0 * c

def response(freq, layers, mode, angle, pol):
    zin_load = np.zeros_like(freq, dtype=complex)
    for _, row in layers.iloc[::-1].iterrows():
        d = float(row["thickness_mm"]) * 1e-3
        if mode == "Paper 2025 optimized stack":
            eps_xy, eps_z, mu_xy, mu_z = paper_material(row["material_number"], freq)
        elif mode == "Custom isotropic":
            eps = cplx(float(row["eps_real"]), float(row["eps_loss"]))
            mu = cplx(float(row["mu_real"]), float(row["mu_loss"]))
            eps_xy = np.full_like(freq, eps, dtype=complex)
            eps_z = np.full_like(freq, eps, dtype=complex)
            mu_xy = np.full_like(freq, mu, dtype=complex)
            mu_z = np.full_like(freq, mu, dtype=complex)
        else:
            eps_xy_v = cplx(float(row["eps_xy_real"]), float(row["eps_xy_loss"]))
            eps_z_v = cplx(float(row["eps_z_real"]), float(row["eps_z_loss"]))
            mu_xy_v = cplx(float(row["mu_xy_real"]), float(row["mu_xy_loss"]))
            mu_z_v = cplx(float(row["mu_z_real"]), float(row["mu_z_loss"]))
            eps_xy = np.full_like(freq, eps_xy_v, dtype=complex)
            eps_z = np.full_like(freq, eps_z_v, dtype=complex)
            mu_xy = np.full_like(freq, mu_xy_v, dtype=complex)
            mu_z = np.full_like(freq, mu_z_v, dtype=complex)

        kz, zc = kz_zc(freq, angle, pol, eps_xy, eps_z, mu_xy, mu_z)
        t = np.tan(kz * d)
        zin = zc * (zin_load + 1j * zc * t) / (zc + 1j * zin_load * t)
        zin_load = zin

    gamma = (zin_load - z_air(angle, pol)) / (zin_load + z_air(angle, pol))
    R = np.abs(gamma) ** 2
    A = np.clip(1 - R, 0, 1)
    RL = 10 * np.log10(np.maximum(R, 1e-15))
    RL = np.where(np.abs(RL) < 1e-10, 0, RL)
    return zin_load, gamma, R.real, A.real, RL.real

def default_df(mode, n):
    if mode == "Paper 2025 optimized stack":
        return pd.DataFrame([
            {"layer": 1, "material_number": 13, "material_name": "Polyaniline", "thickness_mm": 1.996},
            {"layer": 2, "material_number": 10, "material_name": "Fe-Si-Al", "thickness_mm": 0.288},
            {"layer": 3, "material_number": 13, "material_name": "Polyaniline", "thickness_mm": 1.998},
            {"layer": 4, "material_number": 10, "material_name": "Fe-Si-Al", "thickness_mm": 1.219},
            {"layer": 5, "material_number": 8,  "material_name": "Mg-Zn ferrite", "thickness_mm": 1.999},
        ])
    if mode == "Custom isotropic":
        return pd.DataFrame([{
            "layer": i + 1, "thickness_mm": 1.0,
            "eps_real": 4.0, "eps_loss": 0.3,
            "mu_real": 1.0, "mu_loss": 0.0,
        } for i in range(int(n))])
    return pd.DataFrame([{
        "layer": i + 1, "thickness_mm": 1.0,
        "eps_xy_real": 4.0, "eps_xy_loss": 0.3,
        "eps_z_real": 4.0, "eps_z_loss": 0.3,
        "mu_xy_real": 1.0, "mu_xy_loss": 0.0,
        "mu_z_real": 1.0, "mu_z_loss": 0.0,
    } for i in range(int(n))])

def required_cols(mode):
    if mode == "Paper 2025 optimized stack":
        return ["material_number", "thickness_mm"]
    if mode == "Custom isotropic":
        return ["thickness_mm", "eps_real", "eps_loss", "mu_real", "mu_loss"]
    return ["thickness_mm", "eps_xy_real", "eps_xy_loss", "eps_z_real", "eps_z_loss",
            "mu_xy_real", "mu_xy_loss", "mu_z_real", "mu_z_loss"]

st.title("Multilayer Anisotropic Absorber Calculator")
st.caption("Metal-backed multilayer absorber calculator: oblique incidence, TE/TM, and uniaxial anisotropy.")

with st.expander("模型说明", expanded=True):
    st.markdown(
        """
### 1. 层顺序

表格中的第 1 行表示最靠近空气入射侧的材料层，最后 1 行表示最靠近金属背板的材料层。

### 2. 损耗参数约定

本程序内部采用如下复数形式：

ε = ε' - jε''

μ = μ' - jμ''

因此，在表格中填写损耗项时，损耗参数填正数即可。

### 3. 单轴各向异性材料

本程序目前考虑单轴各向异性材料，即：

εx = εy = εxy

εz 单独设置

μx = μy = μxy

μz 单独设置

其中：

- eps_xy_real、eps_xy_loss 表示 x/y 方向介电常数的实部和损耗；
- eps_z_real、eps_z_loss 表示 z 方向介电常数的实部和损耗；
- mu_xy_real、mu_xy_loss 表示 x/y 方向磁导率的实部和损耗；
- mu_z_real、mu_z_loss 表示 z 方向磁导率的实部和损耗。

### 4. 垂直入射和斜入射的区别

在垂直入射时，电磁波主要受到横向参数 εxy 和 μxy 的影响。

当考虑斜入射时，尤其是 TM 极化情况下，z 方向参数 εz 和 μz 的影响会更加明显。
"""
    )


st.sidebar.header("Frequency")
f1 = st.sidebar.number_input("Start frequency / GHz", min_value=0.001, value=0.1, step=0.1, format="%.3f")
f2 = st.sidebar.number_input("Stop frequency / GHz", min_value=0.001, value=10.0, step=0.5, format="%.3f")
nf = st.sidebar.number_input("Number of points", min_value=11, max_value=10001, value=501, step=10)

st.sidebar.header("Incidence")
angle_text = st.sidebar.text_input("Incident angles / degree", "0, 30, 60")
pol_choice = st.sidebar.selectbox("Polarization", ["TE", "TM", "Both"], index=0)

st.sidebar.header("Material")
mode = st.sidebar.selectbox("Input mode", ["Custom isotropic", "Custom uniaxial anisotropic", "Paper 2025 optimized stack"], index=1)
n_layers = 5 if mode == "Paper 2025 optimized stack" else st.sidebar.number_input("Number of layers", 1, 30, 2, 1)

st.subheader("Layer Parameters")
if mode == "Paper 2025 optimized stack":
    st.info("This uses the paper stack: 13 / 10 / 13 / 10 / 8, with thicknesses 1.996, 0.288, 1.998, 1.219, 1.999 mm. Material parameters are frequency-dependent.")
else:
    st.write("Edit the table below. Row 1 is the air side; the last row is the metal-backed side.")

df0 = default_df(mode, n_layers)

up = st.file_uploader("Optional: upload CSV layer table", type=["csv"])
if up is not None:
    try:
        df = pd.read_csv(up)
        if "layer" not in df.columns:
            df.insert(0, "layer", range(1, len(df) + 1))
        st.success("CSV loaded.")
    except Exception as e:
        st.error(f"CSV read failed: {e}")
        df = df0
else:
    df = df0

df = st.data_editor(df, use_container_width=True, hide_index=True, num_rows="fixed" if mode == "Paper 2025 optimized stack" else "dynamic")

st.download_button(
    "Download CSV template",
    data=default_df(mode, n_layers).to_csv(index=False).encode("utf-8-sig"),
    file_name="layer_template_v2.csv",
    mime="text/csv",
)

if st.button("Calculate", type="primary"):
    try:
        if f2 <= f1:
            raise ValueError("Stop frequency must be greater than start frequency.")
        for c in required_cols(mode):
            if c not in df.columns:
                raise ValueError(f"Missing column: {c}")
            if df[c].isna().any():
                raise ValueError(f"Column has empty values: {c}")

        freq = np.linspace(float(f1), float(f2), int(nf))
        angles = parse_angles(angle_text)
        for a in angles:
            if a < 0 or a >= 89.9:
                raise ValueError("Angles must be in [0, 89.9) degrees.")
        pols = ["TE", "TM"] if pol_choice == "Both" else [pol_choice]

        all_rows = []
        summary = []
        figA, axA = plt.subplots()
        figRL, axRL = plt.subplots()

        for pol in pols:
            for angle in angles:
                zin, gamma, R, A, RL = response(freq, df, mode, angle, pol)
                label = f"{pol}, {angle:g} deg"
                axA.plot(freq, A, label=label)
                axRL.plot(freq, RL, label=label)

                idxA = int(np.argmax(A))
                idxR = int(np.argmin(R))
                summary.append({
                    "polarization": pol,
                    "angle_deg": angle,
                    "max_absorptivity_A": A[idxA],
                    "freq_at_max_A_GHz": freq[idxA],
                    "min_reflectivity_R": R[idxR],
                    "freq_at_min_R_GHz": freq[idxR],
                    "max_return_loss_dB": RL[idxR],
                })

                for i in range(len(freq)):
                    all_rows.append({
                        "frequency_GHz": freq[i],
                        "polarization": pol,
                        "angle_deg": angle,
                        "absorptivity_A": A[i],
                        "reflectivity_R": R[i],
                        "return_loss_dB": RL[i],
                        "abs_gamma": abs(gamma[i]),
                        "zin_real_ohm": zin[i].real,
                        "zin_imag_ohm": zin[i].imag,
                    })

        axA.set_xlabel("Frequency (GHz)")
        axA.set_ylabel("Absorptivity")
        axA.set_title("Absorptivity of Metal-backed Multilayer Absorber")
        axA.set_ylim(0, 1.05)
        axA.grid(True)
        axA.legend()

        axRL.set_xlabel("Frequency (GHz)")
        axRL.set_ylabel("S11 Magnitude (dB)")
        axRL.set_title("S11 Magnitude")
        axRL.grid(True)
        axRL.legend()

        summary_df = pd.DataFrame(summary)
        result_df = pd.DataFrame(all_rows)

        st.subheader("Summary")
        st.dataframe(summary_df, use_container_width=True, hide_index=True)

        st.subheader("Absorptivity")
        st.pyplot(figA)

        st.subheader("Return Loss")
        st.pyplot(figRL)

        st.subheader("Full result table")
        st.dataframe(result_df, use_container_width=True, hide_index=True)

        st.download_button("Download full result CSV", result_df.to_csv(index=False).encode("utf-8-sig"), "absorber_result_v2.csv", "text/csv")
        st.download_button("Download summary CSV", summary_df.to_csv(index=False).encode("utf-8-sig"), "absorber_summary_v2.csv", "text/csv")

    except Exception as e:
        st.error(f"Calculation failed: {e}")

st.divider()
st.markdown("""
**CST comparison:** for a metal-backed one-port model, compare with `1 - abs(S11)^2`.
If transmission exists, use `1 - abs(S11)^2 - abs(S21)^2`.
Make sure frequency range, layer order, thickness, material tensor direction, incident angle and polarization are identical.
""")
