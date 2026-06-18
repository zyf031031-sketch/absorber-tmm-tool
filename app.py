import io
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

C0 = 299_792_458.0
Z0 = 376.730313668

st.set_page_config(page_title="Multilayer Absorber TMM Tool", layout="wide")

st.title("Multilayer Absorber Absorptivity Tool")
st.caption("Normal-incidence, metal-backed multilayer absorber calculator using transmission-line / TMM theory.")

with st.expander("使用说明", expanded=True):
    st.markdown(
        """
        **第一版适用范围：**
        - 多层均匀平板材料；
        - 垂直入射；
        - 最后一面为金属背板，因此透射率近似为 0；
        - 每层材料用复介电常数和复磁导率描述。

        **层顺序：**
        表格第 1 行是最靠近空气入射侧的层，最后 1 行是最靠近金属背板的层。

        **损耗参数约定：**
        程序内部使用：
        \[
        \epsilon_r = \epsilon' - j\epsilon'', \quad \mu_r = \mu' - j\mu''
        \]
        所以表格里的 `eps_loss` 和 `mu_loss` 填正数即可。
        """
    )

st.sidebar.header("Frequency Settings")
f_start = st.sidebar.number_input("Start frequency / GHz", min_value=0.001, value=2.0, step=0.5)
f_stop = st.sidebar.number_input("Stop frequency / GHz", min_value=0.001, value=18.0, step=0.5)
n_points = st.sidebar.number_input("Number of points", min_value=11, max_value=5001, value=401, step=10)

st.sidebar.header("Layer Settings")
n_layers = st.sidebar.number_input("Number of layers", min_value=1, max_value=20, value=2, step=1)

default_rows = []
for i in range(int(n_layers)):
    default_rows.append(
        {
            "layer": i + 1,
            "thickness_mm": 1.0,
            "eps_real": 4.0,
            "eps_loss": 0.3,
            "mu_real": 1.0,
            "mu_loss": 0.0,
        }
    )

st.subheader("Layer Parameters")
st.write("Edit the table below. Row 1 is the air side; the last row is the metal-backed side.")

df = st.data_editor(
    pd.DataFrame(default_rows),
    num_rows="fixed",
    use_container_width=True,
    hide_index=True,
)

uploaded = st.file_uploader(
    "Optional: upload CSV layer table",
    type=["csv"],
    help="CSV columns: thickness_mm, eps_real, eps_loss, mu_real, mu_loss. The optional layer column is allowed.",
)

if uploaded is not None:
    try:
        csv_df = pd.read_csv(uploaded)
        required = ["thickness_mm", "eps_real", "eps_loss", "mu_real", "mu_loss"]
        missing = [c for c in required if c not in csv_df.columns]
        if missing:
            st.error(f"CSV missing required columns: {missing}")
        else:
            if "layer" not in csv_df.columns:
                csv_df.insert(0, "layer", range(1, len(csv_df) + 1))
            df = csv_df[["layer", "thickness_mm", "eps_real", "eps_loss", "mu_real", "mu_loss"]]
            st.success("CSV loaded successfully.")
            st.dataframe(df, use_container_width=True, hide_index=True)
    except Exception as exc:
        st.error(f"Failed to read CSV: {exc}")

def metal_backed_absorption(freq_ghz: np.ndarray, layers: pd.DataFrame):
    """
    Calculates reflection and absorption of a metal-backed multilayer absorber
    at normal incidence.

    Uses transmission-line recursion:
    Zin = Zi * (ZL + j Zi tan(ki d)) / (Zi + j ZL tan(ki d))
    Metal backing is modeled as a short circuit: ZL = 0.
    """
    f_hz = np.asarray(freq_ghz, dtype=float) * 1e9
    k0 = 2 * np.pi * f_hz / C0

    # Metal backing: short-circuit load.
    zin_load = np.zeros_like(f_hz, dtype=complex)

    # Recursion from metal side to air side.
    for _, row in layers.iloc[::-1].iterrows():
        d = float(row["thickness_mm"]) * 1e-3
        eps_r = complex(float(row["eps_real"]), -float(row["eps_loss"]))
        mu_r = complex(float(row["mu_real"]), -float(row["mu_loss"]))

        zi = Z0 * np.sqrt(mu_r / eps_r)
        ki = k0 * np.sqrt(mu_r * eps_r)
        tan_term = np.tan(ki * d)

        numerator = zin_load + 1j * zi * tan_term
        denominator = zi + 1j * zin_load * tan_term
        zin = zi * numerator / denominator
        zin_load = zin

    gamma = (zin_load - Z0) / (zin_load + Z0)
    reflectivity = np.abs(gamma) ** 2
    absorptivity = 1 - reflectivity

    return zin_load, gamma, reflectivity.real, np.clip(absorptivity.real, 0, 1)

run = st.button("Calculate", type="primary")

if run:
    try:
        required_cols = ["thickness_mm", "eps_real", "eps_loss", "mu_real", "mu_loss"]
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing column: {col}")
            if df[col].isna().any():
                raise ValueError(f"Column has empty values: {col}")

        freq = np.linspace(float(f_start), float(f_stop), int(n_points))
        zin, gamma, R, A = metal_backed_absorption(freq, df)

        result = pd.DataFrame(
            {
                "frequency_GHz": freq,
                "reflectivity_R": R,
                "absorptivity_A": A,
                "return_loss_dB": np.where(
                    np.abs(-10 * np.log10(np.maximum(R, 1e-15))) < 1e-10,
                    0,
                    -10 * np.log10(np.maximum(R, 1e-15))
                ),
                "reflection_coefficient_abs": np.abs(gamma),
                "input_impedance_real_ohm": zin.real,
                "input_impedance_imag_ohm": zin.imag,
            }
        )

        c1, c2, c3 = st.columns(3)
        c1.metric("Max absorptivity", f"{A.max():.4f}")
        c2.metric("Frequency at max A", f"{freq[A.argmax()]:.4f} GHz")
        c3.metric("Min reflectivity", f"{R.min():.4e}")

        fig, ax = plt.subplots()
        ax.plot(freq, A, label="Absorptivity A")
        ax.plot(freq, R, label="Reflectivity R")
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("Value")
        ax.set_title("Metal-backed Multilayer Absorber Response")
        ax.set_ylim(0, 1.05)
        ax.grid(True)
        ax.legend()
        st.pyplot(fig)

        fig2, ax2 = plt.subplots()
        ax2.plot(freq, result["return_loss_dB"])
        ax2.set_xlabel("Frequency (GHz)")
        ax2.set_ylabel("Return Loss (dB)")
        ax2.set_title("Return Loss")
        ax2.grid(True)
        st.pyplot(fig2)

        st.subheader("Result Table")
        st.dataframe(result, use_container_width=True)

        csv_buffer = io.StringIO()
        result.to_csv(csv_buffer, index=False)
        st.download_button(
            "Download result CSV",
            data=csv_buffer.getvalue(),
            file_name="absorber_result.csv",
            mime="text/csv",
        )

    except Exception as exc:
        st.error(f"Calculation failed: {exc}")

st.divider()
st.markdown(
    """
    **下一版可扩展功能：**
    - 斜入射；
    - TE / TM 极化；
    - 无金属背板时同时计算 R、T、A；
    - 频率相关材料参数；
    - 自动优化层厚度或材料参数。
    """
)
