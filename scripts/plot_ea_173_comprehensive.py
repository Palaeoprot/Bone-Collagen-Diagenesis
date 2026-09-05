import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
cohort_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(cohort_path)

col = df[df["material_category"] == "COLLAGEN"].copy()
ctrl = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()
col["log_thermal_age"] = np.log(col["thermal_age_173"])

r_gas = 8.314462618e-3
t_ref_k = 283.15 # 10 °C
factor = r_gas * (t_ref_k**2)

# Fixed Ea = 173 kJ/mol slope:
b_173 = 173.0 / factor # 0.25952
res_173 = col["log_thermal_age"] - b_173 * col["integrated_temp_c"]

# Calculate A0 across percentiles for 173 kJ/mol:
percentiles = [50, 90, 95, 98, 99, 100]
a0_dict = {}
for p in percentiles:
    ln_a = np.percentile(res_173, p)
    a0_dict[p] = np.exp(ln_a)
    print(f"Ea = 173 kJ/mol, Percentile {p}%: A0 = {a0_dict[p]:.1f} y")

t_grid = np.linspace(-12, 28, 400)

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), dpi=300)

# PANEL A: Linear Scale
ax1.scatter(col["integrated_temp_c"], col["thermal_age_173"],
            c="#980E1E", alpha=0.20, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

ax1.plot(t_grid, a0_dict[50] * np.exp(b_173 * t_grid), color="#4DAF4A", linewidth=2.0, linestyle=":",
         label=f"50% Median ($A_0={a0_dict[50]:.0f}$ y)")
ax1.plot(t_grid, a0_dict[90] * np.exp(b_173 * t_grid), color="#377EB8", linewidth=2.2, linestyle="-.",
         label=f"90% Limit ($A_0={a0_dict[90]:.0f}$ y)")
ax1.plot(t_grid, a0_dict[95] * np.exp(b_173 * t_grid), color="#E41A1C", linewidth=2.8, linestyle="-",
         label=f"95% Upper Limit ($A_0={a0_dict[95]:.0f}$ y)")
ax1.plot(t_grid, a0_dict[98] * np.exp(b_173 * t_grid), color="#984EA3", linewidth=2.2, linestyle="--",
         label=f"98% Upper Limit ($A_0={a0_dict[98]:.0f}$ y)")
ax1.plot(t_grid, a0_dict[99] * np.exp(b_173 * t_grid), color="#FF7F00", linewidth=2.0, linestyle="--",
         label=f"99% Upper Limit ($A_0={a0_dict[99]:.0f}$ y)")
ax1.plot(t_grid, a0_dict[100] * np.exp(b_173 * t_grid), color="#000000", linewidth=2.2, linestyle="-",
         label=f"100% Maximum Envelope ($A_0={a0_dict[100]:.0f}$ y)")

ax1.set_xlim(-12, 28)
ax1.set_ylim(0, 600000)
ax1.set_xlabel("Integrated Paleoclimate Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_ylabel("Equivalent Thermal Age @ 10 °C (yr, Linear Scale)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_title(r"$\mathbf{E_a = 173\ kJ\cdot mol^{-1}}$ Upper Limit Fits (Linear Scale)", fontsize=14, fontweight="bold", pad=12)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax1.legend(loc="upper left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

# PANEL B: Logarithmic Scale
ax2.scatter(col["integrated_temp_c"], col["thermal_age_173"],
            c="#980E1E", alpha=0.20, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

ax2.plot(t_grid, a0_dict[50] * np.exp(b_173 * t_grid), color="#4DAF4A", linewidth=2.0, linestyle=":",
         label=f"50% Median ($A_0={a0_dict[50]:.0f}$ y)")
ax2.plot(t_grid, a0_dict[90] * np.exp(b_173 * t_grid), color="#377EB8", linewidth=2.2, linestyle="-.",
         label=f"90% Limit ($A_0={a0_dict[90]:.0f}$ y)")
ax2.plot(t_grid, a0_dict[95] * np.exp(b_173 * t_grid), color="#E41A1C", linewidth=2.8, linestyle="-",
         label=f"95% Upper Limit ($A_0={a0_dict[95]:.0f}$ y)")
ax2.plot(t_grid, a0_dict[98] * np.exp(b_173 * t_grid), color="#984EA3", linewidth=2.2, linestyle="--",
         label=f"98% Upper Limit ($A_0={a0_dict[98]:.0f}$ y)")
ax2.plot(t_grid, a0_dict[99] * np.exp(b_173 * t_grid), color="#FF7F00", linewidth=2.0, linestyle="--",
         label=f"99% Upper Limit ($A_0={a0_dict[99]:.0f}$ y)")
ax2.plot(t_grid, a0_dict[100] * np.exp(b_173 * t_grid), color="#000000", linewidth=2.2, linestyle="-",
         label=f"100% Maximum Envelope ($A_0={a0_dict[100]:.0f}$ y)")

ax2.set_yscale("log")
ax2.set_xlim(-12, 28)
ax2.set_ylim(20, 10000000)
ax2.set_xlabel("Integrated Paleoclimate Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_ylabel("Equivalent Thermal Age @ 10 °C (yr, Log Scale)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_title(r"$\mathbf{E_a = 173\ kJ\cdot mol^{-1}}$ Upper Limit Fits (Logarithmic Scale)", fontsize=14, fontweight="bold", pad=12)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax2.grid(True, which="both", linestyle="--", alpha=0.4, color="#cccccc")
ax2.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Ea_173_All_Upper_Limits.png")
plt.savefig(out_fig, dpi=300)
plt.close()
print(f"Saved figure to: {out_fig}")
