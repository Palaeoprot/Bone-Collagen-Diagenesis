import os
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
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

# 1. Empirical best-fit 95% curve (Ea = 131.3 kJ/mol)
mod_emp_95 = smf.quantreg("log_thermal_age ~ integrated_temp_c", col).fit(q=0.95)
b_emp = mod_emp_95.params["integrated_temp_c"]
a0_emp = np.exp(mod_emp_95.params["Intercept"])
ea_emp = b_emp * r_gas * (t_ref_k**2)

# 2. Constrained Ea = 173 kJ/mol curve (95% containment)
b_173 = 173.0 / (r_gas * (t_ref_k**2)) # 0.25952
residuals_173 = col["log_thermal_age"] - b_173 * col["integrated_temp_c"]
ln_a_173 = np.percentile(residuals_173, 95)
a0_173 = np.exp(ln_a_173)

# 3. Constrained Ea = 111.6 kJ/mol (Buhr & Gräter 2026 Triple Helix QM/MM)
b_111 = 111.6 / (r_gas * (t_ref_k**2)) # 0.16738
residuals_111 = col["log_thermal_age"] - b_111 * col["integrated_temp_c"]
ln_a_111 = np.percentile(residuals_111, 95)
a0_111 = np.exp(ln_a_111)

print(f"Empirical 95% fit: b = {b_emp:.4f} -> Ea = {ea_emp:.1f} kJ/mol | A0 = {a0_emp:.1f} y")
print(f"Fixed Ea = 173 kJ/mol fit: b = {b_173:.4f} -> Ea = 173.0 kJ/mol | A0 = {a0_173:.1f} y")
print(f"Fixed Ea = 111.6 kJ/mol fit: b = {b_111:.4f} -> Ea = 111.6 kJ/mol | A0 = {a0_111:.1f} y")

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), dpi=300)

temp_grid = np.linspace(-15, 25, 250)
curve_emp = a0_emp * np.exp(b_emp * temp_grid)
curve_173 = a0_173 * np.exp(b_173 * temp_grid)
curve_111 = a0_111 * np.exp(b_111 * temp_grid)

# --- PANEL A: Purified Bone Collagen ---
ax1 = axes[0]
sc1 = ax1.scatter(col["integrated_temp_c"], col["thermal_age_173"],
                  c=col["c14_age"], cmap="viridis", alpha=0.30, s=12, rasterized=True)
cbar1 = plt.colorbar(sc1, ax=ax1, pad=0.02)
cbar1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

# Plot curves
ax1.plot(temp_grid, curve_emp, color="#C41E23", lw=3.0, linestyle="-",
         label=f"Empirical Best Fit 95%: $E_a = {ea_emp:.1f}\\ \\mathrm{{kJ/mol}}$ ($b = {b_emp:.4f}$)")

ax1.plot(temp_grid, curve_173, color="#000000", lw=2.5, linestyle="--",
         label=f"Fixed Lab Estimate: $E_a = 173.0\\ \\mathrm{{kJ/mol}}$ ($b = {b_173:.4f}$ — Too steep!)")

ax1.plot(temp_grid, curve_111, color="#009688", lw=2.2, linestyle="-.",
         label=f"QM/MM Triple Helix: $E_a = 111.6\\ \\mathrm{{kJ/mol}}$ ($b = {b_111:.4f}$)")

# Reference horizontal line
ax1.axhline(26278, color="#888888", linestyle=":", lw=1.5, label="Empirical 10 °C Thermal Ceiling ($26,278\\ \\mathrm{y}$)")

ax1.set_xlim(-15, 25)
ax1.set_ylim(0, 160000)
ax1.set_title(r"$\mathbf{A}$  Purified Bone Collagen: Empirical Fit vs Overestimated $E_a = 173\ \mathrm{kJ\cdot mol^{-1}}$", 
             fontsize=12.5, fontweight="bold", pad=12, loc="left")
ax1.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Linear Scale)", fontsize=11, fontweight="bold")
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.2, loc="upper left")

# --- PANEL B: Logarithmic Comparison Highlighting Slope Discrepancy ---
ax2 = axes[1]
sc2 = ax2.scatter(col["integrated_temp_c"], col["thermal_age_173"],
                  c=col["c14_age"], cmap="viridis", alpha=0.30, s=12, rasterized=True)
cbar2 = plt.colorbar(sc2, ax=ax2, pad=0.02)
cbar2.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

ax2.plot(temp_grid, curve_emp, color="#C41E23", lw=3.0, linestyle="-",
         label=f"Empirical Best Fit: $E_a = {ea_emp:.1f}\\ \\mathrm{{kJ/mol}}$")
ax2.plot(temp_grid, curve_173, color="#000000", lw=2.5, linestyle="--",
         label=f"Fixed Lab Estimate: $E_a = 173.0\\ \\mathrm{{kJ/mol}}$ (Too steep)")
ax2.plot(temp_grid, curve_111, color="#009688", lw=2.2, linestyle="-.",
         label=f"QM/MM Triple Helix: $E_a = 111.6\\ \\mathrm{{kJ/mol}}$")

ax2.set_yscale("log")
ax2.set_xlim(-15, 25)
ax2.set_ylim(10, 500000)
ax2.set_title(r"$\mathbf{B}$  Log-Scale Trajectory Comparison: Demonstrating Why $173\ \mathrm{kJ\cdot mol^{-1}}$ Over-predicts", 
             fontsize=12.5, fontweight="bold", pad=12, loc="left")
ax2.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Log Scale)", fontsize=11, fontweight="bold")
ax2.grid(True, which="both", linestyle=":", alpha=0.5)
ax2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.2, loc="upper left")

# Annotation illustrating mismatch
ax2.annotate("Ea = 173 kJ/mol explodes prematurely\nat temperatures > 15 °C", 
             xy=(18, 300000), xytext=(8, 180000),
             arrowprops=dict(arrowstyle="->", color="black", lw=1.5),
             fontsize=9, fontweight="bold", bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff3cd", edgecolor="#856404"))

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Comparison_Ea_173_vs_Empirical_Fit.png")
plt.savefig(out_fig, dpi=300)
print(f"Saved comparison figure to: {out_fig}")
