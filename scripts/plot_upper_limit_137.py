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
factor = r_gas * (t_ref_k**2)

# Candidate Ea fits:
# 1. Ea = 137.0 kJ/mol (User's observation)
b_137 = 137.0 / factor # 0.20551
res_137 = col["log_thermal_age"] - b_137 * col["integrated_temp_c"]
ln_a0_137_95 = np.percentile(res_137, 95)
a0_137_95 = np.exp(ln_a0_137_95)
# Also 98% upper boundary for 137
ln_a0_137_98 = np.percentile(res_137, 98)
a0_137_98 = np.exp(ln_a0_137_98)

# 2. Empirical 95% Quantile Regression fit (Ea = 131.3 kJ/mol)
mod_95 = smf.quantreg("log_thermal_age ~ integrated_temp_c", col).fit(q=0.95)
b_95 = mod_95.params["integrated_temp_c"]
a0_95 = np.exp(mod_95.params["Intercept"])
ea_95 = b_95 * factor

# 3. Fixed Ea = 173 kJ/mol (Historical lab estimate)
b_173 = 173.0 / factor
res_173 = col["log_thermal_age"] - b_173 * col["integrated_temp_c"]
ln_a0_173 = np.percentile(res_173, 95)
a0_173 = np.exp(ln_a0_173)

# 4. Fixed Ea = 111.6 kJ/mol (Triple Helix QM/MM)
b_111 = 111.6 / factor
res_111 = col["log_thermal_age"] - b_111 * col["integrated_temp_c"]
ln_a0_111 = np.percentile(res_111, 95)
a0_111 = np.exp(ln_a0_111)

print(f"Ea = 137.0 kJ/mol (95% env): b = {b_137:.4f} | A0 = {a0_137_95:.1f} y")
print(f"Ea = 137.0 kJ/mol (98% env): b = {b_137:.4f} | A0 = {a0_137_98:.1f} y")
print(f"Empirical 95% fit: b = {b_95:.4f} -> Ea = {ea_95:.1f} kJ/mol | A0 = {a0_95:.1f} y")
print(f"Ea = 173.0 kJ/mol (95% env): b = {b_173:.4f} | A0 = {a0_173:.1f} y")

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

t_grid = np.linspace(-12, 28, 400)
fit_137_95 = a0_137_95 * np.exp(b_137 * t_grid)
fit_137_98 = a0_137_98 * np.exp(b_137 * t_grid)
fit_emp_95 = a0_95 * np.exp(b_95 * t_grid)
fit_173_95 = a0_173 * np.exp(b_173 * t_grid)
fit_111_95 = a0_111 * np.exp(b_111 * t_grid)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7.5), dpi=300)

# PANEL A: Linear Scale (0 to 600,000 yr)
ax1.scatter(col["integrated_temp_c"], col["thermal_age_173"],
            c="#980E1E", alpha=0.18, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

# Plot curves
ax1.plot(t_grid, fit_137_95, color="#D95F02", linewidth=3.0, linestyle="-",
         label=r"$\mathbf{E_a = 137.0\ kJ/mol}$ (95% boundary: $b=0.2055$, $A_0=4,692$ y)")
ax1.plot(t_grid, fit_137_98, color="#D95F02", linewidth=2.0, linestyle="--", alpha=0.85,
         label=r"$E_a = 137.0\ kJ/mol$ (98% boundary: $A_0=5,648$ y)")
ax1.plot(t_grid, fit_emp_95, color="#1B9E77", linewidth=2.2, linestyle=":",
         label=r"Empirical 95% Fit ($E_a = 131.3\ kJ/mol$, $b=0.1970$, $A_0=4,956$ y)")
ax1.plot(t_grid, fit_173_95, color="#7570B3", linewidth=2.2, linestyle="-.",
         label=r"Lab Estimate ($E_a = 173.0\ kJ/mol$, $b=0.2595$, $A_0=3,412$ y)")
ax1.plot(t_grid, fit_111_95, color="#2B83BA", linewidth=1.8, linestyle="--", alpha=0.75,
         label=r"Triple Helix QM/MM ($E_a = 111.6\ kJ/mol$, $b=0.1674$)")

ax1.set_xlim(-12, 28)
ax1.set_ylim(0, 500000)
ax1.set_xlabel("Integrated Paleoclimate Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_ylabel("Equivalent Thermal Age @ 10 °C (yr, Linear Scale)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_title("Panel A: Upper Limit Fit — Linear Scale", fontsize=14, fontweight="bold", pad=12)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax1.legend(loc="upper left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

# Annotation callout on ax1
ax1.annotate(r"$\mathbf{E_a = 137\ kJ/mol}$ perfectly balances" + "\n" +
             "cold permafrost and warm Mediterranean ceilings\nwithout astronomical runaway ($173\ kJ/mol$)",
             xy=(16, a0_137_95 * np.exp(b_137 * 16)), xytext=(5, 360000),
             arrowprops=dict(arrowstyle="->", color="#D95F02", lw=1.8),
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF3E0", edgecolor="#D95F02", lw=1.2),
             fontsize=9.5, fontweight="bold")

# PANEL B: Logarithmic Scale
ax2.scatter(col["integrated_temp_c"], col["thermal_age_173"],
            c="#980E1E", alpha=0.18, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

ax2.plot(t_grid, fit_137_95, color="#D95F02", linewidth=3.0, linestyle="-",
         label=r"$\mathbf{E_a = 137.0\ kJ/mol}$ (95% envelope)")
ax2.plot(t_grid, fit_137_98, color="#D95F02", linewidth=2.0, linestyle="--", alpha=0.85,
         label=r"$E_a = 137.0\ kJ/mol$ (98% envelope)")
ax2.plot(t_grid, fit_emp_95, color="#1B9E77", linewidth=2.2, linestyle=":",
         label=r"Empirical 95% Fit ($E_a = 131.3\ kJ/mol$)")
ax2.plot(t_grid, fit_173_95, color="#7570B3", linewidth=2.2, linestyle="-.",
         label=r"Lab Estimate ($E_a = 173.0\ kJ/mol$)")
ax2.plot(t_grid, fit_111_95, color="#2B83BA", linewidth=1.8, linestyle="--", alpha=0.75,
         label=r"Triple Helix QM/MM ($E_a = 111.6\ kJ/mol$)")

ax2.set_yscale("log")
ax2.set_xlim(-12, 28)
ax2.set_ylim(50, 3000000)
ax2.set_xlabel("Integrated Paleoclimate Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_ylabel("Equivalent Thermal Age @ 10 °C (yr, Log Scale)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_title("Panel B: Upper Limit Fit — Logarithmic Scale", fontsize=14, fontweight="bold", pad=12)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax2.grid(True, which="both", linestyle="--", alpha=0.4, color="#cccccc")
ax2.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Upper_Limit_Fit_Ea_137.png")
plt.savefig(out_fig, dpi=300)
plt.close()
print(f"Saved figure to: {out_fig}")
