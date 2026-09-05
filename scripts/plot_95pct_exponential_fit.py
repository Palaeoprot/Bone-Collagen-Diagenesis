import os
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
cohort_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(cohort_path)

collagen = df[df["material_category"] == "COLLAGEN"].copy()
controls = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

# Fit the 95% exponential boundary curve: ln(thermal_age) = ln(A) + b * T
collagen["log_thermal_age"] = np.log(collagen["thermal_age_173"])
mod_exp_95 = smf.quantreg("log_thermal_age ~ integrated_temp_c", collagen).fit(q=0.95)

ln_a = mod_exp_95.params["Intercept"]
b = mod_exp_95.params["integrated_temp_c"]
a_0 = np.exp(ln_a)

# Verify containment
pred_log_95 = mod_exp_95.predict(collagen)
pct_contained = (collagen["log_thermal_age"] <= pred_log_95).mean() * 100.0
print(f"95% Exponential Fit: Thermal_Age_95 = {a_0:.2f} * exp({b:.4f} * T)")
print(f"Percentage of collagen points strictly contained: {pct_contained:.2f}%")

# Setup Plot styling
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), dpi=300)

col_color = "#980E1E"    # Deep red / collagen (Palaeoprot palette)
ctrl_color = "#507351"   # Forest green / organic controls (Palaeoprot palette)
curve_color = "#C41E23"  # Vibrant red for boundary curve

t_min, t_max = -15, 30
th_max_linear = 220000

# Temperature grid for plotting the curve
temp_grid = np.linspace(-15, 25, 200)
curve_95 = a_0 * np.exp(b * temp_grid)

# --- PANEL A: Purified Bone Collagen with 95% Exponential Boundary Curve ---
ax1 = axes[0]
sc1 = ax1.scatter(collagen["integrated_temp_c"], collagen["thermal_age_173"],
                  c=collagen["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar1 = plt.colorbar(sc1, ax=ax1, pad=0.02)
cbar1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

# Plot the 95% Exponential Boundary Curve
ax1.plot(temp_grid, curve_95, color=curve_color, lw=3.0, linestyle="-",
         label=f"95% Exponential Upper Envelope:\n$t_{{\\mathrm{{eff, 95\\%}}}} = {a_0:,.0f} \\cdot e^{{{b:.4f} \\cdot T}}$ ({pct_contained:.1f}% contained)")

# Also include horizontal line at 10°C benchmark for comparison
ax1.axhline(26278, color="#555555", linestyle=":", lw=1.5, 
            label="Empirical 10 °C Thermal Ceiling ($26,278\\ \\mathrm{y}$)")

ax1.set_xlim(t_min, t_max)
ax1.set_ylim(0, th_max_linear)
ax1.set_title(r"$\mathbf{A}$  Purified Bone Collagen ($N = 18,101$) with 95% Exponential Fit", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax1.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Linear Scale)", fontsize=11, fontweight="bold")
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

# --- PANEL B: Non-Collagen Controls for Direct Comparison ---
ax2 = axes[1]
sc2 = ax2.scatter(controls["integrated_temp_c"], controls["thermal_age_173"],
                  c=controls["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar2 = plt.colorbar(sc2, ax=ax2, pad=0.02)
cbar2.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

# Plot the same 95% collagen exponential curve onto Panel B to show how controls break through
ax2.plot(temp_grid, curve_95, color=curve_color, lw=2.5, linestyle="--",
         label="Collagen 95% Exponential Envelope (Reference)")

ax2.set_xlim(t_min, t_max)
ax2.set_ylim(0, th_max_linear)
ax2.set_title(r"$\mathbf{B}$  Non-Collagen Controls: Charcoal, Wood, Seeds ($N = 18,101$)", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax2.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Linear Scale)", fontsize=11, fontweight="bold")
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

n_escaped = (controls["thermal_age_173"] > th_max_linear).sum()
ax2.text(18, th_max_linear * 0.90, f"▲ {n_escaped:,} control samples\nextend above 220,000 y\n(up to 763,744 y)",
         fontsize=9, color="#507351", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0fdf4", edgecolor="#507351", alpha=0.9))

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Panel_A_95pct_Exponential_Fit.png")
plt.savefig(out_fig, dpi=300)
print(f"Saved figure with 95% exponential fit to: {out_fig}")
