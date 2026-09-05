import os, time
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Output directory: external data folder per Data Analysis Output Location Directive
output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
os.makedirs(output_dir, exist_ok=True)

data_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
print(f"Loading cohort from {data_path}...")
df = pd.read_parquet(data_path)

# Extract Collagen vs Controls
collagen = df[df["material_category"] == "COLLAGEN"].copy()
controls = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

print(f"Collagen samples: {len(collagen):,}")
print(f"Control samples: {len(controls):,}")

# Quantile Regression at tau = 0.95 (95th percentile upper boundary)
print("\nFitting 95th-percentile Quantile Regression models...")
# 1. Collagen at Ea = 173 kJ/mol
mod_col_173 = smf.quantreg("thermal_age_173 ~ c14_age", collagen).fit(q=0.95)
print("Collagen (173 kJ/mol) 95th % Model:")
print(mod_col_173.summary())

# 2. Controls at Ea = 173 kJ/mol
mod_ctrl_173 = smf.quantreg("thermal_age_173 ~ c14_age", controls).fit(q=0.95)
print("\nControls (173 kJ/mol) 95th % Model:")
print(mod_ctrl_173.summary())

# 3. Collagen at Ea = 100 kJ/mol (Uncatalyzed solution dipeptide)
mod_col_100 = smf.quantreg("thermal_age_100 ~ c14_age", collagen).fit(q=0.95)

# --- Visualization ---
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=300)

# Color scheme from Palaeoprot standardized rules
col_color = "#980E1E"    # Deep red / collagen
ctrl_color = "#507351"   # Forest green / organic controls
line_col = "#C41E23"
line_ctrl = "#30802f"

# Panel A: Ea = 173 kJ/mol (Collins Lab Kinetic Benchmark)
ax1 = axes[0]
ax1.scatter(controls["c14_age"], controls["thermal_age_173"], 
            color=ctrl_color, alpha=0.15, s=8, label="Non-Collagen Controls (Charcoal, Seed, Wood)", rasterized=True)
ax1.scatter(collagen["c14_age"], collagen["thermal_age_173"], 
            color=col_color, alpha=0.25, s=10, label="Purified Bone Collagen (Target Signal)", rasterized=True)

# 95th percentile regression lines
age_grid = np.linspace(500, 42000, 100)
pred_col_173 = mod_col_173.predict({"c14_age": age_grid})
pred_ctrl_173 = mod_ctrl_173.predict({"c14_age": age_grid})

ax1.plot(age_grid, pred_col_173, color=line_col, lw=2.5, linestyle="-", 
         label=f"Collagen 95th % Envelope (Slope = {mod_col_173.params['c14_age']:.2f})")
ax1.plot(age_grid, pred_ctrl_173, color=line_ctrl, lw=2.5, linestyle="--", 
         label=f"Controls 95th % Envelope (Slope = {mod_ctrl_173.params['c14_age']:.2f})")

ax1.set_title(r"$\mathbf{A}$  Thermal Age vs Chronological Age ($E_a = 173\ \mathrm{kJ\cdot mol^{-1}}$)", fontsize=13, pad=12, loc="left")
ax1.set_xlabel("Calendar Radiocarbon Age ($^{14}\mathrm{C}$ BP)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=11, fontweight="bold")
ax1.set_xlim(0, 43000)
ax1.set_ylim(0, 120000)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")

# Panel B: Ea = 100 kJ/mol (Uncatalyzed Solution Dipeptide Benchmark)
ax2 = axes[1]
ax2.scatter(controls["c14_age"], controls["thermal_age_100"], 
            color=ctrl_color, alpha=0.15, s=8, label="Non-Collagen Controls", rasterized=True)
ax2.scatter(collagen["c14_age"], collagen["thermal_age_100"], 
            color=col_color, alpha=0.25, s=10, label="Purified Bone Collagen", rasterized=True)

pred_col_100 = mod_col_100.predict({"c14_age": age_grid})
ax2.plot(age_grid, pred_col_100, color=line_col, lw=2.5, linestyle="-", 
         label=f"Collagen 95th % Envelope (Slope = {mod_col_100.params['c14_age']:.2f})")

ax2.set_title(r"$\mathbf{B}$  Thermal Age vs Chronological Age ($E_a = 100\ \mathrm{kJ\cdot mol^{-1}}$ — Solution)", fontsize=13, pad=12, loc="left")
ax2.set_xlabel("Calendar Radiocarbon Age ($^{14}\mathrm{C}$ BP)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=11, fontweight="bold")
ax2.set_xlim(0, 43000)
ax2.set_ylim(0, 60000)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")

plt.tight_layout()
plot_path = os.path.join(output_dir, "Figure_Thermal_Age_Collagen_vs_Controls_Ea_Fitting.png")
plt.savefig(plot_path, dpi=300)
print(f"\nFigure saved to: {plot_path}")

# Output summary table
summary_stats = {
    "Model": ["Collagen (173 kJ/mol)", "Controls (173 kJ/mol)", "Collagen (100 kJ/mol)"],
    "Slope": [mod_col_173.params["c14_age"], mod_ctrl_173.params["c14_age"], mod_col_100.params["c14_age"]],
    "Intercept": [mod_col_173.params["Intercept"], mod_ctrl_173.params["Intercept"], mod_col_100.params["Intercept"]],
    "p_value_slope": [mod_col_173.pvalues["c14_age"], mod_ctrl_173.pvalues["c14_age"], mod_col_100.pvalues["c14_age"]]
}
summary_df = pd.DataFrame(summary_stats)
summary_path = os.path.join(output_dir, "quantile_regression_summary.csv")
summary_df.to_csv(summary_path, index=False)
print(f"Summary table saved to: {summary_path}")
