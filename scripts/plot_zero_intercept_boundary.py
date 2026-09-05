import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import os

output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
data_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(data_path)

col = df[df["material_category"] == "COLLAGEN"]
ctrl = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"]

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, ax = plt.subplots(figsize=(10, 7), dpi=300)

col_color = "#980E1E"    # Deep red / collagen
ctrl_color = "#507351"   # Forest green / organic controls
line_col = "#C41E23"

ax.scatter(ctrl["c14_age"], ctrl["thermal_age_173"], 
           color=ctrl_color, alpha=0.12, s=8, label=f"Non-Collagen Controls (N = {len(ctrl):,})", rasterized=True)
ax.scatter(col["c14_age"], col["thermal_age_173"], 
           color=col_color, alpha=0.25, s=10, label=f"Purified Bone Collagen (N = {len(col):,})", rasterized=True)

# 95th percentile forced through (0,0)
# Slope = 2.1255
age_grid = np.linspace(0, 42000, 100)
slope_95 = 2.1255
th_line = slope_95 * age_grid

ax.plot(age_grid, th_line, color=line_col, lw=3.0, linestyle="-", 
        label=f"Zero-Intercept 95th % Envelope (Slope = {slope_95:.2f})\n[1% Residual Collagen Horizon: ~26,000–30,000 y @ 10°C]")

# Add reference horizontal line at 95th percentile thermal age (26,278 y)
ax.axhline(26278, color="#111111", linestyle=":", lw=1.8, label="Empirical 95% Thermal Horizon ($t_{\\mathrm{eff}} = 26,278\\ \\mathrm{y}$ @ 10 °C)")

ax.set_title(r"Zero-Intercept Envelope Fitting: 1% Residual Collagen Boundary ($E_a = 173\ \mathrm{kJ\cdot mol^{-1}}$)", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax.set_xlabel(r"Calendar Radiocarbon Age ($^{14}\mathrm{C}$ BP)", fontsize=11, fontweight="bold")
ax.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference)", fontsize=11, fontweight="bold")
ax.set_xlim(0, 43000)
ax.set_ylim(0, 120000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

plt.tight_layout()
out_png = os.path.join(output_dir, "Figure_Zero_Intercept_1pct_Collagen_Boundary.png")
plt.savefig(out_png, dpi=300)
print(f"Plot saved to: {out_png}")
