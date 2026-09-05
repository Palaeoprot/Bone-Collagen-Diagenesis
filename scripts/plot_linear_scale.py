import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
cohort_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(cohort_path)

collagen = df[df["material_category"] == "COLLAGEN"].copy()
controls = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

# Generate two side-by-side plots with normal (linear) y-axis
fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), dpi=300)

col_color = "#980E1E"    # Deep red / collagen (Palaeoprot palette)
ctrl_color = "#507351"   # Forest green / organic controls (Palaeoprot palette)

# Common limits on linear scale
t_min, t_max = -15, 30
th_max_linear = 220000  # Linear upper limit showing collagen envelope and contrast

# Panel A: Purified Bone Collagen
ax1 = axes[0]
sc1 = ax1.scatter(collagen["integrated_temp_c"], collagen["thermal_age_173"],
                  c=collagen["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar1 = plt.colorbar(sc1, ax=ax1, pad=0.02)
cbar1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

# 95th percentile thermal age ceiling line (~26,278 y)
ax1.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
            label="95th % Thermal Ceiling ($t_{\\mathrm{eff}} = 26,278\\ \\mathrm{y}$ @ 10 °C)")

ax1.set_xlim(t_min, t_max)
ax1.set_ylim(0, th_max_linear)
ax1.set_title(r"$\mathbf{A}$  Purified Bone Collagen ($N = 18,101$)", fontsize=13, fontweight="bold", pad=12, loc="left")
ax1.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Linear Scale)", fontsize=11, fontweight="bold")
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

# Panel B: Non-Collagen Controls (Charcoal, Seed, Wood)
ax2 = axes[1]
sc2 = ax2.scatter(controls["integrated_temp_c"], controls["thermal_age_173"],
                  c=controls["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar2 = plt.colorbar(sc2, ax=ax2, pad=0.02)
cbar2.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

# Reference line from collagen panel for direct visual comparison
ax2.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
            label="Collagen 95% Thermal Ceiling (Reference)")

ax2.set_xlim(t_min, t_max)
ax2.set_ylim(0, th_max_linear)
ax2.set_title(r"$\mathbf{B}$  Non-Collagen Controls: Charcoal, Wood, Seeds ($N = 18,101$)", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax2.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C Reference — Linear Scale)", fontsize=11, fontweight="bold")
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

# Annotation pointing out points in controls escaping beyond the top of the linear frame
n_escaped = (controls["thermal_age_173"] > th_max_linear).sum()
ax2.text(22, th_max_linear * 0.92, f"▲ {n_escaped:,} control samples\nextend above 220,000 y\n(up to 763,744 y)",
         fontsize=9, color="#507351", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0fdf4", edgecolor="#507351", alpha=0.9))

plt.tight_layout()
out_plot = os.path.join(output_dir, "Figure_Integrated_Temp_vs_Thermal_Age_Linear_Scale.png")
plt.savefig(out_plot, dpi=300)
print(f"Saved linear scale figure to: {out_plot}")
