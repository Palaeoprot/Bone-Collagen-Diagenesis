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

# 2x2 grid: Top row zoomed to 0-80,000 years; Bottom row full range (0-800,000 years)
fig, axes = plt.subplots(2, 2, figsize=(16, 13), dpi=300)

t_min, t_max = -15, 30

# --- ROW 1: Focus / Zoom (0 to 80,000 years) ---
# Panel A1: Collagen
ax_a1 = axes[0, 0]
sc_a1 = ax_a1.scatter(collagen["integrated_temp_c"], collagen["thermal_age_173"],
                      c=collagen["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar_a1 = plt.colorbar(sc_a1, ax=ax_a1, pad=0.02)
cbar_a1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=9, fontweight="bold")
cbar_a1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

ax_a1.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
              label="95th % Thermal Ceiling ($t_{\\mathrm{eff}} = 26,278\\ \\mathrm{y}$ @ 10 °C)")
ax_a1.set_xlim(t_min, t_max)
ax_a1.set_ylim(0, 80000)
ax_a1.set_title(r"$\mathbf{A1}$  Purified Bone Collagen (0 – 80,000 y view)", fontsize=12, fontweight="bold", pad=10, loc="left")
ax_a1.set_xlabel("Integrated Paleotemperature (°C)", fontsize=10, fontweight="bold")
ax_a1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=10, fontweight="bold")
ax_a1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax_a1.grid(True, linestyle=":", alpha=0.5)
ax_a1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")

# Panel B1: Controls (0 to 80,000 years)
ax_b1 = axes[0, 1]
sc_b1 = ax_b1.scatter(controls["integrated_temp_c"], controls["thermal_age_173"],
                      c=controls["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar_b1 = plt.colorbar(sc_b1, ax=ax_b1, pad=0.02)
cbar_b1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=9, fontweight="bold")
cbar_b1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

ax_b1.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
              label="Collagen 95% Thermal Ceiling (Reference)")
ax_b1.set_xlim(t_min, t_max)
ax_b1.set_ylim(0, 80000)
ax_b1.set_title(r"$\mathbf{B1}$  Non-Collagen Controls (0 – 80,000 y view)", fontsize=12, fontweight="bold", pad=10, loc="left")
ax_b1.set_xlabel("Integrated Paleotemperature (°C)", fontsize=10, fontweight="bold")
ax_b1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=10, fontweight="bold")
ax_b1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax_b1.grid(True, linestyle=":", alpha=0.5)
ax_b1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")


# --- ROW 2: Full Linear Range (0 to 800,000 years) ---
# Panel A2: Collagen Full
ax_a2 = axes[1, 0]
sc_a2 = ax_a2.scatter(collagen["integrated_temp_c"], collagen["thermal_age_173"],
                      c=collagen["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar_a2 = plt.colorbar(sc_a2, ax=ax_a2, pad=0.02)
cbar_a2.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=9, fontweight="bold")
cbar_a2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

ax_a2.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
              label="Collagen 95% Thermal Ceiling ($26,278\\ \\mathrm{y}$)")
ax_a2.set_xlim(t_min, t_max)
ax_a2.set_ylim(0, 800000)
ax_a2.set_title(r"$\mathbf{A2}$  Purified Bone Collagen (Full Scale 0 – 800,000 y)", fontsize=12, fontweight="bold", pad=10, loc="left")
ax_a2.set_xlabel("Integrated Paleotemperature (°C)", fontsize=10, fontweight="bold")
ax_a2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=10, fontweight="bold")
ax_a2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax_a2.grid(True, linestyle=":", alpha=0.5)
ax_a2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")

# Panel B2: Controls Full
ax_b2 = axes[1, 1]
sc_b2 = ax_b2.scatter(controls["integrated_temp_c"], controls["thermal_age_173"],
                      c=controls["c14_age"], cmap="viridis", alpha=0.35, s=12, rasterized=True)
cbar_b2 = plt.colorbar(sc_b2, ax=ax_b2, pad=0.02)
cbar_b2.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=9, fontweight="bold")
cbar_b2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

ax_b2.axhline(26278, color="#980E1E", linestyle="--", lw=2.2, 
              label="Collagen 95% Thermal Ceiling (Reference)")
ax_b2.set_xlim(t_min, t_max)
ax_b2.set_ylim(0, 800000)
ax_b2.set_title(r"$\mathbf{B2}$  Non-Collagen Controls (Full Scale 0 – 800,000 y)", fontsize=12, fontweight="bold", pad=10, loc="left")
ax_b2.set_xlabel("Integrated Paleotemperature (°C)", fontsize=10, fontweight="bold")
ax_b2.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C)", fontsize=10, fontweight="bold")
ax_b2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax_b2.grid(True, linestyle=":", alpha=0.5)
ax_b2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9, loc="upper left")

plt.tight_layout()
out_quad = os.path.join(output_dir, "Figure_Integrated_Temp_vs_Thermal_Age_Linear_MultiScale.png")
plt.savefig(out_quad, dpi=300)
print(f"Saved multi-scale linear figure to: {out_quad}")
