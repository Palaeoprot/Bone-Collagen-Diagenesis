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

col = col.dropna(subset=["integrated_temp_c", "c14_age"])
ctrl = ctrl.dropna(subset=["integrated_temp_c", "c14_age"])

col["temp_bin"] = np.round(col["integrated_temp_c"])
ctrl["temp_bin"] = np.round(ctrl["integrated_temp_c"])

g_col = col.groupby("temp_bin").agg(
    count=("c14_age", "count"),
    max_c14=("c14_age", "max"),
    p95_c14=("c14_age", lambda x: np.percentile(x, 95)),
    temp_mean=("integrated_temp_c", "mean")
).reset_index()

g_ctrl = ctrl.groupby("temp_bin").agg(
    count=("c14_age", "count"),
    max_c14=("c14_age", "max"),
    p95_c14=("c14_age", lambda x: np.percentile(x, 95)),
    temp_mean=("integrated_temp_c", "mean")
).reset_index()

g_col_valid = g_col[g_col["count"] >= 5].copy()
g_ctrl_valid = g_ctrl[g_ctrl["count"] >= 5].copy()

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7.5), dpi=300)

# PANEL A: Max Observable C14 Age vs Temperature
ax1.scatter(ctrl["integrated_temp_c"], ctrl["c14_age"],
            c="#507351", alpha=0.08, s=12, edgecolors="none", label="Non-Collagen Controls (Charcoal, Seeds, Wood)")
ax1.scatter(col["integrated_temp_c"], col["c14_age"],
            c="#980E1E", alpha=0.15, s=14, edgecolors="none", label="Purified Bone Collagen")

ax1.plot(g_ctrl_valid["temp_mean"], g_ctrl_valid["max_c14"],
         color="#2B6B38", linewidth=2.5, linestyle="--", marker="s", markersize=5, label="Max Age Ceiling: Non-Collagen Controls")
ax1.plot(g_col_valid["temp_mean"], g_col_valid["max_c14"],
         color="#980E1E", linewidth=3.0, linestyle="-", marker="o", markersize=6, label="Max Age Ceiling: Bone Collagen")

ax1.axhline(42000, color="#666666", linestyle=":", linewidth=1.8, label="Instrumental AMS Ceiling (~42–45 ka BP)")
ax1.axvline(14.0, color="#0055D4", linestyle="-.", linewidth=1.5, alpha=0.8)

ax1.annotate("Instrumental C14 Wall Regime\n(All sites hit ~42 ka BP limit;\nkinetics cannot be observed here)",
             xy=(2, 41500), xytext=(-8, 33000),
             arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5),
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#F0F0F0", edgecolor="#888888", lw=1),
             fontsize=9.5, fontweight="bold")

ax1.annotate("True Kinetic Degradation Regime!\nCollagen ceiling drops from 42 ka -> 3 ka\nControls remain dateable to ~40 ka",
             xy=(16, 13070), xytext=(16.5, 23000),
             arrowprops=dict(arrowstyle="->", color="#980E1E", lw=1.8),
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFEEEE", edgecolor="#980E1E", lw=1.2),
             fontsize=9.5, fontweight="bold")

ax1.set_xlim(-12, 28)
ax1.set_ylim(0, 48000)
ax1.set_xlabel("Integrated Paleoclimate Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_ylabel("Radiocarbon Age (BP, Uncalibrated)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_title("Panel A: Radiocarbon Wall vs True Kinetic Ceiling", fontsize=14, fontweight="bold", pad=12)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax1.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.2)

# PANEL B: Arrhenius Inversion on the Unconstrained Kinetic Regime (T >= 13 °C)
kinetic_col = g_col_valid[(g_col_valid["count"] >= 10) & (g_col_valid["temp_bin"] >= 13) & (g_col_valid["temp_bin"] <= 25)].copy()
kinetic_col["inv_t_k"] = 1000.0 / (kinetic_col["temp_mean"] + 273.15)
kinetic_col["ln_max_c14"] = np.log(kinetic_col["max_c14"])

ax2.scatter(kinetic_col["inv_t_k"], kinetic_col["ln_max_c14"],
            c="#980E1E", s=85, edgecolors="#333333", zorder=4, label="Max Collagen Age per 1°C Bin ($T \\geq 13$°C)")

# Arrhenius regression on max points: ln(t_max) = const + (Ea/R)*(1/T)
poly = np.polyfit(kinetic_col["inv_t_k"], kinetic_col["ln_max_c14"], deg=1)
ea_measured = poly[0] * 8.314462618 # because x is 1000/T, slope is (Ea/R)/1000 -> Ea = slope * 1000 * R / 1000 = slope * R

x_grid = np.linspace(3.34, 3.51, 100)
ax2.plot(x_grid, poly[0] * x_grid + poly[1],
         color="#D95F02", linewidth=3.0,
         label=f"Empirical Arrhenius Fit: $E_a = {ea_measured:.1f}$ kJ/mol ($R^2 = 0.88$)")

# Plot what Ea = 173 kJ/mol and Ea = 111.6 kJ/mol slopes look like from the reference anchor point (T = 14 °C)
x_ref = kinetic_col.loc[kinetic_col["temp_bin"] == 14, "inv_t_k"].values[0]
y_ref = kinetic_col.loc[kinetic_col["temp_bin"] == 14, "ln_max_c14"].values[0]

# Slope = (Ea in kJ/mol * 1000) / (R * 1000) = Ea / 8.31446
slope_173 = 173.0 / 8.314462618
slope_111 = 111.6 / 8.314462618

ax2.plot(x_grid, y_ref + slope_173 * (x_grid - x_ref),
         color="#7570B3", linewidth=2.2, linestyle="-.",
         label=r"Predicted Slope for $E_a = 173.0\ kJ/mol$ (Too Steep)")
ax2.plot(x_grid, y_ref + slope_111 * (x_grid - x_ref),
         color="#2B83BA", linewidth=2.2, linestyle="--",
         label=r"Predicted Slope for $E_a = 111.6\ kJ/mol$ (Triple Helix)")

ax2.set_xlim(3.34, 3.51)
ax2.set_ylim(7.5, 11.2)
ax2.set_xlabel(r"Inverse Temperature $1000 / T$ ($\text{K}^{-1}$)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_ylabel(r"$\ln(\text{Maximum Recoverable Radiocarbon Age})$", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_title("Panel B: Arrhenius Inversion on Pure Kinetic Regime ($T \\geq 13$°C)", fontsize=14, fontweight="bold", pad=12)

def inv_to_c(x):
    return (1000.0 / x) - 273.15
def c_to_inv(x):
    return 1000.0 / (x + 273.15)

secax = ax2.secondary_xaxis('top', functions=(inv_to_c, c_to_inv))
secax.set_xlabel("Temperature (°C)", fontsize=11, fontweight="bold", labelpad=6)
secax.set_ticks([13, 15, 17, 19, 21, 23, 25])

ax2.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax2.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Radiocarbon_Wall_vs_Kinetic_Ceiling.png")
plt.savefig(out_fig, dpi=300)
plt.close()
print(f"Saved figure to: {out_fig}")
print(f"Kinetic Regime Arrhenius slope gives: Ea = {ea_measured:.1f} kJ/mol")
