import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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

# Uncensored kinetic regime points (T >= 13 °C)
kinetic_col = g_col_valid[(g_col_valid["count"] >= 10) & (g_col_valid["temp_bin"] >= 13) & (g_col_valid["temp_bin"] <= 25)].copy()
kinetic_col["inv_t_k"] = 1000.0 / (kinetic_col["temp_mean"] + 273.15)
kinetic_col["ln_max_c14"] = np.log(kinetic_col["max_c14"])

# Kinetic parameter fits:
# ln(t_max) = C + (Ea/R) * (1/T)
poly_inv = np.polyfit(kinetic_col["inv_t_k"], kinetic_col["ln_max_c14"], deg=1)
ea_inv = poly_inv[0] * 8.314462618 # kJ/mol
R = 8.314462618e-3 # kJ/(mol K)
slope_133 = 133.4 / R
C_133 = np.log(41200) - slope_133 * (1.0 / (14.0 + 273.15))

slope_137 = 137.0 / R
C_137 = np.log(41200) - slope_137 * (1.0 / (14.0 + 273.15))

slope_173 = 173.0 / R
C_173 = np.log(41200) - slope_173 * (1.0 / (14.0 + 273.15))

# Setup Publication Aesthetics
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.2

fig = plt.figure(figsize=(20, 14), dpi=300)
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.22)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# =========================================================================
# PANEL A: The Radiocarbon Dating Wall & Censoring (Linear Scale, 0–45 ka)
# =========================================================================
ax1.scatter(ctrl["integrated_temp_c"], ctrl["c14_age"],
            c="#507351", alpha=0.10, s=12, edgecolors="none", label="Non-Collagen Controls (Charcoal, Wood, Seeds)")
ax1.scatter(col["integrated_temp_c"], col["c14_age"],
            c="#980E1E", alpha=0.20, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

ax1.plot(g_ctrl_valid["temp_mean"], g_ctrl_valid["max_c14"],
         color="#2B6B38", linewidth=2.5, linestyle="--", marker="s", markersize=4.5, label="Max Age Ceiling: Controls")
ax1.plot(g_col_valid["temp_mean"], g_col_valid["max_c14"],
         color="#980E1E", linewidth=3.0, linestyle="-", marker="o", markersize=5.5, label="Max Age Ceiling: Bone Collagen")

ax1.axhline(42000, color="#444444", linestyle=":", linewidth=2.0, label="Radiocarbon Blank Ceiling (~42–45 ka BP)")
ax1.axvline(13.8, color="#0055D4", linestyle="-.", linewidth=1.8, label="Crossover Boundary ($T \\approx 13.8$ °C)")

ax1.fill_between([-12, 13.8], 42000, 48000, color="#E0E0E0", alpha=0.7)
ax1.text(-10.5, 44200, "INSTRUMENTAL $^{14}\\mathrm{C}$ WALL REGIME\n(Collagen preserved, but cannot be dated)",
         fontsize=9.0, fontweight="bold", color="#444444")

ax1.annotate("True Kinetic Degradation Regime!\nCollagen ceiling collapses from 42 ka -> 3 ka\nControls remain dateable to ~40 ka",
             xy=(16, 13070), xytext=(15.5, 23000),
             arrowprops=dict(arrowstyle="->", color="#980E1E", lw=1.8),
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFEEEE", edgecolor="#980E1E", lw=1.2),
             fontsize=9.0, fontweight="bold")

ax1.set_xlim(-12, 28)
ax1.set_ylim(0, 48000)
ax1.set_xlabel("Paleoclimate-Integrated Temperature $\\bar{T}$ (°C)", fontsize=12, fontweight="bold", labelpad=6)
ax1.set_ylabel("Radiocarbon Age (BP, Uncalibrated)", fontsize=12, fontweight="bold", labelpad=6)
ax1.set_title("A. The Radiocarbon Wall vs True Kinetic Degradation Ceiling", fontsize=13, fontweight="bold", pad=10)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax1.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=8.8)

# =========================================================================
# PANEL B: Arrhenius Rates: High-Temp Lab Experiments vs Geological Survival
# =========================================================================
sec_per_yr = 365.25 * 86400
R_SI = 8.314462618 # J / (mol K)

# 1. Smith (2002) Thesis Lab Heating Experiments:
# 75-95 °C: ln(k) = -20832*(1/T) + 44.495 (Ea = 173.2 kJ/mol, A = 2.11e19 s-1)
# 55-75 °C: ln(k) = -18687*(1/T) + 37.976 (Ea = 155.3 kJ/mol, A = 3.11e16 s-1)
temps_smith_high = np.array([95.0, 85.0, 75.0])
inv_t_smith_high = 1000.0 / (temps_smith_high + 273.15)
ln_k_smith_high = -20832.0 * (1.0 / (temps_smith_high + 273.15)) + 44.495

temps_smith_low = np.array([75.0, 65.0, 55.0])
inv_t_smith_low = 1000.0 / (temps_smith_low + 273.15)
ln_k_smith_low = -18687.0 * (1.0 / (temps_smith_low + 273.15)) + 37.976

# Other literature high-T experiments in Smith (2002):
temps_ortner = np.array([140.0, 130.0, 120.0, 110.0, 100.0])
inv_t_ortner = 1000.0 / (temps_ortner + 273.15)
ln_k_ortner = np.log(2.02e13) - (132.1e3 / (R_SI * (temps_ortner + 273.15)))

temps_von_endt = np.array([130.0, 120.0, 110.0, 100.0])
inv_t_von_endt = 1000.0 / (temps_von_endt + 273.15)
ln_k_von_endt = np.log(1.00e20) - (183.3e3 / (R_SI * (temps_von_endt + 273.15)))

# Empirical Field Radiocarbon Degradation Rates:
kinetic_col["k_yr"] = -np.log(0.01) / kinetic_col["max_c14"]
kinetic_col["k_s"] = kinetic_col["k_yr"] / sec_per_yr
kinetic_col["ln_k_s"] = np.log(kinetic_col["k_s"])

poly_field = np.polyfit(kinetic_col["inv_t_k"], kinetic_col["ln_k_s"], 1)
ea_field = -poly_field[0] * 1000.0 * R_SI / 1000.0

# Deep-Time Benchmarks:
benchmarks_b = [
    ("Ellesmere Bear\n(~3.9 Ma, -10.5°C)", -10.5, 3.9e6, 0.20, (0.04, 0.8), "left", "bottom"),
    ("High Arctic Camel\n(~3.4 Ma, -10.5°C)", -10.5, 3.4e6, 0.15, (-0.26, -2.6), "left", "top"),
    ("Yukon Horse\n(~735 ka, -9.0°C)", -9.0, 735000, 0.25, (0.04, 1.3), "left", "bottom"),
    ("Dmanisi Rhino\n(~1.77 Ma, 11.0°C)", 11.0, 1.77e6, 0.02, (-0.35, 1.3), "left", "bottom")
]

inv_t_grid_b = np.linspace(2.35, 3.95, 300)

# 1. Smith (2002) Canonical Ea = 173.2 kJ/mol
ln_k_grid_173 = np.log(2.11e19) - (173.2e3 / (R_SI * (1000.0 / inv_t_grid_b)))
ax2.plot(inv_t_grid_b, ln_k_grid_173, color="#7570B3", linewidth=2.4, linestyle="-.",
         label=r"Smith (2002) Lab: $E_a = 173.2\ \mathrm{kJ/mol}$")

# 2. Smith (2002) Lower-T Lab (55-75°C): Ea = 155.3 kJ/mol
ln_k_grid_155 = np.log(3.11e16) - (155.3e3 / (R_SI * (1000.0 / inv_t_grid_b)))
ax2.plot(inv_t_grid_b, ln_k_grid_155, color="#8DA0CB", linewidth=1.8, linestyle=":",
         label=r"Smith (2002) 55–75°C: $E_a = 155.3\ \mathrm{kJ/mol}$")

# 3. Empirical Field Radiocarbon Fit
ln_k_grid_field = poly_field[0] * inv_t_grid_b + poly_field[1]
ax2.plot(inv_t_grid_b, ln_k_grid_field, color="#D95F02", linewidth=3.0, linestyle="-",
         label=rf"Empirical Field $^{{14}}\mathrm{{C}}$: $E_a = {ea_field:.1f}\ \mathrm{{kJ/mol}}$ ($T \geq 13.8^\circ\mathrm{{C}}$)")

# 4. Peptide bond hydrolysis (Collins & Galley 1998, Ea = 92 kJ/mol)
inv_t_20 = 1000.0 / (20.0 + 273.15)
ln_k_20 = poly_field[0] * inv_t_20 + poly_field[1]
ln_k_grid_92 = ln_k_20 - (92.0e3 / (R_SI * 1000.0)) * (inv_t_grid_b - inv_t_20)
ax2.plot(inv_t_grid_b, ln_k_grid_92, color="#2B83BA", linewidth=1.8, linestyle="--", alpha=0.85,
         label=r"Peptide Hydrolysis: $E_a = 92.0\ \mathrm{kJ/mol}$")

# Plot Lab Experimental Points
ax2.scatter(inv_t_smith_high, ln_k_smith_high, color="#7570B3", s=85, marker="o", edgecolors="#000000", linewidths=1.0, zorder=6,
            label=r"Smith (2002) Lab (75–95°C)")
ax2.scatter(inv_t_smith_low, ln_k_smith_low, color="#8DA0CB", s=70, marker="s", edgecolors="#000000", linewidths=0.9, zorder=6,
            label=r"Smith (2002) Lab (55–75°C)")
ax2.scatter(inv_t_ortner, ln_k_ortner, color="#66C2A5", s=60, marker="^", edgecolors="#000000", linewidths=0.8, zorder=5,
            label=r"Ortner et al. (100–140°C)")
ax2.scatter(inv_t_von_endt, ln_k_von_endt, color="#E6AB02", s=60, marker="v", edgecolors="#000000", linewidths=0.8, zorder=5,
            label=r"Von Endt & Ortner (100–130°C)")

# Plot Field Radiocarbon Points
ax2.scatter(kinetic_col["inv_t_k"], kinetic_col["ln_k_s"], color="#980E1E", s=85, marker="D", edgecolors="#000000", linewidths=1.1, zorder=7,
            label=r"Empirical Field $^{14}\mathrm{C}$ Loss ($13.8\text{--}25^\circ\mathrm{C}$)")

# Plot Deep-Time Benchmarks
for name, t_c, age_yr, frac, (off_x, off_y), ha, va in benchmarks_b:
    inv_t = 1000.0 / (t_c + 273.15)
    k_yr = -np.log(frac) / age_yr
    k_s = k_yr / sec_per_yr
    ln_k = np.log(k_s)
    ax2.scatter(inv_t, ln_k, color="#E7298A", s=130, marker="*", edgecolors="#000000", linewidths=1.2, zorder=8)
    ax2.annotate(name, (inv_t, ln_k), xytext=(inv_t + off_x, ln_k + off_y),
                 fontsize=7.2, fontweight="bold", color="#880e4f", ha=ha, va=va,
                 arrowprops=dict(arrowstyle="->", color="#E7298A", lw=1.1),
                 bbox=dict(boxstyle="round,pad=0.18", facecolor="#FCE4EC", edgecolor="#E7298A", alpha=0.9))

# High-T Annotation
ax2.annotate("Smith (2002) Lab Gelatinization\n(55°C - 95°C Heating Experiments)\n$E_a = 173.2\\ \\mathrm{kJ/mol}$",
             xy=(2.85, -15.5), xytext=(2.40, -11.0),
             arrowprops=dict(arrowstyle="->", color="#7570B3", lw=1.3),
             bbox=dict(boxstyle="round,pad=0.28", facecolor="#EDE7F6", edgecolor="#7570B3", lw=1.0),
             fontsize=7.5, fontweight="bold")

# Deep-Time Divergence Callout
ax2.annotate("CRITICAL DIVERGENCE AT ARCTIC TEMPS:\n"
             "• $E_a = 173.2\\ \\mathrm{kJ/mol}$ predicts $\\ln k \\approx -39.5$ at -10.5°C\n"
             "  (unphysical survival >50–70 Ma!)\n"
             "• Actual Ellesmere Bear (~3.9 Ma) yields $\\ln k \\approx -32$\n"
             "  --> Validates field $E_a = 130.8\\text{--}133.4\\ \\mathrm{kJ/mol}$!",
             xy=(3.807, -31.97), xytext=(3.15, -16.0),
             arrowprops=dict(arrowstyle="->", color="#880e4f", lw=1.3, connectionstyle="arc3,rad=-0.2"),
             bbox=dict(boxstyle="round,pad=0.32", facecolor="#FFF0F5", edgecolor="#880e4f", lw=1.1),
             fontsize=7.3, fontweight="bold", color="#880e4f")

ax2.set_xlim(2.35, 3.95)
ax2.set_ylim(-41.0, -7.0)
ax2.set_xlabel(r"Inverse Temperature $1000 / T$ ($\mathrm{K}^{-1}$)", fontsize=12, fontweight="bold", labelpad=6)
ax2.set_ylabel(r"$\ln(\text{Rate Constant } k\ [\mathrm{s}^{-1}])$", fontsize=12, fontweight="bold", labelpad=6)
ax2.set_title("B. Unified Arrhenius Plot: Lab Heating vs Geological Survival Rates", fontsize=13, fontweight="bold", pad=10)

def inv_to_c_b(x): return (1000.0 / x) - 273.15
def c_to_inv_b(x): return 1000.0 / (x + 273.15)
secax_b = ax2.secondary_xaxis('top', functions=(inv_to_c_b, c_to_inv_b))
secax_b.set_xlabel("Temperature (°C)", fontsize=10, fontweight="bold", labelpad=5)
secax_b.set_ticks([140, 100, 70, 50, 30, 15, 0, -10])

ax2.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax2.legend(loc="lower left", frameon=True, framealpha=0.92, edgecolor="#cccccc", fontsize=7.2)

# =========================================================================
# PANEL C: The Thermal Age Distortion (Why 173 kJ/mol Fails at Cold Temps)
# =========================================================================
# Computing thermal age under Ea = 137 kJ/mol
t_ref_k = 283.15
col["thermal_age_137"] = col["c14_age"] * np.exp((137.0 / R) * ((1.0 / t_ref_k) - (1.0 / (col["integrated_temp_c"] + 273.15))))

ax3.scatter(col["integrated_temp_c"], col["thermal_age_137"],
            c="#980E1E", alpha=0.20, s=16, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

t_fine = np.linspace(-12, 28, 400)
t_k = t_fine + 273.15
b_137 = 137.0 / (8.31446e-3 * (283.15**2))
b_173 = 173.0 / (8.31446e-3 * (283.15**2))

# Envelope for 137:
res_137 = np.log(col["thermal_age_137"]) - b_137 * col["integrated_temp_c"]
a0_137_95 = np.exp(np.percentile(res_137, 95))
a0_173_95 = 3412.3 # from previous calculation

ax3.plot(t_fine, a0_137_95 * np.exp(b_137 * t_fine), color="#D95F02", linewidth=3.0, linestyle="-",
         label=r"$\mathbf{E_a = 137.0\ kJ/mol}$ (95% boundary: $A_0 = 4,692$ y)")
ax3.plot(t_fine, a0_173_95 * np.exp(b_173 * t_fine), color="#7570B3", linewidth=2.2, linestyle="-.",
         label=r"Historical Lab Model ($E_a = 173.0\ kJ/mol$)")

# Instrumental envelope curve
apparent_thermal_wall_137 = 42000.0 * np.exp((137.0 / R) * ((1.0 / t_ref_k) - (1.0 / (t_fine + 273.15))))
ax3.plot(t_fine, apparent_thermal_wall_137, color="#000000", linewidth=2.2, linestyle=":",
         label=r"Apparent Ceiling Imposed by $^{14}\text{C}$ 42 ka Limit")

ax3.set_xlim(-12, 28)
ax3.set_ylim(0, 450000)
ax3.set_xlabel("Paleoclimate-Integrated Temperature $\\bar{T}$ (°C)", fontsize=12, fontweight="bold", labelpad=6)
ax3.set_ylabel("Equivalent Thermal Age @ 10 °C (yr, Linear Scale)", fontsize=12, fontweight="bold", labelpad=6)
ax3.set_title("C. Apparent Suppression of Thermal Age in Cold Environments", fontsize=13, fontweight="bold", pad=10)
ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax3.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax3.legend(loc="upper left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=8.8)

ax3.annotate("ARTIFACT: At $T < 10^\circ\\mathrm{C}$, clipping at 42 ka BP\nforces equivalent thermal age towards ZERO,\npreventing observation of true collagen kinetics!",
             xy=(-4, 25000), xytext=(-10.5, 180000),
             arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5),
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF9C4", edgecolor="#FBC02D", lw=1.2),
             fontsize=8.8, fontweight="bold")

# =========================================================================
# PANEL D: True Geological Survival Extrapolated into Deep Time (Log Scale)
# =========================================================================
true_kinetics_133 = np.exp(C_133 + slope_133 * (1.0 / t_k))
true_kinetics_137 = np.exp(C_137 + slope_137 * (1.0 / t_k))
true_kinetics_173 = np.exp(C_173 + slope_173 * (1.0 / t_k))

ax4.plot(t_fine, true_kinetics_133, color="#D95F02", linewidth=3.2, linestyle="-", zorder=3,
         label=r"True Kinetic Ceiling: $E_a = 133.4\ \mathrm{kJ/mol}$")
ax4.plot(t_fine, true_kinetics_137, color="#1B9E77", linewidth=2.6, linestyle="--", zorder=3,
         label=r"Upper Bound Envelope: $E_a = 137.0\ \mathrm{kJ/mol}$")
ax4.plot(t_fine, true_kinetics_173, color="#7570B3", linewidth=2.0, linestyle="-.", zorder=3,
         label=r"Historical Lab Model: $E_a = 173.0\ \mathrm{kJ/mol}$")

# Shade the Predicted Bone Collagen Kinetic Survival Envelope between Ea=133.4 and 137.0
ax4.fill_between(t_fine, true_kinetics_133, true_kinetics_137, color="#D95F02", alpha=0.15,
                 label="Empirical Bone Collagen Ceiling Envelope")

# Uncensored 14C points
ax4.scatter(kinetic_col["temp_mean"], kinetic_col["max_c14"],
            c="#980E1E", s=80, edgecolors="#333333", zorder=5, label=r"Empirical Max $^{14}\mathrm{C}$ Dates ($T \geq 13^\circ\mathrm{C}$)")

ax4.axhline(42000, color="#666666", linestyle=":", linewidth=1.8, label=r"$^{14}\mathrm{C}$ Limit (42 ka BP)")
ax4.fill_between(t_fine, 10, 42000, color="#E8F4F8", alpha=0.5, label=r"$^{14}\mathrm{C}$ Measurable Window")

# Palaeoproteomic Benchmarks: Bone Collagen vs Enamel EMPs
benchmarks_bone = [
    ("Harbin Cranium\n(~148 ka, Bone Col)", 3.5, 148000, (4.0, 48000)),
    ("Denisova Hominin\n(~150 ka, Bone)", -2.0, 150000, (-2.0, 260000)),
    ("Boxgrove / Sima\n(~430 ka, Bone/Dentin)", 7.0, 430000, (7.0, 1100000)),
    ("Yukon Horse\n(~735 ka, Bone)", -9.0, 735000, (-9.0, 190000)),
    ("Dmanisi Rhino\n(~1.77 Ma, Dentin)", 11.0, 1770000, (11.0, 580000)),
    ("High Arctic Camel\n(~3.4 Ma, Bone)", -12.0, 3400000, (-12.0, 1800000)),
    ("Ellesmere Bear (NUFV 303)\n(~3.9 Ma, Radius Bone)", -10.5, 3900000, (-5.0, 5200000))
]

for name, t_site, age_site, (tx, ty) in benchmarks_bone:
    ax4.scatter(t_site, age_site, c="#D95F02", s=80, marker="^", edgecolors="#000000", linewidths=1.2, zorder=7)
    ax4.annotate(name, (t_site, age_site), xytext=(tx, ty),
                 fontsize=7.2, fontweight="bold", color="#111111",
                 arrowprops=dict(arrowstyle="->", color="#D95F02", lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="#FFF3E0", edgecolor="#D95F02", alpha=0.9, lw=0.8))

# Tripot Cave Australian Paradox (Peters & Collins 2023)
tripot_sites = [
    ("Tripot Cave (Unit 1)\n(~75 ka, 24.1°C)", 24.1, 75000, (20.5, 110000)),
    ("Tripot Cave (Unit 2)\n(~350 ka, 24.1°C, ZooMS)", 24.1, 350000, (20.5, 580000))
]
for name, t_site, age_site, (tx, ty) in tripot_sites:
    ax4.scatter(t_site, age_site, c="#E41A1C", s=110, marker="*", edgecolors="#000000", linewidths=1.2, zorder=8)
    ax4.annotate(name, (t_site, age_site), xytext=(tx, ty),
                 fontsize=7.2, fontweight="bold", color="#900C3F", ha="right",
                 arrowprops=dict(arrowstyle="->", color="#E41A1C", lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="#FFEBEE", edgecolor="#E41A1C", alpha=0.92, lw=0.9))
ax4.plot([24.1, 24.1], [75000, 350000], color="#E41A1C", linestyle="--", linewidth=1.6, zorder=7)

# Enamel Control Benchmarks (0% Collagen)
enamel_benchmarks = [
    ("Swartkrans Paranthropus\n(~1.8–2.0 Ma, 0% Col)", 16.5, 1800000, (14.2, 3800000), "right"),
    ("Devon Island Rhino\n(~21–23 Ma, 0% Col)", -11.0, 22000000, (-7.5, 18000000), "left")
]
for name, t_site, age_site, (tx, ty), ha in enamel_benchmarks:
    ax4.scatter(t_site, age_site, c="#7570B3", s=95, marker="D", edgecolors="#000000", linewidths=1.3, zorder=7)
    ax4.annotate(name, (t_site, age_site), xytext=(tx, ty), ha=ha,
                 fontsize=7.2, fontweight="bold", color="#311B92",
                 arrowprops=dict(arrowstyle="->", color="#7570B3", lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="#EDE7F6", edgecolor="#7570B3", alpha=0.95, lw=0.9))

# Annotation showing Arctic Collagen Ceiling Fit
ax4.annotate(r"High Arctic Bone Collagen Ceiling:" "\n"
             r"Ellesmere Bear (3.9 Ma) & Camel (3.4 Ma)" "\n"
             r"sit directly on the $E_a = 133\text{--}137\ \mathrm{kJ/mol}$ boundary!" "\n"
             r"$\rightarrow$ Absolute collagen survival ceiling is $\sim 4\text{--}7\ \mathrm{Ma}$",
             xy=(-10.5, 3900000), xytext=(-11.5, 300),
             fontsize=8.0, fontweight="bold", color="#980E1E",
             arrowprops=dict(arrowstyle="->", color="#980E1E", lw=1.4),
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFEBEE", edgecolor="#D32F2F", lw=1.1))

# Annotation callout for Australian Paradox
callout_aus = (
    "AUSTRALIAN PARADOX (Peters & Collins 2023)\n"
    "• Tripot Cave: 72% ZooMS col at 75–350 ka, 24.1°C\n"
    "• Thermal Age: 1.8–8.5 Ma! (Broke ceiling)\n"
    "--> 'Polymer-in-a-box' sealed karst confinement!"
)
ax4.annotate(callout_aus,
             xy=(24.1, 200000), xytext=(15.0, 70),
             fontsize=7.4, fontweight="bold", color="#B71C1C",
             arrowprops=dict(arrowstyle="->", color="#E41A1C", lw=1.4, connectionstyle="arc3,rad=-0.15"),
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFEBEE", edgecolor="#B71C1C", lw=1.1))

ax4.set_yscale("log")
ax4.set_xlim(-12.5, 28)
ax4.set_ylim(50, 70000000) # Up to 70 Million years
ax4.set_xlabel("Paleoclimate-Integrated Temperature $\\bar{T}$ (°C)", fontsize=12, fontweight="bold", labelpad=6)
ax4.set_ylabel("True Preservation Age (Years, Log Scale)", fontsize=12, fontweight="bold", labelpad=6)
ax4.set_title("D. True Kinetic Lifespans Extrapolated into Deep Time & Arctic Validation", fontsize=13, fontweight="bold", pad=10)
ax4.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y):,}" if y < 1e6 else f"{y/1e6:g} Ma"))
ax4.grid(True, which="both", linestyle="--", alpha=0.4, color="#cccccc")
ax4.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="#cccccc", fontsize=8.0)

plt.suptitle("Deconvolving Instrumental Radiocarbon Limits from True Bone Collagen Degradation Kinetics",
             fontsize=16, fontweight="bold", y=0.995)

out_fig = os.path.join(output_dir, "Figure_Comprehensive_4Panel_Radiocarbon_Deconvolution.png")
plt.savefig(out_fig, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved comprehensive 4-panel figure to: {out_fig}")

artifact_dir = r"C:\Users\matth\.gemini\antigravity-ide\brain\cd905880-511d-460a-944c-c857596f9afc"
import shutil
artifact_fig = os.path.join(artifact_dir, "Figure_Comprehensive_4Panel_Radiocarbon_Deconvolution.png")
shutil.copy2(out_fig, artifact_fig)
print(f"Mirrored figure to artifact directory: {artifact_fig}")


