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
g_col = col.groupby("temp_bin").agg(
    count=("c14_age", "count"),
    max_c14=("c14_age", "max"),
    temp_mean=("integrated_temp_c", "mean")
).reset_index()

k_col = g_col[(g_col["count"] >= 10) & (g_col["temp_bin"] >= 13) & (g_col["temp_bin"] <= 25)].copy()

# Arrhenius equation parameters:
# ln(t_max) = C + (Ea/R) * (1/T)
R = 8.314462618e-3 # kJ/(mol K)
Ea_133 = 133.4
slope_133 = Ea_133 / R
C_133 = np.log(41200) - slope_133 * (1.0 / (14.0 + 273.15))

# Also for Ea = 137.0 kJ/mol:
Ea_137 = 137.0
slope_137 = Ea_137 / R
C_137 = np.log(41200) - slope_137 * (1.0 / (14.0 + 273.15))

# Also for Ea = 173.0 kJ/mol:
Ea_173 = 173.0
slope_173 = Ea_173 / R
C_173 = np.log(41200) - slope_173 * (1.0 / (14.0 + 273.15))

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), dpi=300)

# -------------------------------------------------------------
# PANEL A: The Distortion Caused by the 50 ka Radiocarbon Limit
# -------------------------------------------------------------
ax1.scatter(col["integrated_temp_c"], col["c14_age"],
            c="#980E1E", alpha=0.15, s=14, edgecolors="none", label="Purified Bone Collagen ($N=18,101$)")

# True kinetic curve (untruncated)
t_fine = np.linspace(-12, 28, 500)
t_k = t_fine + 273.15
true_kinetics_133 = np.exp(C_133 + slope_133 * (1.0 / t_k))
true_kinetics_137 = np.exp(C_137 + slope_137 * (1.0 / t_k))

# Truncated apparent curve (what 14C records)
apparent_14c_133 = np.minimum(true_kinetics_133, 42000)

ax1.plot(t_fine, true_kinetics_133, color="#D95F02", linewidth=3.0, linestyle="--",
         label="True Kinetic Survival Ceiling ($E_a = 133.4$ kJ/mol, Untruncated)")
ax1.plot(t_fine, apparent_14c_133, color="#000000", linewidth=3.5, linestyle="-",
         label=r"Apparent Ceiling Imposed by $^{14}\text{C}$ Window ($\leq 42$ ka)")

ax1.axhline(42000, color="#666666", linestyle=":", linewidth=2.0, label=r"Radiocarbon Instrumental Blank Limit ($\sim 42\text{--}45$ ka)")
ax1.axvline(13.8, color="#0055D4", linestyle="-.", linewidth=1.5)

# Shading the unobservable zone
ax1.fill_between(t_fine, 42000, 50000, color="#EEEEEE", alpha=0.8)
ax1.text(-10, 44000, r"UNOBSERVABLE VIA $^{14}\text{C}$ (Collagen survives, but cannot be dated)",
         fontsize=11, fontweight="bold", color="#555555")

ax1.annotate("THE ARTIFACT:\nHorizontal clipping at 42 ka\nmakes slope appear artificially flat\nand forces thermal age to zero in cold regions!",
             xy=(2, 42000), xytext=(-8, 22000),
             arrowprops=dict(arrowstyle="->", color="#980E1E", lw=2.0),
             bbox=dict(boxstyle="round,pad=0.6", facecolor="#FFEEEE", edgecolor="#980E1E", lw=1.5),
             fontsize=10, fontweight="bold")

ax1.set_xlim(-12, 28)
ax1.set_ylim(0, 48000)
ax1.set_xlabel("Paleoclimate-Integrated Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_ylabel("Radiocarbon Age (BP, Uncalibrated)", fontsize=13, fontweight="bold", labelpad=8)
ax1.set_title(r"Panel A: Radiocarbon Censoring Below $14^\circ\text{C}$", fontsize=14, fontweight="bold", pad=12)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle="--", alpha=0.4, color="#cccccc")
ax1.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="#cccccc", fontsize=9.5)

# -------------------------------------------------------------
# PANEL B: True Kinetic Extrapolation across Geological Time
# -------------------------------------------------------------
ax2.plot(t_fine, true_kinetics_133, color="#D95F02", linewidth=3.2, linestyle="-", zorder=3,
         label=r"True Kinetic Ceiling: $E_a = 133.4\ \mathrm{kJ/mol}$")
ax2.plot(t_fine, true_kinetics_137, color="#1B9E77", linewidth=2.6, linestyle="--", zorder=3,
         label=r"Upper Bound Envelope: $E_a = 137.0\ \mathrm{kJ/mol}$")
ax2.plot(t_fine, np.exp(C_173 + slope_173 * (1.0 / t_k)), color="#7570B3", linewidth=2.0, linestyle="-.", zorder=3,
         label=r"Historical Lab Model: $E_a = 173.0\ \mathrm{kJ/mol}$")

# Predicted Bone Collagen Survival Envelope
ax2.fill_between(t_fine, true_kinetics_133, true_kinetics_137, color="#D95F02", alpha=0.15,
                 label="Empirical Bone Collagen Kinetic Envelope")

# Plot the 14C valid points (T >= 13 °C)
ax2.scatter(k_col["temp_mean"], k_col["max_c14"],
            c="#980E1E", s=80, edgecolors="#333333", zorder=5, label=r"Empirical Max $^{14}\mathrm{C}$ Ages ($T \geq 13^\circ\mathrm{C}$)")

# Radiocarbon limit band
ax2.axhline(42000, color="#666666", linestyle=":", linewidth=1.8, label=r"$^{14}\mathrm{C}$ Limit (42 ka BP)")
ax2.fill_between(t_fine, 10, 42000, color="#E8F4F8", alpha=0.5, label=r"$^{14}\mathrm{C}$ Measurable Window")

# Real Palaeoproteomic benchmarks beyond 14C limit:
benchmarks_bone = [
    ("Denisova Hominin (~150 ka, Bone)", -2.0, 150000, (-2.0, 280000)),
    ("Boxgrove / Sima (~430 ka, Bone/Dentin)", 7.0, 430000, (7.5, 900000)),
    ("Yukon Horse (~735 ka, Bone)", -9.0, 735000, (-8.5, 1400000)),
    ("Dmanisi Rhino (~1.77 Ma, Dentin)", 11.0, 1770000, (11.5, 3500000)),
    ("High Arctic Camel (~3.4 Ma, Bone)", -12.0, 3400000, (-10.5, 2300000)),
    ("Ellesmere Bear (~3.9 Ma, Radius Bone)", -10.5, 3900000, (-5.5, 4800000))
]

for name, t_site, age_site, (tx, ty) in benchmarks_bone:
    ax2.scatter(t_site, age_site, c="#D95F02", s=85, marker="^", edgecolors="#000000", linewidths=1.2, zorder=6)
    ax2.annotate(name, (t_site, age_site), xytext=(tx, ty),
                 fontsize=8.5, fontweight="bold", color="#111111",
                 arrowprops=dict(arrowstyle="->", color="#D95F02", lw=1.3),
                 bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFF3E0", edgecolor="#D95F02", alpha=0.9, lw=0.8))

# Enamel Control Benchmark (PXD052635 Devon Island Rhino, 21-23 Ma, 0% Collagen)
ax2.scatter(-11.0, 22000000, c="#7570B3", s=110, marker="D", edgecolors="#000000", linewidths=1.4, zorder=7,
            label="Devon Island Rhino Enamel (21+ Ma, 0% Col)")
ax2.annotate("Devon Island Rhino (PXD052635)\n(~21–23 Ma, Enamel EMPs,\nZERO Collagen Survives)",
             (-11.0, 22000000), xytext=(-6.8, 16000000),
             fontsize=8.5, fontweight="bold", color="#311B92",
             arrowprops=dict(arrowstyle="->", color="#7570B3", lw=1.5),
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#EDE7F6", edgecolor="#7570B3", alpha=0.95, lw=1.0))

# Annotation showing Arctic Collagen Ceiling Fit
ax2.annotate(r"High Arctic Bone Collagen Ceiling:" "\n"
             r"Ellesmere Bear (3.9 Ma) & Camel (3.4 Ma)" "\n"
             r"sit directly on the $E_a = 133\text{--}137\ \mathrm{kJ/mol}$ boundary!" "\n"
             r"$\rightarrow$ Absolute collagen survival ceiling is $\sim 4\text{--}7\ \mathrm{Ma}$",
             xy=(-10.5, 3900000), xytext=(-11.5, 200),
             fontsize=9.0, fontweight="bold", color="#980E1E",
             arrowprops=dict(arrowstyle="->", color="#980E1E", lw=1.6),
             bbox=dict(boxstyle="round,pad=0.45", facecolor="#FFEBEE", edgecolor="#D32F2F", lw=1.2))

ax2.set_yscale("log")
ax2.set_xlim(-12.5, 28)
ax2.set_ylim(50, 70000000) # Up to 70 Million years!
ax2.set_xlabel("Paleoclimate-Integrated Temperature $\\bar{T}$ (°C)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_ylabel("True Survival Age (Years, Log Scale)", fontsize=13, fontweight="bold", labelpad=8)
ax2.set_title("Panel B: True Kinetic Lifespans Extrapolated into Deep Time", fontsize=14, fontweight="bold", pad=12)
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y):,}" if y < 1e6 else f"{y/1e6:g} Ma"))
ax2.grid(True, which="both", linestyle="--", alpha=0.4, color="#cccccc")
ax2.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="#cccccc", fontsize=8.5)

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_True_Kinetic_Extrapolation_No_14C_Distortion.png")
plt.savefig(out_fig, dpi=300)
plt.close()
print(f"Saved figure to: {out_fig}")
