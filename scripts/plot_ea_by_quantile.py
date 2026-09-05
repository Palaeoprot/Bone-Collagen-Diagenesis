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

# Model fits across 90%, 95%, 98%, 99%
models = {}
for q in [0.90, 0.95, 0.98, 0.99]:
    m = smf.quantreg("log_thermal_age ~ integrated_temp_c", col).fit(q=q)
    b = m.params["integrated_temp_c"]
    b_se = m.bse["integrated_temp_c"]
    a0 = np.exp(m.params["Intercept"])
    ea = b * r_gas * (t_ref_k**2)
    models[q] = {"model": m, "b": b, "b_se": b_se, "a0": a0, "ea": ea}

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0

fig, axes = plt.subplots(1, 2, figsize=(16, 7.5), dpi=300)

t_min, t_max = -15, 25
th_max_linear = 220000
temp_grid = np.linspace(t_min, t_max, 200)

# Colors for the curves
colors_dict = {
    0.90: "#1b5e20",  # Dark Green (90%)
    0.95: "#C41E23",  # Vibrant Red (95% primary)
    0.98: "#ff9800",  # Orange (98%)
    0.99: "#9c27b0"   # Purple (99%)
}

# --- PANEL A: Purified Bone Collagen with 90%, 95%, 98%, 99% Curves ---
ax1 = axes[0]
sc1 = ax1.scatter(col["integrated_temp_c"], col["thermal_age_173"],
                  c=col["c14_age"], cmap="viridis", alpha=0.30, s=12, rasterized=True)
cbar1 = plt.colorbar(sc1, ax=ax1, pad=0.02)
cbar1.set_label("Calendar Age ($^{14}\\mathrm{C}$ BP)", fontsize=10, fontweight="bold")
cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

for q in [0.90, 0.95, 0.98, 0.99]:
    d = models[q]
    curve_y = d["a0"] * np.exp(d["b"] * temp_grid)
    ls = "-" if q == 0.95 else ("-." if q == 0.90 else "--")
    lw = 3.0 if q == 0.95 else 2.0
    label = f"{int(q*100)}% Boundary: $E_a \\approx {d['ea']:.1f}\\ \\mathrm{{kJ/mol}}$ ($b = {d['b']:.4f}$)"
    ax1.plot(temp_grid, curve_y, color=colors_dict[q], lw=lw, linestyle=ls, label=label)

ax1.set_xlim(t_min, 30)
ax1.set_ylim(0, th_max_linear)
ax1.set_title(r"$\mathbf{A}$  Purified Bone Collagen: $E_a$ Variation Across 90%–99% Envelopes", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax1.set_xlabel("Integrated Paleotemperature (°C over sample lifetime)", fontsize=11, fontweight="bold")
ax1.set_ylabel("Arrhenius Thermal Age (Years @ 10 °C — Linear Scale)", fontsize=11, fontweight="bold")
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

# --- PANEL B: Activation Energy vs Quantile Progression ---
ax2 = axes[1]
all_qs = [0.50, 0.70, 0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99]
ea_list = []
ea_err = []
b_list = []

for q in all_qs:
    mq = smf.quantreg("log_thermal_age ~ integrated_temp_c", col).fit(q=q)
    bq = mq.params["integrated_temp_c"]
    bq_se = mq.bse["integrated_temp_c"]
    ea_val = bq * r_gas * (t_ref_k**2)
    ea_se_val = bq_se * r_gas * (t_ref_k**2)
    ea_list.append(ea_val)
    ea_err.append(ea_se_val)
    b_list.append(bq)

q_pcts = [q * 100 for q in all_qs]

# Plot Ea vs Quantile
ax2.errorbar(q_pcts, ea_list, yerr=ea_err, fmt='o-', color="#980E1E", ecolor="#C41E23", 
             elinewidth=2.0, capsize=4, markersize=7, lw=2.2, label=r"Empirical Field $E_a(q)$")

# Reference lines
ax2.axhline(173.0, color="#111111", linestyle=":", lw=1.8, label="Collins Lab Kinetic Benchmark (173 kJ/mol)")
ax2.axhline(111.6, color="#009688", linestyle="--", lw=1.8, label="Triple Helix QM/MM Barrier (Buhr & Gräter 2026, 111.6 kJ/mol)")
ax2.axhline(100.0, color="#607d8b", linestyle="-.", lw=1.5, label="Uncatalyzed Solution Dipeptide (100 kJ/mol)")

# Highlight 90%, 95%, 99% points
for target_q, colr in [(90, "#1b5e20"), (95, "#C41E23"), (99, "#9c27b0")]:
    idx = q_pcts.index(target_q)
    val = ea_list[idx]
    ax2.plot(target_q, val, 'o', color=colr, markersize=10, zorder=5)
    ax2.annotate(f"{target_q}%: {val:.1f} kJ/mol", (target_q, val),
                 textcoords="offset points", xytext=(-25, 12 if target_q!=99 else -20),
                 fontsize=9.5, fontweight="bold", color=colr,
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=colr, alpha=0.9))

ax2.set_xlim(48, 100)
ax2.set_ylim(60, 190)
ax2.set_title(r"$\mathbf{B}$  Effective Activation Energy ($E_a$) Progression Across Quantiles", 
             fontsize=13, fontweight="bold", pad=12, loc="left")
ax2.set_xlabel("Data Containment Quantile (%)", fontsize=11, fontweight="bold")
ax2.set_ylabel(r"Effective Activation Energy $E_a\ (\mathrm{kJ\cdot mol^{-1}}\ \mathrm{around}\ 10^\circ\mathrm{C})$", fontsize=11, fontweight="bold")
ax2.grid(True, linestyle=":", alpha=0.5)
ax2.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.0, loc="lower left")

plt.tight_layout()
out_fig = os.path.join(output_dir, "Figure_Ea_Variation_90_95_99_Quantiles.png")
plt.savefig(out_fig, dpi=300)
print(f"Saved figure to: {out_fig}")

# Save numerical table
df_summary = pd.DataFrame({
    "Quantile_Pct": q_pcts,
    "Slope_b": b_list,
    "Ea_kJ_mol": ea_list,
    "Ea_SE": ea_err
})
df_summary.to_csv(os.path.join(output_dir, "ea_by_quantile_summary.csv"), index=False)
print("Saved summary table to ea_by_quantile_summary.csv")
