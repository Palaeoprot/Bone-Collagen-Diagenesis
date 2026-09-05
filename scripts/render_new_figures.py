import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Set high-quality publication styling
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

fig_dir = r"C:\Users\matth\Documents\GitHub\Palaeoprot-Publications\projects\Ea Collagen Hydrolysis\figures"
os.makedirs(fig_dir, exist_ok=True)

# ==============================================================================
# FIGURE 1: BIOPOLYMER DEGRADATION AXIS (COLLAGEN VS DNA ARRHENIUS KINETICS)
# ==============================================================================
fig, ax = plt.subplots(figsize=(10, 7), dpi=300)

r_gas = 8.314462618e-3 # kJ / (mol * K)
t_c_range = np.linspace(-15, 45, 200)
t_k_range = t_c_range + 273.15
inv_t_range = 1000.0 / t_k_range

# 1. Smith (2002) Powdered Bone Lab Heating (Ea = 173.2 kJ/mol, ln A = 44.49)
ea_smith = 173.2
a_smith = np.exp(44.49) # per year
ln_k_smith = np.log(a_smith) - (ea_smith / (r_gas * t_k_range))
half_life_smith = np.log(2) / np.exp(ln_k_smith)

# 2. Ortner et al. (1972) Intact Bone Lab Benchmark (Ea = 132.1 kJ/mol)
# ln A ~ 31.0 per year
ea_ortner = 132.1
a_ortner = np.exp(31.2)
ln_k_ortner = np.log(a_ortner) - (ea_ortner / (r_gas * t_k_range))
half_life_ortner = np.log(2) / np.exp(ln_k_ortner)

# 3. Allentoft et al. (2012) Ancient DNA in Bone (Ea = 127.3 kJ/mol, half-life = 521 y @ 13.1°C)
ea_dna = 127.3
k_dna_13 = np.log(2) / 521.0
a_dna = k_dna_13 / np.exp(-ea_dna / (r_gas * (13.1 + 273.15)))
ln_k_dna = np.log(a_dna) - (ea_dna / (r_gas * t_k_range))
half_life_dna = np.log(2) / np.exp(ln_k_dna)

# 4. Global Radiocarbon 1% Empirical Field Barrier (Ea = 82.2 kJ/mol)
ea_field = 82.2
# 42,000 yr @ 14°C -> k_limit = 4.60517 / 42000
k_field_14 = 4.60517 / 42000.0
a_field = k_field_14 / np.exp(-ea_field / (r_gas * (14.0 + 273.15)))
ln_k_field = np.log(a_field) - (ea_field / (r_gas * t_k_range))
half_life_field = np.log(2) / np.exp(ln_k_field)

# Plot lines
ax.plot(t_c_range, half_life_smith, label='Powdered Bone Denaturation (Smith 2002; $E_a = 173.2$ kJ/mol)', color='#8e44ad', linewidth=2.5, linestyle='--')
ax.plot(t_c_range, half_life_ortner, label='Intact Bone Scaffolding (Ortner 1972; $E_a = 132.1$ kJ/mol)', color='#2980b9', linewidth=2.5)
ax.plot(t_c_range, half_life_dna, label='Ancient DNA in Bone (Allentoft 2012; $E_a = 127.3$ kJ/mol)', color='#e67e22', linewidth=2.5, linestyle='-.')
ax.plot(t_c_range, half_life_field, label='Global Empirical 1% Bone Collagen Limit ($E_a = 82.2$ kJ/mol)', color='#27ae60', linewidth=2.5)

# Add deep time empirical survivors
benchmarks = [
    ("Ellesmere Bear (3.9 Ma, -10.5°C)", -10.5, 3.9e6, '#2c3e50'),
    ("High Arctic Camel (3.4 Ma, -10.5°C)", -10.5, 3.4e6, '#2c3e50'),
    ("Dmanisi Dentin (1.77 Ma, 11°C)", 11.0, 1.77e6, '#c0392b'),
    ("Yukon Horse (735 ka, -9.0°C)", -9.0, 735e3, '#2c3e50'),
    ("Sima de los Huesos (430 ka, 7°C)", 7.0, 430e3, '#2980b9'),
    ("Harbin Cranium (148 ka, 3.5°C)", 3.5, 148e3, '#27ae60'),
    ("Denisova Bone (45 ka, -3.3°C)", -3.3, 45e3, '#27ae60'),
]

for name, tc, age, col in benchmarks:
    ax.scatter(tc, age, color=col, s=70, zorder=5, edgecolors='black')
    ax.annotate(name, (tc, age), textcoords="offset points", xytext=(8, -3), fontsize=8.5, fontweight='bold', color=col)

# Formatting
ax.set_yscale('log')
ax.set_xlim(-15, 40)
ax.set_ylim(1e1, 1e9)
ax.set_xlabel('Effective Degradation Temperature (°C)', fontsize=12, fontweight='bold')
ax.set_ylabel('Biopolymer Half-Life / Survival Horizon (Years, Log Scale)', fontsize=12, fontweight='bold')
ax.set_title('Biopolymer Degradation Kinetics: Bone Collagen vs Ancient DNA in Deep Time', fontsize=13, fontweight='bold', pad=12)
ax.grid(True, which='both', alpha=0.3)
ax.axhline(50000, color='gray', linestyle=':', label='Radiocarbon Instrumental Limit (~50 ka)')
ax.axvspan(37, 40, color='#fadbd8', alpha=0.3, label='Human Body Temperature (37°C)')

ax.legend(loc='upper right', framealpha=0.95, fontsize=9)
plt.tight_layout()

f1_png = os.path.join(fig_dir, "Figure_Biopolymer_Kinetics_Collagen_vs_DNA_Arrhenius.png")
f1_svg = os.path.join(fig_dir, "Figure_Biopolymer_Kinetics_Collagen_vs_DNA_Arrhenius.svg")
fig.savefig(f1_png, dpi=300)
fig.savefig(f1_svg)
plt.close(fig)
print("Figure 1 saved:", f1_png)

# ==============================================================================
# FIGURE 2: SEASONAL ACCELERATION & CAVE BUFFERING
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# Panel A: Seasonal Wave & Arrhenius Amplification
theta = np.linspace(0, 2*np.pi, 365)
mat = 0.0 # MAT = 0°C (Siberian / Continental Plain)
amp = 18.0 # Annual range = 36°C (-18°C to +18°C)
t_wave = mat + amp * np.sin(theta)
t_wave_k = t_wave + 273.15

# Arrhenius rate profile
rate_wave = np.exp(-82.2 / (r_gas * t_wave_k))
mean_rate = np.mean(rate_wave)
t_eff_k = -82.2 / (r_gas * np.log(mean_rate))
t_eff_c = t_eff_k - 273.15

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
day_x = np.linspace(0, 12, 365)

ax1.plot(day_x, t_wave, color='#2980b9', linewidth=2.5, label=r'Surface Seasonal Temperature ($T_{\mathrm{range}} = 36$ °C)')
ax1.axhline(mat, color='#2c3e50', linestyle='--', linewidth=2, label=f'Cave / Buffered MAT ({mat:.1f} °C)')
ax1.axhline(t_eff_c, color='#e74c3c', linestyle='-', linewidth=2.5, label=fr'Effective Surface Kinetic $T_{{\mathrm{{eff}}}}$ (+{t_eff_c:.1f} °C)')

ax1.fill_between(day_x, mat, t_wave, where=(t_wave >= mat), color='#e74c3c', alpha=0.15, label='Summer Kinetic Surge (Dominates Decay)')
ax1.fill_between(day_x, mat, t_wave, where=(t_wave < mat), color='#3498db', alpha=0.15, label='Winter Retardation')

ax1.set_xticks(np.linspace(0.5, 11.5, 12))
ax1.set_xticklabels(months)
ax1.set_ylabel('Temperature (°C)', fontsize=11, fontweight='bold')
ax1.set_title('A: Arrhenius Exponential Summer Acceleration', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='lower center', fontsize=8.5, framealpha=0.95)

# Panel B: Benchmark Sites Cave vs Surface Comparison
sites = ['Denisova\n(Altai)', 'Sunghir\n(Russian Pl.)', 'Vindija\n(Croatia)', 'Guattari\n(Italy)', 'Hayonim\n(Levant)']
cave_temps = [-3.3, -0.7, 7.3, 13.3, 16.2]
surf_temps = [11.8, 10.3, 12.6, 17.1, 20.5]
accel_factors = [6.99, 4.07, 1.91, 1.57, 1.65]

x_pos = np.arange(len(sites))
width = 0.35

rects1 = ax2.bar(x_pos - width/2, cave_temps, width, label=r'Cave Buffered MAT ($T_{\mathrm{cave}}$)', color='#2980b9', alpha=0.9)
rects2 = ax2.bar(x_pos + width/2, surf_temps, width, label=r'Surface Seasonal $T_{\mathrm{eff}}$', color='#e74c3c', alpha=0.9)

for i, factor in enumerate(accel_factors):
    y_top = surf_temps[i] + 1.0
    ax2.annotate(f"{factor:.1f}×\nacceleration", (x_pos[i], y_top), ha='center', fontsize=8.5, fontweight='bold', color='#8e44ad')

ax2.set_xticks(x_pos)
ax2.set_xticklabels(sites, fontsize=9.5)
ax2.set_ylabel('Effective Degradation Temperature (°C)', fontsize=11, fontweight='bold')
ax2.set_ylim(-6, 25)
ax2.set_title('B: Paleoclimate Thermal Brackets for Benchmark Sites', fontsize=12, fontweight='bold')
ax2.grid(True, axis='y', alpha=0.3)
ax2.legend(loc='upper left', fontsize=9, framealpha=0.95)

plt.tight_layout()
f2_png = os.path.join(fig_dir, "Figure_Seasonal_Acceleration_and_Cave_Buffering.png")
f2_svg = os.path.join(fig_dir, "Figure_Seasonal_Acceleration_and_Cave_Buffering.svg")
fig.savefig(f2_png, dpi=300)
fig.savefig(f2_svg)
plt.close(fig)
print("Figure 2 saved:", f2_png)

# ==============================================================================
# FIGURE 3: TAPHONOMIC CAVE FILTERING CHRONOLOGY
# ==============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

windows = ['< 5 ka\n(Late Holocene)', '5–15 ka\n(Early Holo/Lategl.)', '15–30 ka\n(LGM & Upper Pal.)', '> 30 ka\n(Mid/Early Upper Pal.)']
total_counts = [35802, 19019, 6456, 5980]
cave_counts = [1110, 2114, 1144, 1830]
open_counts = [34692, 16905, 5312, 4150]
cave_pcts = [3.1, 11.1, 17.7, 30.6]

# Panel A: Stacked determination counts
x_w = np.arange(len(windows))
w_bar = 0.55

ax1.bar(x_w, open_counts, w_bar, label='Open-Air Determinations', color='#7f8c8d', alpha=0.85)
ax1.bar(x_w, cave_counts, w_bar, bottom=open_counts, label='Cave / Karst Determinations', color='#2980b9', alpha=0.95)

for i in range(len(windows)):
    tot = total_counts[i]
    ax1.annotate(f"N = {tot:,}", (x_w[i], tot + 800), ha='center', fontsize=9, fontweight='bold')

ax1.set_xticks(x_w)
ax1.set_xticklabels(windows, fontsize=9.5)
ax1.set_ylabel('Number of Radiocarbon Determinations', fontsize=11, fontweight='bold')
ax1.set_title('A: Global Radiocarbon Cohort Composition ($N = 67,257$)', fontsize=12, fontweight='bold')
ax1.grid(True, axis='y', alpha=0.3)
ax1.legend(loc='upper right', fontsize=9.5)

# Panel B: Cave Proportion Rise
ax2.plot(x_w, cave_pcts, marker='o', markersize=9, color='#e74c3c', linewidth=3, label='Cave Proportion (%)')
for i, pct in enumerate(cave_pcts):
    ax2.annotate(f"{pct:.1f}%", (x_w[i], pct + 1.2), ha='center', fontsize=10.5, fontweight='bold', color='#c0392b')

ax2.set_xticks(x_w)
ax2.set_xticklabels(windows, fontsize=9.5)
ax2.set_ylabel('Proportion of Dates from Caves (%)', fontsize=11, fontweight='bold')
ax2.set_ylim(0, 38)
ax2.set_title('B: Progressive Deep-Time Taphonomic Cave-Filtering', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)

# Highlight tenfold increase
ax2.annotate('10-fold increase in cave reliance\ndue to open-air collagen destruction',
             xy=(3, 30.6), xytext=(1.5, 26),
             arrowprops=dict(facecolor='black', shrink=0.08, width=1.5, headwidth=8),
             fontsize=9.5, fontweight='bold', bbox=dict(boxstyle="round,pad=0.4", fc="#fef9e7", ec="#f39c12", lw=1.5))

plt.tight_layout()
f3_png = os.path.join(fig_dir, "Figure_Taphonomic_Cave_Filtering_Chronology.png")
f3_svg = os.path.join(fig_dir, "Figure_Taphonomic_Cave_Filtering_Chronology.svg")
fig.savefig(f3_png, dpi=300)
fig.savefig(f3_svg)
plt.close(fig)
print("Figure 3 saved:", f3_png)
