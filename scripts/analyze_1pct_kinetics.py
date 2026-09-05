import pandas as pd
import numpy as np

p = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(p)
col = df[df["material_category"] == "COLLAGEN"].copy()
ctrl = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

print("=== INVERSE FORMULATION: CALENDAR AGE VS THERMAL AGE THROUGH (0, 0) ===")
# Calendar Age (x) as a function of Thermal Age (t_eff)
# Or fit Cal Age ~ Thermal Age - 1
import statsmodels.formula.api as smf

# Fit how calendar age relates to thermal age
m_inv = smf.quantreg("c14_age ~ thermal_age_173 - 1", col).fit(q=0.95)
print(f"95th % Slope (Calendar Age / Thermal Age): {m_inv.params['thermal_age_173']:.4f}")

# Distribution of thermal ages across collagen vs controls:
print("\nThermal Age Percentiles for Collagen:")
for q in [50, 75, 90, 95, 99]:
    print(f"  {q}th percentile: {np.percentile(col['thermal_age_173'], q):,.1f} years (@ 10°C)")

print("\nThermal Age Percentiles for Controls:")
for q in [50, 75, 90, 95, 99]:
    print(f"  {q}th percentile: {np.percentile(ctrl['thermal_age_173'], q):,.1f} years (@ 10°C)")

# For various thermal age cutoffs, what is the rate constant k(10°C) for 1% residual collagen?
# ln(1% / 100%) = -4.60517 -> k = 4.60517 / t_eff
r_gas = 8.314462618e-3
t_ref_k = 283.15
arrh_term = np.exp(173.0 / (r_gas * t_ref_k))

print("\n=== KINETIC PARAMETERS FOR 1% RESIDUAL COLLAGEN AT DIFFERENT BOUNDARY PERCENTILES ===")
print("Quantile | Thermal Horizon (y @ 10°C) | k(10°C) [yr^-1] | k(10°C) [s^-1] | Half-life [yr] | Pre-exponential A [s^-1]")
for q in [90, 95, 98, 99]:
    t_h = np.percentile(col["thermal_age_173"], q)
    k_yr = 4.60517 / t_h
    k_s = k_yr / (365.25 * 86400)
    t_half = np.log(2) / k_yr
    A_s = k_s * arrh_term
    print(f" {q:2d}%    | {t_h:15,.1f}            | {k_yr:12.4e}    | {k_s:12.4e}   | {t_half:10,.1f}     | {A_s:12.4e}")
