import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

p = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(p)
col = df[df["material_category"] == "COLLAGEN"].copy()
ctrl = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

print("=== ZERO-INTERCEPT (0, 0) QUANTILE REGRESSION: THERMAL AGE VS CALENDAR AGE ===")

# For Ea = 173 kJ/mol
for q in [0.50, 0.75, 0.90, 0.95, 0.98, 0.99]:
    m_col = smf.quantreg("thermal_age_173 ~ c14_age - 1", col).fit(q=q)
    m_ctrl = smf.quantreg("thermal_age_173 ~ c14_age - 1", ctrl).fit(q=q)
    slope_col = m_col.params["c14_age"]
    slope_ctrl = m_ctrl.params["c14_age"]
    print(f"Quantile {q*100:4.1f}% | Collagen Slope: {slope_col:7.4f} | Control Slope: {slope_ctrl:7.4f}")

# At calendar age = 40,000 BP (the practical radiocarbon detection horizon for collagen):
m95 = smf.quantreg("thermal_age_173 ~ c14_age - 1", col).fit(q=0.95)
s95 = m95.params["c14_age"]
th_limit_40k = s95 * 40000.0
print(f"\nAt 95th % envelope through (0,0):")
print(f"Slope: {s95:.4f}")
print(f"Thermal Age at 40,000 BP = {th_limit_40k:,.1f} years (@ 10°C)")

# In first-order kinetics: ln(C/C0) = -k * t_eff
# If C/C0 = 0.01 (1% residual collagen), then -ln(0.01) = 4.60517
# k_eff(10°C) = 4.60517 / t_thermal_max
k_10c = 4.60517 / th_limit_40k
t_half_10c = np.log(2) / k_10c
print(f"\nFirst-order kinetic parameters implied by 1% threshold at this boundary:")
print(f"k(10°C) = {k_10c:.4e} yr^-1 = {k_10c / (365.25*86400):.4e} s^-1")
print(f"Half-life at 10°C = {t_half_10c:,.1f} years")

# What pre-exponential factor A does this imply for Ea = 173 kJ/mol?
# k = A * exp(-Ea / (R * T)) -> A = k * exp(Ea / (R * T))
r_gas = 8.314462618e-3 # kJ / (mol * K)
t_ref_k = 283.15 # 10 °C
a_factor = k_10c * np.exp(173.0 / (r_gas * t_ref_k))
a_factor_s = (k_10c / (365.25*86400)) * np.exp(173.0 / (r_gas * t_ref_k))
print(f"Implied Pre-exponential factor A:")
print(f"A = {a_factor:.4e} yr^-1 = {a_factor_s:.4e} s^-1")
