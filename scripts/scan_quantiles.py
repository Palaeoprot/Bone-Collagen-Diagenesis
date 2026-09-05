import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

p = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(p)
col = df[df["material_category"] == "COLLAGEN"].copy()
col["log_thermal_age"] = np.log(col["thermal_age_173"])
r_gas = 8.314462618e-3

qs = [0.50, 0.70, 0.80, 0.85, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99]
for q in qs:
    mod = smf.quantreg("log_thermal_age ~ integrated_temp_c", col).fit(q=q)
    b = mod.params["integrated_temp_c"]
    b_se = mod.bse["integrated_temp_c"]
    a0 = np.exp(mod.params["Intercept"])
    ea = b * r_gas * (283.15**2)
    ea_se = b_se * r_gas * (283.15**2)
    print(f"Quantile {int(q*100):2d}%: Slope b = {b:.4f} +/- {b_se:.4f} | A0 = {a0:7.1f} | Ea = {ea:5.1f} +/- {ea_se:.1f} kJ/mol")
