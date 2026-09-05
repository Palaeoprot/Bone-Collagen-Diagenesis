import os, time
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

output_dir = r"D:\26 Modelling Collagen Hydrolysis\outputs"
data_path = r"D:\26 Modelling Collagen Hydrolysis\data\collagen_vs_control_thermal_cohort.parquet"
df = pd.read_parquet(data_path)

collagen = df[df["material_category"] == "COLLAGEN"].copy()
controls = df[df["material_category"] == "NON_COLLAGEN_ORGANIC_CONTROL"].copy()

from paleoclimate_engine import PaleoclimateEngine
engine = PaleoclimateEngine()

# Pre-cache unique coordinates
coords = collagen[["latitude", "longitude"]].drop_duplicates().values
coord_cache = {}
for lat, lon in coords:
    t_bp, tc, tk = engine.get_temperature_series(lat, lon)
    coord_cache[(lat, lon)] = (t_bp, tk)

def compute_thermal_ages_fast(df_sub, ea_val, t_ref_c=10.0):
    t_ref_k = t_ref_c + 273.15
    r_gas = 8.314462618e-3
    ages = []
    for idx, r in df_sub.iterrows():
        lat = r["latitude"]
        lon = r["longitude"]
        cal_age = r["c14_age"]
        t_bp, tk = coord_cache[(lat, lon)]
        mask = t_bp <= cal_age
        st = t_bp[mask].tolist()
        sk = tk[mask].tolist()
        if len(st) == 0 or st[-1] < cal_age:
            st.append(cal_age)
            sk.append(np.interp(cal_age, t_bp, tk))
        st = np.array(st)
        sk = np.nan_to_num(np.array(sk), nan=t_ref_k)
        arrh = np.exp(-(ea_val / r_gas) * (1.0 / sk - 1.0 / t_ref_k))
        if len(st) >= 2:
            th = np.trapezoid(arrh, st) if hasattr(np, "trapezoid") else np.trapz(arrh, st)
        else:
            th = arrh[0] * cal_age
        ages.append(float(th))
    return np.array(ages)

# Grid search Ea from 100 to 220 kJ/mol in steps of 10 kJ/mol
ea_grid = [100.0, 120.0, 140.0, 150.0, 160.0, 170.0, 173.0, 180.0, 190.0, 200.0]
sweep_results = []

print("Running empirical Ea parameter sweep on collagen boundary dispersion...")
for ea in ea_grid:
    th_ages = compute_thermal_ages_fast(collagen, ea)
    collagen[f"th_{int(ea)}"] = th_ages
    
    # 95th percentile upper boundary
    q95 = np.percentile(th_ages, 95)
    q50 = np.percentile(th_ages, 50)
    iqr = np.percentile(th_ages, 75) - np.percentile(th_ages, 25)
    
    # Quantile regression slope
    mod = smf.quantreg(f"th_{int(ea)} ~ c14_age", collagen).fit(q=0.95)
    slope = mod.params["c14_age"]
    
    # Normalized boundary sharpness (Coefficient of variation of upper 5% boundary)
    upper_tail = th_ages[th_ages >= q95]
    cv_tail = np.std(upper_tail) / np.mean(upper_tail)
    
    sweep_results.append({
        "Ea_kJ_mol": ea,
        "Slope_95": slope,
        "Q95_Thermal_Age": q95,
        "CV_Upper_Tail": cv_tail
    })
    print(f"Ea = {ea:5.1f} kJ/mol | 95th Slope: {slope:6.3f} | Q95: {q95:8.1f} y | CV Tail: {cv_tail:6.3f}")

sweep_df = pd.DataFrame(sweep_results)
sweep_path = os.path.join(output_dir, "ea_grid_search_boundary_sharpness.csv")
sweep_df.to_csv(sweep_path, index=False)
print(f"\nSaved grid search results to: {sweep_path}")
