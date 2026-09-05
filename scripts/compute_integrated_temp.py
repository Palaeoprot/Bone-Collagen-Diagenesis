import os, time
import pandas as pd
import numpy as np
from paleoclimate_engine import PaleoclimateEngine

data_dir = r"D:\26 Modelling Collagen Hydrolysis\data"
cohort_path = os.path.join(data_dir, "collagen_vs_control_thermal_cohort.parquet")
print(f"Loading cohort from {cohort_path}...")
df = pd.read_parquet(cohort_path)

engine = PaleoclimateEngine()

coords = df[["latitude", "longitude"]].drop_duplicates().values
print(f"Caching coordinates across {len(coords)} locations...")
coord_cache = {}
for lat, lon in coords:
    t_bp, tc, tk = engine.get_temperature_series(lat, lon)
    coord_cache[(lat, lon)] = (t_bp, tc)

# Calculate mean integrated temperature over sample lifetime for each record:
# T_mean = (1 / t_cal) * \int_0^{t_cal} T(t) dt
mean_temps = []
t0 = time.time()
print("Calculating mean integrated temperature for each sample...")
for idx, r in df.iterrows():
    lat = r["latitude"]
    lon = r["longitude"]
    cal_age = r["c14_age"]
    
    times_bp, tc = coord_cache[(lat, lon)]
    mask = times_bp <= cal_age
    st = times_bp[mask].tolist()
    sc = tc[mask].tolist()
    
    if len(st) == 0 or st[-1] < cal_age:
        t_interp = np.interp(cal_age, times_bp, tc)
        st.append(cal_age)
        sc.append(t_interp)
        
    st = np.array(st)
    sc = np.nan_to_num(np.array(sc), nan=10.0)
    
    if len(st) >= 2:
        integ_t = np.trapezoid(sc, st) if hasattr(np, "trapezoid") else np.trapz(sc, st)
        mean_t = integ_t / cal_age
    else:
        mean_t = sc[0]
        
    mean_temps.append(float(mean_t))

df["integrated_temp_c"] = mean_temps
print(f"Completed in {(time.time()-t0):.1f}s.")

df.to_parquet(cohort_path, index=False)
print(f"Saved updated cohort with integrated_temp_c to: {cohort_path}")
print(df[["material_category", "integrated_temp_c", "thermal_age_173"]].groupby("material_category").describe())
