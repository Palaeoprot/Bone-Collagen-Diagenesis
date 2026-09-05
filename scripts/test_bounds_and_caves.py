from paleoclimate_engine import PaleoclimateEngine
import pandas as pd
import numpy as np

engine = PaleoclimateEngine()
benchmarks = [
    {'name': 'Denisova Cave (Altai)', 'lat': 51.3975, 'lon': 84.6761, 'age': 45000, 'elev': 700},
    {'name': 'Vindija Cave (Croatia)', 'lat': 46.2992, 'lon': 16.0711, 'age': 38000, 'elev': 275},
    {'name': 'Sunghir (Russia, Open-Air)', 'lat': 56.1756, 'lon': 40.5058, 'age': 34000, 'elev': 170},
    {'name': 'Grotta Guattari (Italy)', 'lat': 41.2286, 'lon': 13.0903, 'age': 50000, 'elev': 5},
    {'name': 'Hayonim Cave (Levant)', 'lat': 32.9233, 'lon': 35.2167, 'age': 28000, 'elev': 250}
]

print("\n### Benchmark Seasonal Bounds (Ea = 82.2 kJ/mol, T_ref = 10 °C)\n")
print("| Site | Context | Coordinates | Calendar Age (BP) | Lower Bound Thermal Age (Cave/MAT) | Upper Bound Thermal Age (Surface Seasonal) | Seasonal Acceleration | Cave T_eff | Surface T_eff | ΔT |")
print("| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
for b in benchmarks:
    res = engine.calculate_thermal_age_bounds(b['lat'], b['lon'], b['age'], ea_kj=82.2, t_ref_c=10.0, elevation_m=b['elev'])
    ctx = 'Cave' if ('Cave' in b['name'] or 'Grotta' in b['name']) else 'Open-Air'
    delta_t = res['eff_temp_upper_surface_c'] - res['eff_temp_lower_cave_c']
    print(f"| **{b['name']}** | {ctx} | {b['lat']:.1f}°N, {b['lon']:.1f}°E | {b['age']:,} BP | {res['th_age_lower_cave']:,.0f} yr | {res['th_age_upper_surface']:,.0f} yr | **{res['seasonal_acceleration_factor']:.2f}×** | {res['eff_temp_lower_cave_c']:.1f} °C | {res['eff_temp_upper_surface_c']:.1f} °C | +{delta_t:.1f} °C |")

# Also check cave vs open-air statistics across the radiocarbon dataset
# Also check cave vs open-air statistics across the radiocarbon dataset
df = pd.read_parquet(r"D:\26 Modelling Collagen Hydrolysis\data\harmonized_c14_master_expanded.parquet")
collagen_df = df[df['material_category'].isin(['COLLAGEN', 'BONE_UNDIFFERENTIATED'])].copy()

# Cave indicator
cave_keywords = ['cave', 'grotte', 'hohle', 'höhle', 'grotta', 'cueva', 'abri', 'shelter', 'karst', 'cavern', 'spel']
def is_cave(row):
    text = f"{str(row.get('site_name', ''))} {str(row.get('material_raw', ''))}".lower()
    return any(k in text for k in cave_keywords)

collagen_df['is_cave'] = collagen_df.apply(is_cave, axis=1)

print("\n### Cave vs Open-Air Distribution Across Age Bins\n")
print("| Chronological Window | Total Collagen Determinations | Cave Determinations (*N*) | Open-Air Determinations (*N*) | Cave Proportion (%) |")
print("| :--- | ---: | ---: | ---: | ---: |")

windows = [
    ('< 5 ka BP (Late Holocene)', 0, 5000),
    ('5–15 ka BP (Early Holocene & Lateglacial)', 5000, 15000),
    ('15–30 ka BP (LGM & Upper Palaeolithic)', 15000, 30000),
    ('> 30 ka BP (Middle/Early Upper Palaeolithic)', 30000, 50000),
    ('Total Cohort (0–50 ka BP)', 0, 50000)
]

for label, min_age, max_age in windows:
    sub = collagen_df[(collagen_df['c14_age'] >= min_age) & (collagen_df['c14_age'] < max_age)]
    total = len(sub)
    n_cave = sub['is_cave'].sum()
    n_open = total - n_cave
    pct_cave = (n_cave / total * 100) if total > 0 else 0
    print(f"| **{label}** | {total:,} | {n_cave:,} | {n_open:,} | **{pct_cave:.1f}%** |")
