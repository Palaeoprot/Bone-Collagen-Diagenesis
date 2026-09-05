import xarray as xr
import numpy as np
import pandas as pd
import os

class PaleoclimateEngine:
    def __init__(self, nc_path: str = r"D:\26 Modelling Collagen Hydrolysis\pastclim\Beyer2020_annual_vars_v1.2.2.nc"):
        if not os.path.exists(nc_path):
            raise FileNotFoundError(f"Beyer2020 NetCDF not found at: {nc_path}")
        self.nc_path = nc_path
        # Must open with decode_times=False because time units are 'years since 1950-01-01'
        self.ds = xr.open_dataset(nc_path, decode_times=False)
        
        # In Beyer2020, time is stored as negative BP values (e.g. -120000 to 0)
        # We convert to positive years BP (0 to 120000 BP)
        raw_times = self.ds["time"].values
        self.times_bp = np.abs(raw_times)
        
        # Sort so time runs from 0 BP (present) to 120000 BP (past)
        sort_idx = np.argsort(self.times_bp)
        self.times_bp = self.times_bp[sort_idx]
        self.ds_sorted = self.ds.isel(time=sort_idx)
        
        # Pre-extract coordinate arrays for fast indexing
        self.lats = self.ds["latitude"].values
        self.lons = self.ds["longitude"].values

    def get_temperature_series(self, lat: float, lon: float):
        """
        Extract mean annual temperature series (in Celsius and Kelvin) for a given lat/lon.
        Handles spatial nearest-neighbor search.
        """
        # Longitude normalization to [-180, 180]
        if lon > 180.0:
            lon -= 360.0
            
        pt = self.ds_sorted["bio01"].sel(latitude=lat, longitude=lon, method="nearest")
        temps_c = pt.values.astype(float)
        
        # If oceanic/masked NaN, search 3x3 surrounding box
        if np.isnan(temps_c).any():
            lat_idx = np.argmin(np.abs(self.lats - lat))
            lon_idx = np.argmin(np.abs(self.lons - lon))
            sub = self.ds_sorted["bio01"][:, max(0, lat_idx-2):lat_idx+3, max(0, lon_idx-2):lon_idx+3].values
            # Mean across non-nan spatial neighbors for each time slice
            fallback_c = np.nanmean(sub, axis=(1, 2))
            if not np.isnan(fallback_c).all():
                temps_c = np.where(np.isnan(temps_c), fallback_c, temps_c)

        temps_k = temps_c + 273.15
        return self.times_bp, temps_c, temps_k

    def calculate_thermal_age(self, lat: float, lon: float, cal_age_bp: float,
                              ea_kj: float = 173.0, t_ref_c: float = 10.0, delta_t_micro: float = 0.0):
        """
        Numerically integrate effective thermal age (normalised to t_ref_c):
        Age_therm = \int_0^{t_cal} exp(-Ea/R * (1/(T(t) + delta_t) - 1/T_ref)) dt
        """
        times_bp, temps_c, temps_k = self.get_temperature_series(lat, lon)
        
        # Temperature with microclimate adjustment
        t_hist_k = temps_k + delta_t_micro
        t_ref_k = t_ref_c + 273.15
        r_gas = 8.314462618e-3 # kJ / (mol * K)
        
        # Bound to [0, cal_age_bp]
        mask = times_bp <= cal_age_bp
        sub_times = times_bp[mask].tolist()
        sub_temps = t_hist_k[mask].tolist()
        
        # Ensure exact end point at cal_age_bp
        if len(sub_times) == 0 or sub_times[-1] < cal_age_bp:
            # Interpolate temperature at cal_age_bp
            t_interp = np.interp(cal_age_bp, times_bp, t_hist_k)
            sub_times.append(cal_age_bp)
            sub_temps.append(t_interp)

        sub_times = np.array(sub_times)
        sub_temps = np.array(sub_temps)
        
        # Clean any remaining NaNs
        if np.isnan(sub_temps).any():
            valid_mean = np.nanmean(sub_temps)
            if np.isnan(valid_mean):
                valid_mean = 273.15 + 10.0 # safe fallback
            sub_temps = np.nan_to_num(sub_temps, nan=valid_mean)

        # Arrhenius acceleration factor at each time step
        arrh_factors = np.exp(-(ea_kj / r_gas) * (1.0 / sub_temps - 1.0 / t_ref_k))
        
        # Trapezoidal integration (NumPy 2.x compatible)
        if len(sub_times) >= 2:
            if hasattr(np, "trapezoid"):
                thermal_age = np.trapezoid(arrh_factors, sub_times)
            else:
                thermal_age = np.trapz(arrh_factors, sub_times)
        else:
            thermal_age = arrh_factors[0] * cal_age_bp
            
        return float(thermal_age)

if __name__ == "__main__":
    engine = PaleoclimateEngine()
    
    # Smoketest on benchmark archaeological sites:
    benchmarks = [
        {"name": "Denisova Cave (Altai)", "lat": 51.3975, "lon": 84.6761, "age": 45000},
        {"name": "Vindija Cave (Croatia)", "lat": 46.2992, "lon": 16.0711, "age": 38000},
        {"name": "Sunghir (Russia)", "lat": 56.1756, "lon": 40.5058, "age": 34000},
        {"name": "Grotta Guattari (Italy)", "lat": 41.2286, "lon": 13.0903, "age": 50000},
        {"name": "Hayonim Cave (Levant)", "lat": 32.9233, "lon": 35.2167, "age": 28000}
    ]
    
    print("\n=== Benchmark Sites Paleoclimate & Thermal Age Smoketest (Ea = 173 kJ/mol, T_ref = 10°C) ===")
    for b in benchmarks:
        t_bp, tc, tk = engine.get_temperature_series(b["lat"], b["lon"])
        th_age_173 = engine.calculate_thermal_age(b["lat"], b["lon"], b["age"], ea_kj=173.0)
        th_age_100 = engine.calculate_thermal_age(b["lat"], b["lon"], b["age"], ea_kj=100.0)
        
        print(f"\nSite: {b['name']}")
        print(f"  Coordinates: ({b['lat']:.2f}, {b['lon']:.2f}) | Calendar Age: {b['age']:,} BP")
        print(f"  Modern MAT: {tc[0]:.2f} °C | LGM (~21 ka) MAT: {tc[t_bp==21000][0]:.2f} °C")
        print(f"  Thermal Age @ 173 kJ/mol: {th_age_173:,.1f} y (Relative to 10°C)")
        print(f"  Thermal Age @ 100 kJ/mol: {th_age_100:,.1f} y")
