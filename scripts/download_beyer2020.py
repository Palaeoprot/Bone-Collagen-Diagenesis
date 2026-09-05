import os
import urllib.request
import sys

dest_dir = r"D:\26 Modelling Collagen Hydrolysis\pastclim"
os.makedirs(dest_dir, exist_ok=True)
dest_file = os.path.join(dest_dir, "Beyer2020_annual_vars_v1.2.2.nc")

url = "https://zenodo.org/record/7388091/files/Beyer2020_annual_vars_v1.2.2.nc?download=1"

if os.path.exists(dest_file) and os.path.getsize(dest_file) > 320000000:
    print(f"File already exists and is complete: {dest_file} ({os.path.getsize(dest_file)} bytes)")
    sys.exit(0)

print(f"Downloading Beyer2020 NetCDF to {dest_file} (~323 MB)...")
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})

with urllib.request.urlopen(req) as resp, open(dest_file, 'wb') as out_f:
    total_size = int(resp.headers.get('Content-Length', 0))
    downloaded = 0
    chunk_size = 1024 * 1024 * 4 # 4MB chunks
    
    while True:
        chunk = resp.read(chunk_size)
        if not chunk:
            break
        out_f.write(chunk)
        downloaded += len(chunk)
        mb = downloaded / (1024 * 1024)
        pct = (downloaded / total_size * 100) if total_size else 0
        sys.stdout.write(f"\rDownloaded {mb:.1f} MB / {total_size/(1024*1024):.1f} MB ({pct:.1f}%)")
        sys.stdout.flush()

print(f"\nDownload completed successfully: {os.path.getsize(dest_file)} bytes")
