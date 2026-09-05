import os
import sys
import ssl
import time
import urllib.request

DATA_DIR = r"D:\26 Modelling Collagen Hydrolysis\data"
RAW_DIR = os.path.join(DATA_DIR, "raw_sources")
os.makedirs(RAW_DIR, exist_ok=True)

SOURCES = [
    {
        "name": "14cpalaeolithic",
        "filename": "14cpalaeolithic_v33.xlsx",
        "url": "https://ees.kuleuven.be/en/geography/projects/14c-palaeolithic/download/radiocarbon-palaeolithic-europe-database-v33-extract.xlsx"
    },
    {
        "name": "calpal",
        "filename": "CalPal_2020_08_20.tsv",
        "url": "https://raw.githubusercontent.com/nevrome/CalPal-Database/master/CalPal_2020_08_20.tsv"
    },
    {
        "name": "radonb",
        "filename": "radonb_daily.txt",
        "url": "https://radonb.ufg.uni-kiel.de/radondaily.txt",
        "insecure_ssl": True
    },
    {
        "name": "rxpand_south_america",
        "filename": "rxpand_south_america.csv",
        "url": "https://raw.githubusercontent.com/jgregoriods/rxpand/master/data/expand.csv"
    },
    {
        "name": "austarch",
        "filename": "austarch_australia.csv",
        "url": "https://archaeologydataservice.ac.uk/catalogue/adsdata/arch-1661-1/dissemination/csv/Austarch_1-3_and_IDASQ_28Nov13-1.csv"
    },
    {
        "name": "14sea",
        "filename": "14sea_southeast_asia.xlsx",
        "url": "http://www.14sea.org/img/14SEA_Full_Dataset_2017-01-29.xlsx"
    },
    {
        "name": "euroevol",
        "filename": "euroevol_samples.csv",
        "url": "https://discovery.ucl.ac.uk/id/eprint/1469811/7/EUROEVOL09-07-201516-34_C14Samples.csv"
    },
    {
        "name": "sard_south_africa",
        "filename": "sard_south_africa.csv",
        "url": "https://raw.githubusercontent.com/emmaloftus/Southern-African-Radiocarbon-Database/main/SARD_Mar2021_14C.csv"
    },
    {
        "name": "aida",
        "filename": "aida_near_east.csv",
        "url": "https://raw.githubusercontent.com/apalmisano82/AIDA/main/dates.csv"
    },
    {
        "name": "nerd",
        "filename": "nerd_near_east.csv",
        "url": "https://raw.githubusercontent.com/apalmisano82/NERD/main/nerd.csv"
    },
    {
        "name": "adrac",
        "filename": "adrac_central_africa.csv",
        "url": "https://raw.githubusercontent.com/dirkseidensticker/aDRAC/master/aDRAC.csv"
    },
    {
        "name": "caribbean",
        "filename": "caribbean_14c.csv",
        "url": "https://raw.githubusercontent.com/philriris/caribbean-14C/main/data/caribbean_14C.csv"
    },
    {
        "name": "jomon",
        "filename": "jomon_japan.csv",
        "url": "https://raw.githubusercontent.com/ercrema/jomonPhasesAndPopulation/master/data/c14dates.csv"
    },
    {
        "name": "p3k14c_latest",
        "filename": "p3k14c_scrubbed_fuzzed.csv",
        "url": "https://raw.githubusercontent.com/people3k/p3k14c/refs/tags/2025.07/inst/p3k14c_scrubbed_fuzzed.csv"
    },
    {
        "name": "xronos_csv",
        "filename": "xronos_data.csv",
        "url": "https://xronos.ch/data.csv",
        "retries": 3,
        "timeout": 30
    }
]

def download_file(src):
    dest_path = os.path.join(RAW_DIR, src["filename"])
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 1000:
        print(f"[CACHED] {src['name']} -> {dest_path} ({os.path.getsize(dest_path)} bytes)")
        return dest_path

    url = src["url"]
    retries = src.get("retries", 2)
    timeout = src.get("timeout", 15)
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Palaeoproteomics Research Crawler"}
    
    ctx = None
    if src.get("insecure_ssl"):
        ctx = ssl._create_unverified_context()

    print(f"[DOWNLOADING] {src['name']} from {url}...")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                data = resp.read()
                with open(dest_path, "wb") as f:
                    f.write(data)
                print(f"[SUCCESS] {src['name']} -> Saved {len(data):,} bytes to {dest_path}")
                return dest_path
        except Exception as e:
            print(f"[ATTEMPT {attempt} FAILED] {src['name']}: {e}")
            if attempt < retries:
                time.sleep(2)
    print(f"[WARNING] Could not download {src['name']}. Will skip or proceed with mirrors.")
    return None

if __name__ == "__main__":
    print(f"=== Downloading Radiocarbon Datasets into {RAW_DIR} ===")
    results = {}
    for src in SOURCES:
        path = download_file(src)
        results[src["name"]] = path is not None
    
    print("\n=== Download Summary ===")
    for k, v in results.items():
        status = "OK" if v else "FAILED/SKIPPED"
        print(f"  {k:25s}: {status}")
