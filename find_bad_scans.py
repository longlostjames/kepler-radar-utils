import os
import argparse
from netCDF4 import Dataset

def find_short_time_dim_files(directory, threshold=15):
    short_files = []

    for fname in os.listdir(directory):
        if fname.endswith('.nc'):
            path = os.path.join(directory, fname)
            try:
                with Dataset(path, 'r') as nc:
                    if 'time' in nc.dimensions:
                        time_len = len(nc.dimensions['time'])
                        if time_len < threshold:
                            short_files.append((fname, time_len))
                    else:
                        print(f"⚠️  No 'time' dimension in: {fname}")
            except Exception as e:
                print(f"❌ Error reading {fname}: {e}")
    
    return short_files

def main():
    parser = argparse.ArgumentParser(description="Find NetCDF files with short 'time' dimensions.")
    parser.add_argument("-d", "--date", required=True, help="Date in YYYYMMDD format")
    args = parser.parse_args()

    # Construct the directory path based on the date
    directory = f'/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/cobalt/L1b/{args.date}'
    
    if not os.path.isdir(directory):
        print(f"❌ Directory does not exist: {directory}")
        return

    short_time_files = find_short_time_dim_files(directory)

    print("\nFiles with 'time' dimension < 15:")
    for fname, tlen in short_time_files:
        print(f"{fname} — time length: {tlen}")

if __name__ == "__main__":
    main()