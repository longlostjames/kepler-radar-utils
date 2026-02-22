import os
import shutil
import re
import pandas as pd

def read_bad_scans(file_path):
    bad_scans = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            match = re.match(r'(.+\.nc)\s+—\s+time length: (\d+)', line)
            if match:
                filename = match.group(1)
                time_len = int(match.group(2))
                bad_scans.append({
                    'filename': filename,
                    'time_length': time_len
                })
            else:
                print(f"⚠️ Could not parse line: {line}")
    
    return pd.DataFrame(bad_scans)

def move_bad_files(bad_df, search_root, target_dir='BAD_DATA'):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    
    moved = []

    for fname in bad_df['filename']:
        found = False
        for root, dirs, files in os.walk(search_root):
            if fname in files:
                src = os.path.join(root, fname)
                dst = os.path.join(target_dir, fname)
                shutil.move(src, dst)
                moved.append(fname)
                print(f"✅ Moved: {fname}")
                found = True
                break
        if not found:
            print(f"❌ File not found: {fname}")
    
    return moved

# Example usage:
bad_df = read_bad_scans('BAD_SCANS')
search_root = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/cobalt/L1b'
moved_files = move_bad_files(bad_df, search_root)