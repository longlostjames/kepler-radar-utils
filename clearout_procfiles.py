#!/usr/bin/env python

import os
import time
from datetime import datetime

# Get the current date in the form YYYYMMDD
current_date = datetime.now().strftime("%Y%m%d");

# Define the directory path labelled by the current date
procpath = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1/cobalt/L1b'
#procpath = '/gws/nopw/j04/ncas_radar_vol2/cjw/projects/cobalt/L1'
directory_path = os.path.join(procpath,current_date);

print(directory_path);

# Get the current time in seconds since the epoch
current_time = time.time();

# Define the time threshold (1 minute ago)
time_threshold = current_time - (1 * 60);

# Check if the directory exists
if os.path.exists(directory_path):
    # Iterate over each file in the directory
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename);
        # Check if it is a file (not a directory)
        if os.path.isfile(file_path):
            # Get the creation time of the file
            file_time = os.path.getctime(file_path);
            # If the file is older than 1 minute, delete it
            if file_time < time_threshold:
                print(filename.endswith('.nc'));
                if filename.endswith('.nc'):
                    os.remove(file_path);
                    print(f"Deleted: {file_path}");
                                                                                    
