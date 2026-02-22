#!/usr/bin/env python

import getopt, sys, os

import datetime

import netCDF4 as nc4


import shutil
import glob
import gzip

import getpass, socket


sys.path.append('/home/users/cjwalden/git/kepler-radar-utils')

import kepler_utils_old

from pathlib import Path
homepath = Path.home()


try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:", ["date=","inpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)


campaign = 'ccrest-m';

data_version = "1.0.0"


for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    else:
        assert False, "unhandled option"

amof_proc_path = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-mobile-ka-band-1'

procpath = os.path.join(amof_proc_path,campaign,'L1b',datestr);


def find_netcdf_files(root_directory):  
    netcdf_files = []  
    for dirpath, dirnames, filenames in os.walk(root_directory):  
        for filename in filenames:  
            if filename.endswith('.nc'):  # NetCDF files usually have .nc extension  
                full_path = os.path.join(dirpath, filename)  
                netcdf_files.append(full_path)  
    
    return netcdf_files  

cfradfiles = find_netcdf_files(procpath)  
    
print("NetCDF files found:")  
for f in cfradfiles:  
    print(f)
    kepler_utils_old.cfradial_add_processing_software_metadata(f)


