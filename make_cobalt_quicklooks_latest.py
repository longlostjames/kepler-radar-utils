#!/usr/bin/env python

import getopt, sys, os

import datetime

import netCDF4 as nc4

import pyart
import numpy as np
import numpy.ma as ma
import shutil
import glob
import gzip

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cmocean
import getpass, socket

import pandas as pd

import cftime

version = 0.1

from pathlib import Path
homepath = Path.home()


try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:bl", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

tracking_tag = 'AMOF_20231120125118';

campaign = 'cobalt';


yaml_project_file = os.path.join(homepath,'amof_campaigns',f'{campaign}_project.yml')
yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

blflag = False;

for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    elif o == "-b":
        blflag = True;
    elif o == "-l":
        latest = True
    else:
        assert False, "unhandled option"



ncas_radar_path = '/gws/nopw/j04/ncas_radar_vol2';
ncas_obs_proc_path = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1'
#inpath = os.path.join(ncas_radar_path,'cjw','projects',campaign,'L1');                                         
inpath = os.path.join(ncas_obs_proc_path,campaign,'L1b')

outpath = os.path.join(ncas_radar_path,'cjw','projects',campaign,'L1')
figpath = os.path.join(outpath,'quicklooks');


#inpath = os.path.join(ncas_radar_vol1_path,'cjw','projects',campaign,'kepler','L1a');
#figpath = os.path.join(ncas_radar_vol1_path,'public','ccrest_kepler');
cobalt_cmd_path = os.path.join(ncas_radar_path,'cjw','projects',campaign,'cobalt-command')

def read_cobalt_cmd_file(file_path):
    with open(file_path, "r") as file:
        line = file.readline()
        newline = line.rstrip().replace(': ', ':') 
        pairs = newline.split()
        data = {}

        for pair in pairs:
            key, value = pair.split(':')
            data[key.strip()] = value.strip()

    return data

import re
log_file_pattern = re.compile(r"kepler-(\d{8})-(\d{6})\.log")


#def find_nearest_file(directory, netcdf_time):
#    nearest_file = None
#    nearest_time_diff = datetime.timedelta.max
#
#    for file_name in os.listdir(directory):
#        match = log_file_pattern.match(file_name)
#        if match:
#            # Extract date and time from the file name
#            date_part = match.group(1)
#            time_part = match.group(2)
#            file_time_str = f"{date_part}{time_part}"
#            file_time = datetime.datetime.strptime(file_time_str, "%Y%m%d%H%M%S")

            # Calculate time difference
#            time_diff = abs(netcdf_time - file_time)

            # Update nearest file if this one is closer
#            if time_diff < nearest_time_diff:
#                nearest_file = file_name
#                nearest_time_diff = time_diff

#    return nearest_file


def extract_azimuth_from_file(file_path):  
    # Implement a way to extract the azimuth from the file.  
    # This is a placeholder; adjust according to your file's format.  
    data = read_cobalt_cmd_file(file_path)  
    return float(data.get('azimuth', 0))  # Assuming azimuth is stored under the key 'azimuth'  


def find_nearest_file(directory, netcdf_time, target_azimuth):  
    nearest_file = None  
    nearest_file_time = None  
    nearest_azimuth_diff = float('inf')  
    candidates = []  

    # First pass: find files with the closest azimuth  
    for file_name in os.listdir(directory):  
        match = log_file_pattern.match(file_name)  
        if match:  
            # Extract date and time from the file name  
            date_part = match.group(1)  
            time_part = match.group(2)  
            file_time_str = f"{date_part}{time_part}"  
            file_time = datetime.datetime.strptime(file_time_str, "%Y%m%d%H%M%S")  

            # Calculate time difference  
            time_diff = abs(netcdf_time - file_time)  

            # Only consider files with a time difference of less than 600 seconds  
            if time_diff < datetime.timedelta(seconds=600):  
                # Extract azimuth from the file  
                azimuth = extract_azimuth_from_file(os.path.join(directory, file_name))  
                azimuth_diff = abs(target_azimuth - azimuth)  

                # If this file's azimuth is closer, reset candidates  
                if azimuth_diff < nearest_azimuth_diff:  
                    nearest_azimuth_diff = azimuth_diff  
                    candidates = [(file_name, file_time, time_diff)]  
                elif azimuth_diff == nearest_azimuth_diff:  
                    candidates.append((file_name, file_time, time_diff))  

    # Second pass: among candidates, find the one with the closest time  
    nearest_time_diff = datetime.timedelta.max  

    for file_name, file_time, time_diff in candidates:  
        if time_diff < nearest_time_diff:  
            nearest_file = file_name  
            nearest_file_time = file_time  
            nearest_time_diff = time_diff  

    return nearest_file, nearest_file_time  # Return both the file name and file time  




def get_time_coverage_start(netcdf_file_path):
    try:
        # Open the NetCDF file
        with nc4.Dataset(netcdf_file_path, 'r') as DS:
            # Extract the global attribute 'time_coverage_start'
            if 'time_coverage_start' in DS.ncattrs():
                time_coverage_start = DS.getncattr('time_coverage_start')
                return time_coverage_start
            else:
                raise ValueError("Attribute 'time_coverage_start' not found in the NetCDF file.")
    except Exception as e:
        print(f"Error reading NetCDF file: {e}")
        return None

def get_target_altitude(earth_radius, target_range, zenith_angle):
    """
    Calculate the altitude of the radar target above the Earth's surface.

    Parameters:
        earth_radius (float): Radius of the Earth in meters.
        target_range (float): Range to the radar target in meters.
        zenith_angle (float): Angle between the local zenith and the target in radians.

    Returns:
        float: Altitude of the target in meters.
    """
    # Use the law of cosines to calculate the distance from Earth's center to the target
    r_t = np.sqrt(
        earth_radius**2 + target_range**2 + 2 * earth_radius * target_range * np.cos(zenith_angle)
    )


    # Altitude is the distance from Earth's surface to the target
    altitude = np.abs(r_t - earth_radius)

    return altitude


def great_circle_distance(earth_radius, target_range, zenith_angle):
    """
    Calculate the great-circle distance from the radar location to the sub-target point.

    Parameters:
        earth_radius (float): Radius of the Earth in meters.
        target_range (float): Range to the radar target in meters.
        zenith_angle (float): Angle between the local zenith and the target in radians.

    Returns:
        float: Great-circle distance in meters.
    """
    # Compute the cosine of the great-circle angular distance


    target_altitude = get_target_altitude(earth_radius, target_range, zenith_angle)

    # Calculate the angular distance in radians
    delta = np.arcsin(-target_range*np.sin(zenith_angle)/(earth_radius+target_altitude))

    # Convert angular distance to great-circle distance in meters
    great_circle_dist = earth_radius * delta

    return great_circle_dist




def make_cobalt_rhi_plot(ncfile,figpath,blflag=False):

    # Get the time_coverage_start
    time_coverage_start = get_time_coverage_start(ncfile);
    netcdf_time = datetime.datetime.strptime(time_coverage_start, "%Y-%m-%dT%H:%M:%SZ")

    print(netcdf_time);

    DS = nc4.Dataset(ncfile);
    scan_azimuth = DS['fixed_angle'][0];
    #print(f'product version = {product_version}')
    scan_span = np.abs(DS['elevation'][-1]-DS['elevation'][0])
    DS.close();

    print(f"scan_azimuth = {scan_azimuth}")
    print(f"scan_span = {scan_span}")


    if scan_span > 60:
        wind_scan = True
    else:
        wind_scan = False

    if not wind_scan:
        # Regular expression pattern for log files

        # Find the nearest file
        nearest_file, logfile_time = find_nearest_file(cobalt_cmd_path, netcdf_time, scan_azimuth)

        if nearest_file and logfile_time:
            logfile_time_str = datetime.datetime.strftime(logfile_time, '%Y-%m-%dT%H:%M:%SZ')
        else:
            print("No matching log file found.")
            logfile_time_str = "N/A"  # Or handle this case appropriately

        if nearest_file:
            print(f"The nearest file is: {nearest_file}")
        else:
            print("No matching files found.")

        cmd_data = read_cobalt_cmd_file(os.path.join(cobalt_cmd_path,nearest_file));

        print(cmd_data);

    if blflag:
        hmax = 4;
        xmin = 0;
        xmax = 20;
    else:
        hmax = 14;
        xmin = 0;
        xmax = 40;

    if wind_scan:
        hmax = 14;
        xmin = -14;
        xmax = 14;
        
    dbz_cmap = 'pyart_HomeyerRainbow';
    vel_cmap = 'pyart_balance';
    ldr_cmap = 'pyart_SpectralExtended';
    ldr_cmap = 'viridis';
    spw_cmap = 'pyart_SpectralExtended';
    #spw_cmap = 'pyart_ChaseSpectral';
    #spw_cmap = 'pyart_NWS_SPW';

    from matplotlib import colors

    DS = nc4.Dataset(ncfile);
    product_version = DS.product_version;
    #print(f'product version = {product_version}')
    DS.close();

    RadarDS = pyart.io.read_cfradial(ncfile);

    dtime0 = cftime.num2pydate(RadarDS.time['data'][0],RadarDS.time['units']);
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S");
    nsweeps = RadarDS.nsweeps;

    vel_field = RadarDS.fields['VEL']

    vel_limit_lower = vel_field['field_limit_lower']
    vel_limit_upper = vel_field['field_limit_upper']

#    fig, ax = plt.subplots(nsweeps,4,figsize=(24,nsweeps*4),constrained_layout=True)
#    fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)

    gatefilter = pyart.correct.GateFilter(RadarDS)
    gatefilter.exclude_below('SNR', -20)

    #pyart.correct.despeckle_field(RadarDS, "SNR", size=3, threshold=-20, gatefilter=gatefilter, delta=5.0)

    display = pyart.graph.RadarDisplay(RadarDS);

    figpath = os.path.join(figpath,'rhi',datestr);
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);



    if nsweeps>1:

        for s in range(nsweeps):
            print(f"sweep {s}/{nsweeps}");
            rhi_az = RadarDS.get_azimuth(s)[0];
            fig, ax = plt.subplots(2,2,figsize=(15,15),constrained_layout=True);
            fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)
            display.plot_rhi("DBZ", ax=ax[0,0], sweep=s, vmin=-60, vmax=40, norm=None, gatefilter=gatefilter,
                             cmap=dbz_cmap, colorbar_orient='horizontal',filter_transitions=True);         
            ax[0,0].set_ylim(0,hmax)
            ax[0,0].set_xlim(xmin,xmax)
            ax[0,0].grid(True)
            #ax[0,0].invert_xaxis();
            ax[0,0].set_aspect('equal');
#            display.plot_colorbar(ax=ax[0],shrink=0.8);
            display.plot_rhi("VEL", ax=ax[1,0], sweep=s, vmin=vel_limit_lower, vmax=vel_limit_upper, gatefilter=gatefilter,
                             norm=None, cmap=vel_cmap, colorbar_orient='horizontal',filter_transitions=True)
            ax[1,0].set_ylim(0,hmax)
            ax[1,0].set_xlim(xmin,xmax)
            ax[1,0].grid(True)
            #ax[1,0].invert_xaxis();
            ax[1,0].set_aspect('equal','box');
            display.plot_rhi("WIDTH", ax=ax[1,1], sweep=s, 
                             norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), gatefilter=gatefilter,
                             cmap=spw_cmap, colorbar_orient='horizontal',filter_transitions=True)
            ax[1,1].set_ylim(0,hmax)
            ax[1,1].set_xlim(xmin,xmax)
            ax[1,1].grid(True)
            #ax[1,1].invert_xaxis();
            ax[1,1].set_aspect('equal','box');
            display.plot_rhi("LDR", ax=ax[0,1], sweep=s, vmin=-35, vmax=5, norm=None,  gatefilter=gatefilter,
                             cmap=ldr_cmap, colorbar_orient='horizontal',filter_transitions=True)
            ax[0,1].set_ylim(0,hmax)
            ax[0,1].set_xlim(xmin,xmax)
            ax[0,1].grid(True)
            #ax[0,1].invert_xaxis();
            ax[0,1].set_aspect('equal','box');
            sweep_start_index = RadarDS.get_start(s);
            dtime_sweep = cftime.num2pydate(RadarDS.time['data'][sweep_start_index],RadarDS.time['units']);
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
            figname = f'ncas-radar-mobile-ka-band-1_chilbolton_cobalt_{dtime_sweep_str}_rhi_az{rhi_az:0.2f}_l1_{product_version}.png';
            plt.savefig(os.path.join(figpath,figname),dpi=300);
            plt.close();
    else:
        fig, ax = plt.subplots(2,2,figsize=(14,8),constrained_layout=True)
        fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)
        display.plot_rhi("DBZ", ax=ax[0,0], sweep=0, vmin=-60, vmax=40, 
                         norm=None, cmap=dbz_cmap, colorbar_orient='horizontal')
        ax[0,0].set_ylim(0,hmax)
        ax[0,0].set_xlim(xmin,xmax)
        ax[0,0].grid(True)
        #ax[0].invert_xaxis();
        ax[0,0].set_aspect('equal','box');
        display.plot_rhi("VEL", ax=ax[1,0], sweep=0, vmin=vel_limit_lower, vmax=vel_limit_upper,  gatefilter=gatefilter,
                         norm=None, cmap=vel_cmap, colorbar_orient='horizontal',filter_transitions=True)
        ax[1,0].set_ylim(0,hmax)
        ax[1,0].set_xlim(xmin,xmax)
        ax[1,0].grid(True)
        #ax[1].invert_xaxis();
        ax[1,0].set_aspect('equal','box');
        display.plot_rhi("WIDTH", ax=ax[0,1], sweep=0, gatefilter=gatefilter,
                         norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), 
                         cmap=spw_cmap, colorbar_orient='horizontal',filter_transitions=True)
        ax[0,1].set_ylim(0,hmax)
        ax[0,1].set_xlim(xmin,xmax)
        ax[0,1].grid(True)
        #ax[2].invert_xaxis();
        ax[0,1].set_aspect('equal','box'); 
        display.plot_rhi("LDR", ax=ax[1,1], sweep=0,vmin=-35, vmax=5, norm=None, gatefilter=gatefilter,
                         cmap=ldr_cmap, colorbar_orient='horizontal',filter_transitions=True)
        ax[1,1].set_ylim(0,hmax)
        ax[1,1].set_xlim(xmin,xmax)
        ax[1,1].grid(True)
        #ax[3].invert_xaxis();
        ax[1,1].set_aspect('equal','box');

        orig_title = ax[0,0].get_title();
        title_lines = orig_title.split("\n");

        if not wind_scan:
            target_range = float(cmd_data['range'])/1000.0;
            target_elev = np.deg2rad(float(cmd_data['elevation']));

            target_height = target_range * np.sin(target_elev);
            target_zenith_angle = np.pi/2.0-target_elev;

            print(f'Target zenith angle = {target_zenith_angle}')

            target_horiz_dist = target_range * np.cos(target_elev);
            target_height = target_range * np.sin(target_elev);

            target_arc_dist = great_circle_distance(6371.0, target_range, target_zenith_angle)
            target_altitude = get_target_altitude(6371.0, target_range, np.pi-target_zenith_angle)

            print(f'Target arc dist = {target_arc_dist}');
            print(f'Target altitude = {target_altitude}');
            print(f'Target horiz dist = {target_horiz_dist}');
            print(f'Target height = {target_height}')

            ax[0,0].plot(target_arc_dist, target_altitude, 'k+')  
            ax[1,0].plot(target_arc_dist, target_altitude, 'k+')  
            ax[0,1].plot(target_arc_dist, target_altitude, 'k+')  
            ax[1,1].plot(target_arc_dist, target_altitude, 'k+')  

            ax[0,0].plot(target_horiz_dist, target_height, 'b+')  
            ax[1,0].plot(target_horiz_dist, target_height, 'b+')  
            ax[0,1].plot(target_horiz_dist, target_height, 'b+')  
            ax[1,1].plot(target_horiz_dist, target_height, 'b+')  

            target_str = f"{logfile_time_str} {cmd_data['flight_id']} {cmd_data['icao_id']} {cmd_data['aircraft_type']} lean_burn={cmd_data['lean_burn']} Az={float(cmd_data['azimuth']):.1f} El={float(cmd_data['elevation']):.1f}"
            if 'advected' in cmd_data:
                target_str += f" advected {cmd_data['advected']}"

            title_lines[1] = target_str;
        else:
            title_lines[1] = "Wind Scan"
        
        new_title = "\n".join(title_lines)  # Rejoin the lines

        ax[0,0].set_title("",fontsize=8)
        ax[1,0].set_title("",fontsize=8)
        ax[0,1].set_title("",fontsize=8)
        ax[1,1].set_title("",fontsize=8)

        fig.suptitle(new_title, fontsize=11)  

        figname = f'ncas-mobile-ka-band-radar-1_chilbolton_cobalt_{dtime0_str}_l1_{product_version}.png'
        plt.savefig(os.path.join(figpath,figname),dpi=300);
        plt.close();



def make_cobalt_vpt_plot_day(datestr,inpath,figpath,blflag=False):
    
    if blflag:
        hmax = 4;
    else:
        hmax = 14;

    dbz_cmap = 'pyart_HomeyerRainbow';
    vel_cmap = 'pyart_balance';
    ldr_cmap = 'pyart_SpectralExtended';
    ldr_cmap = 'viridis';
    spw_cmap = 'pyart_SpectralExtended';
    #spw_cmap = 'pyart_ChaseSpectral';
    #spw_cmap = 'pyart_NWS_SPW';

    velmin = -5.0;
    velmax = 5.0;

    from matplotlib import colors


    current_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    prev_date = current_date - datetime.timedelta(days=1);
    prevstr = prev_date.strftime('%Y%m%d');

    inpath_date = os.path.join(inpath,datestr);
    
    os.chdir(inpath_date);
    vpt_file = [os.path.join(inpath_date,f) for f in glob.glob('*{}*vpt*.nc'.format(datestr))][0]


    DS = nc4.Dataset(vpt_file);
    product_version = DS.product_version;
    print(f'product version = {product_version}')
    DS.close();


    RadarDS_VPT = pyart.io.read_cfradial(vpt_file);
    nsweeps = RadarDS_VPT.nsweeps;
    

    vel_field = RadarDS_VPT.fields['VEL']

    vel_limit_lower = vel_field['field_limit_lower']
    vel_limit_upper = vel_field['field_limit_upper']

    nsweeps_prev =0;
    try:
        inpath_prev = os.path.join(inpath,prevstr);
        os.chdir(inpath_prev);

        vpt_file_prev = [os.path.join(inpath_prev,f) for f in glob.glob('*{}*vpt*.nc'.format(prevstr))][0]

        RadarDS_VPT_prev = pyart.io.read_cfradial(vpt_file_prev);
        nsweeps_prev = RadarDS_VPT_prev.nsweeps;
    except:
        pass
    
    fig, ax = plt.subplots(4,1,figsize=(12,18),constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)

    
    #pyart.correct.despeckle_field(RadarDS, "SNR", size=3, threshold=-20, gatefilter=gatefilter, delta=5.0)

    #display = pyart.graph.RadarDisplay(RadarDS_VPT);
    RadarSweepDS = RadarDS_VPT.extract_sweeps([0])
    print("sweep = 0");
    gatefilter = pyart.correct.GateFilter(RadarSweepDS)
    #gatefilter.exclude_below('SNR', -20)
    display = pyart.graph.RadarDisplay(RadarSweepDS);

    import cftime
    dtime = cftime.num2pydate(RadarDS_VPT.time['data'],RadarDS_VPT.time['units'])

    dt_min = dtime[0].replace(hour=0,minute=0,second=0); #datetime.strptime(datestr, '%Y%m%d')
    dt_max = dt_min + datetime.timedelta(days=1)

    time_str = dtime[0].strftime("%Y-%m-%d");
    
   
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "DBZ")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter,
                     vmin=-60, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, colorbar_orient='horizontal')
    ax[0].set_ylim(0,hmax)
    ax[0].grid(True);
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "VEL")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter,
                     vmin=vel_limit_lower, vmax=vel_limit_upper, norm=None, cmap=vel_cmap, colorbar_orient='horizontal')
    ax[1].set_ylim(0,hmax)
    ax[1].grid(True)
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "WIDTH")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter, 
                     norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), cmap=spw_cmap, colorbar_orient='horizontal')
    ax[2].set_ylim(0,hmax)
    ax[2].grid(True)
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "LDR")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter,
                     vmin=-35, vmax=5, norm=None, cmap=ldr_cmap, colorbar_orient='horizontal')
    ax[3].set_ylim(0,hmax)
    ax[3].grid(True)

    ax[0].set_xlim(dt_min,dt_max);
    ax[1].set_xlim(dt_min,dt_max);
    ax[2].set_xlim(dt_min,dt_max);
    ax[3].set_xlim(dt_min,dt_max);

    ax[0].set_xlabel('Time (UTC)');
    ax[1].set_xlabel('Time (UTC)');
    ax[2].set_xlabel('Time (UTC)');
    ax[3].set_xlabel('Time (UTC)');

    nsweeps4now = nsweeps-1;

    for s in range(1,nsweeps4now):
        RadarSweepDS = RadarDS_VPT.extract_sweeps([s])
        print(f"sweep = {s}");
        gatefilter = pyart.correct.GateFilter(RadarSweepDS)
        gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(RadarSweepDS);
        display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                         vmin=-60, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, 
                         colorbar_flag=False, title_flag=False)
        display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                         vmin=vel_limit_lower, vmax=vel_limit_upper, norm=None, cmap=vel_cmap, colorbar_flag=False, title_flag=False)
        display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, edges=False, gatefilter=gatefilter, 
                         norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), cmap=spw_cmap, 
                         colorbar_flag=False, title_flag=False)
        display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                         vmin=-35, vmax=5, norm=None, cmap=ldr_cmap, colorbar_flag=False, title_flag=False)

    if (nsweeps_prev>0):
        for s in range(nsweeps_prev-1,nsweeps_prev):
            RadarSweepDS = RadarDS_VPT_prev.extract_sweeps([s])
            print(f"sweep = {s}");
            gatefilter = pyart.correct.GateFilter(RadarSweepDS)
            gatefilter.exclude_below('SNR', -20)
            display = pyart.graph.RadarDisplay(RadarSweepDS);
            display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                             vmin=-60, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, 
                             colorbar_flag=False, title_flag=False)
            display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                             vmin=vel_limit_lower, vmax=vel_limit_upper, norm=None, cmap=vel_cmap, colorbar_flag=False, title_flag=False)
            display.plot_vpt("WIDTH", ax=ax[2], time_axis_flag=True, edges=False, gatefilter=gatefilter, 
                             norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), cmap=spw_cmap, 
                             colorbar_flag=False, title_flag=False)
            display.plot_vpt("LDR", ax=ax[3], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                             vmin=-35, vmax=5, norm=None, cmap=ldr_cmap, colorbar_flag=False, title_flag=False)
     
    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)
    ax[3].grid(True)

    ax[0].set_xlabel('Time (UTC)');
    ax[1].set_xlabel('Time (UTC)');
    ax[2].set_xlabel('Time (UTC)');
    ax[3].set_xlabel('Time (UTC)');

    figname = f'ncas-mobile-ka-band-radar-1_chilbolton_cobalt_{datestr}_vpt_l1_{product_version}.png';
    
    if blflag:
        figname = figname.replace('.png','_bl.png');

    figpath = os.path.join(figpath,'vpt');
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    plt.savefig(os.path.join(figpath,figname),dpi=300);

    plt.close();

def make_ccrest_vad_plot_day(datestr,inpath,figpath,zlevels,blflag=False):

    if blflag:
        hmax = 4;
    else:
        hmax = 12;

    inpath_date = os.path.join(inpath,datestr);

    os.chdir(inpath_date);
    vad_file = [os.path.join(inpath_date,f) for f in glob.glob('*{}*vad*.nc'.format(datestr))][0];

    DS = nc4.Dataset(vad_file);
    product_version = DS.product_version;
    print(f'product version = {product_version}')
    DS.close();

    radar = pyart.io.read_cfradial(vad_file);

    vad_ray_index_start = [];
    vad_ray_index_end = [];

    for s in range(radar.nsweeps):
        vad_ray_index_start.append(radar.sweep_start_ray_index['data'][s]);
        vad_ray_index_end.append(radar.sweep_end_ray_index['data'][s]);
    
    dt_vad_start = cftime.num2pydate(radar.time['data'][vad_ray_index_start[:]],radar.time['units']);
    dt_vad_end = cftime.num2pydate(radar.time['data'][vad_ray_index_end[:]],radar.time['units']);

    u_allsweeps = []
    v_allsweeps = []

    for s in range(radar.nsweeps):
        radar_1sweep = radar.extract_sweeps([s])
        #vel_texture = pyart.retrieve.calculate_velocity_texture(radar_1sweep, vel_field='VEL', wind_size=6, nyq=4.7, check_nyq_uniform=True)
        #radar_1sweep.add_field('txtVEL',vel_texture)
        gatefilter = pyart.correct.GateFilter(radar_1sweep)
        #gatefilter.exclude_below('SNR', -5.8)
        gatefilter.exclude_below('SNR', 0)
        vad = pyart.retrieve.vad_browning(radar_1sweep, "VEL", z_want=zlevels,gatefilter=gatefilter) 
        u_allsweeps.append(vad.u_wind)
        v_allsweeps.append(vad.v_wind)

    u_vel = np.array(u_allsweeps);
    v_vel = np.array(v_allsweeps);
    # Average U and V over all sweeps and compute magnitude and angle
    #u_avg = np.nanmean(np.array(u_allsweeps), axis=0)
    #v_avg = np.nanmean(np.array(v_allsweeps), axis=0)
    #orientation = np.rad2deg(np.arctan2(-u_avg, -v_avg)) % 360
    #speed = np.sqrt(u_avg**2 + v_avg**2)
    orientation = np.rad2deg(np.arctan2(-u_vel, -v_vel)) % 360;
    speed = np.sqrt(u_vel**2 + v_vel**2)

    speed = ma.masked_where(speed>100.,speed);
    orientation = ma.masked_where(speed>100.,orientation);

    vad_duration = (dt_vad_end-dt_vad_start);
    dt_vad_mid = dt_vad_start + 0.5*vad_duration; 
    vad_duration[:] = datetime.timedelta(minutes=12);

    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%H:%M');

    fig, ax = plt.subplots(2,1,figsize=(12,8),constrained_layout=True);
    fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2);

    dtime = cftime.num2pydate(radar.time['data'],radar.time['units'])

    dt_min = dtime[0].replace(hour=0,minute=0,second=0); #datetime.strptime(datestr, '%Y%m%d')
    dt_max = dt_min + datetime.timedelta(days=1)

    ax[0].xaxis.set_major_formatter(myFmt);
    ax[0].grid(True);
    h1 = ax[0].pcolormesh([dt_vad_start[0]-0.5*vad_duration[0],dt_vad_end[0]+0.5*vad_duration[0]],
                          zlevels/1000.,orientation[0:1,:-1].transpose(),cmap='twilight_shifted',vmin=0,vmax=360);
    ax[0].set_ylim(0,12);
    ax[0].set_ylabel('Distance above radar [km]');
    cb0 = plt.colorbar(h1,ax=ax[0],orientation='horizontal',shrink=0.8);
    cb0.ax.set_xlabel("Wind from direction (deg)");
    ax[0].set_facecolor('white');
    for s in range(1,radar.nsweeps):
        ax[0].pcolormesh([dt_vad_start[s]-0.5*vad_duration[s],dt_vad_end[s]+0.5*vad_duration[s]],
                         zlevels/1000.,orientation[s:s+1,:-1].transpose(),cmap='twilight_shifted',vmin=0,vmax=360);

    ax[0].set_xlabel("Time (UTC)");

    ax[1].xaxis.set_major_formatter(myFmt);
    ax[1].grid(True);
    h2 = ax[1].pcolormesh([dt_vad_start[0]-0.5*vad_duration[0],dt_vad_end[0]+0.5*vad_duration[0]],
                          zlevels/1000.,speed[0:1,:-1].transpose(),cmap='viridis',vmin=0,vmax=50);
    ax[1].grid(True);
    ax[1].set_ylim(0,12);
    ax[1].set_ylabel('Distance above radar (km)');
    cb1 = plt.colorbar(h2,ax=ax[1],orientation='horizontal',shrink=0.8);
    cb1.ax.set_xlabel("Wind speed (m/s)");
    ax[1].set_facecolor('gainsboro');


    for s in range(1,radar.nsweeps):
        ax[1].pcolormesh([dt_vad_start[s]-0.5*vad_duration[s],dt_vad_end[s]+0.5*vad_duration[s]],zlevels/1000.,speed[s:s+1,:-1].transpose(),cmap='viridis',vmin=0,vmax=50);

    nlevels = zlevels.shape[0];
    nsweeps = radar.nsweeps

    x = np.tile(dt_vad_mid,[nlevels,1]).transpose();
    y = np.tile(zlevels/1000.,[nsweeps,1]);
    print(x.shape,y.shape)
#    ax[1].barbs(x[:,5::10],y[:,5::10],u_vel[:,5::10], v_vel[:,5::10],sizes=dict(emptybarb=0.),length=6,barbcolor='sandybrown')
    ax[1].barbs(x[:,5::10],y[:,5::10],u_vel[:,5::10], v_vel[:,5::10],np.sqrt(u_vel[:,5::10]**2+v_vel[:,5::10]),sizes=dict(emptybarb=0.),length=6,cmap='Grays',clim=[0,50])
    ax[1].set_facecolor('gainsboro');
    ax[1].set_xlabel("Time (UTC)");

    dtime = cftime.num2pydate(radar.time['data'],radar.time['units'])[0];
    time_str = dtime.strftime("%Y-%m-%d");
    title0 =  f"{pyart.graph.common.generate_radar_name(radar)} {time_str}" + "\n"; 
    title0 += f"Wind direction from VAD at elevation {radar.fixed_angle['data'][0]}";
    title1 =  f"{pyart.graph.common.generate_radar_name(radar)} {time_str}" + "\n"; 
    title1 += f"Wind speed from VAD at elevation {radar.fixed_angle['data'][0]}";

    ax[0].set_xlim(dt_min,dt_max);
    ax[1].set_xlim(dt_min,dt_max);
    ax[0].grid(True);
    ax[1].grid(True);

    ax[0].set_title(title0);
    ax[1].set_title(title1);

    figpath = os.path.join(figpath,'vad');
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    plt.savefig(os.path.join(figpath,f'ncas-radar-mobile-ka-band-1_chilbolton_ccrest-m_{datestr}_vad_l1_{product_version}.png'),dpi=300)

    plt.close();



def make_cobalt_rhi_plots_day(datestr,inpath,figpath,blflag=False):
    inpath_date = os.path.join(inpath,datestr);
    
    os.chdir(inpath_date);
    rhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*rhi*.nc'.format(datestr))]

    print(f'rhi files = ',rhi_files);

    for f in rhi_files:
        try:
            make_cobalt_rhi_plot(f,figpath,blflag=blflag);
        except:
            print('Exception in RHI plotting')
        finally:
            pass

    return

make_cobalt_rhi_plots_day(datestr,inpath,figpath,blflag=False);

try:
    make_cobalt_vpt_plot_day(datestr,inpath,figpath,blflag=False);
finally:
    pass


#make_ccrest2_rhi_plots_day(datestr,inpath,figpath,blflag=False);
#zlevels = np.arange(100, 15000, 100);  # height above radar
#make_ccrest_vad_plot_day(datestr,inpath,figpath,zlevels);

