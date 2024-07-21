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
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:b", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

tracking_tag = 'AMOF_20230201132601';

campaign = 'ccrest-m';


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
    else:
        assert False, "unhandled option"



ncas_radar_vol1_path = '/gws/nopw/j04/ncas_radar_vol1';
ncas_obs_proc_path = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-mobile-ka-band-1'
#inpath = os.path.join(ncas_radar_vol1_path,'cjw','projects',campaign,'kepler','L1a');                                         
inpath = os.path.join(ncas_obs_proc_path,campaign,'L1b')

figpath = os.path.join(inpath,'quicklooks');


#inpath = os.path.join(ncas_radar_vol1_path,'cjw','projects',campaign,'kepler','L1a');
#figpath = os.path.join(ncas_radar_vol1_path,'public','ccrest_kepler');


def make_ccrest_rhi_plot(ncfile,figpath,ccrest_az,blflag=False):

    if blflag:
        hmax = 4;
        xmin = 0;
        xmax = 20;
    else:
        hmax = 12;
        xmin = 0;
        xmax = 25;
    
    dbz_cmap = 'pyart_HomeyerRainbow';
    vel_cmap = 'pyart_balance';
    ldr_cmap = 'pyart_SpectralExtended';
    #ldr_cmap = 'pyart_ChaseSpectral';
    spw_cmap = 'pyart_SpectralExtended';
    #spw_cmap = 'pyart_ChaseSpectral';
    #spw_cmap = 'pyart_NWS_SPW';

    from matplotlib import colors

    DS = nc4.Dataset(ncfile);
    product_version = DS.product_version;
    print(f'product version = {product_version}')
    DS.close();

    RadarDS = pyart.io.read_cfradial(ncfile);

    dtime0 = cftime.num2pydate(RadarDS.time['data'][0],RadarDS.time['units']);
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S");
    nsweeps = RadarDS.nsweeps;

#    fig, ax = plt.subplots(nsweeps,4,figsize=(24,nsweeps*4),constrained_layout=True)
#    fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)

    gatefilter = pyart.correct.GateFilter(RadarDS)
    gatefilter.exclude_below('SNR', -8)

    #pyart.correct.despeckle_field(RadarDS, "SNR", size=3, threshold=-20, gatefilter=gatefilter, delta=5.0)

    display = pyart.graph.RadarDisplay(RadarDS);

    figpath = os.path.join(figpath,'rhi',datestr);
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    figpath_ccrest_az = os.path.join(figpath,ccrest_az);
    if not os.path.isdir(figpath_ccrest_az):
        os.makedirs(figpath_ccrest_az);


    if nsweeps>1:

        for s in range(nsweeps):
            print(f"sweep {s}/{nsweeps}");
            rhi_az = RadarDS.get_azimuth(s)[0];
            fig, ax = plt.subplots(2,2,figsize=(15,15)); #,constrained_layout=True)
            #fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)
            display.plot_rhi("DBZ", ax=ax[0,0], sweep=s, vmin=-40, vmax=40, norm=None, gatefilter=gatefilter,
                             cmap=dbz_cmap, colorbar_orient='horizontal',reverse_xaxis=True,filter_transitions=True);         
            ax[0,0].set_ylim(0,hmax)
            ax[0,0].set_xlim(xmin,xmax)
            ax[0,0].grid(True)
            ax[0,0].invert_xaxis();
            ax[0,0].set_aspect('equal');
#            display.plot_colorbar(ax=ax[0],shrink=0.8);
            display.plot_rhi("VEL", ax=ax[1,0], sweep=s, vmin=-10.66, vmax=10.66, gatefilter=gatefilter,
                             norm=None, cmap=vel_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
            ax[1,0].set_ylim(0,hmax)
            ax[1,0].set_xlim(xmin,xmax)
            ax[1,0].grid(True)
            ax[1,0].invert_xaxis();
            ax[1,0].set_aspect('equal','box');
            display.plot_rhi("WIDTH", ax=ax[1,1], sweep=s, 
                             norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), gatefilter=gatefilter,
                             cmap=spw_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
            ax[1,1].set_ylim(0,hmax)
            ax[1,1].set_xlim(xmin,xmax)
            ax[1,1].grid(True)
            ax[1,1].invert_xaxis();
            ax[1,1].set_aspect('equal','box');
            display.plot_rhi("LDR", ax=ax[0,1], sweep=s, vmin=-35, vmax=5, norm=None,  gatefilter=gatefilter,
                             cmap=ldr_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
            ax[0,1].set_ylim(0,hmax)
            ax[0,1].set_xlim(xmin,xmax)
            ax[0,1].grid(True)
            ax[0,1].invert_xaxis();
            ax[0,1].set_aspect('equal','box');
            sweep_start_index = RadarDS.get_start(s);
            dtime_sweep = cftime.num2pydate(RadarDS.time['data'][sweep_start_index],RadarDS.time['units']);
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
            figname = f'ncas-radar-mobile-ka-band-1_chilbolton_ccrest-m_{dtime_sweep_str}_rhi_az{rhi_az:0.2f}_l1_{product_version}.png';
            plt.savefig(os.path.join(figpath_ccrest_az,figname),dpi=300);
            plt.close();

    else:
        fig, ax = plt.subplots(4,1,figsize=(20,20)); #,constrained_layout=True)
        #fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)
        display.plot_rhi("DBZ", ax=ax[0], sweep=0, vmin=-40, vmax=40, 
                         norm=None, cmap=dbz_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
        ax[0].set_ylim(0,hmax)
        ax[0].set_xlim(xmin,xmax)
        ax[0].grid(True)
        ax[0].invert_xaxis();
        ax[0].set_aspect('equal','box');
        display.plot_rhi("VEL", ax=ax[1], sweep=0, vmin=-10.66, vmax=10.66,  gatefilter=gatefilter,
                         norm=None, cmap=vel_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
        ax[1].set_ylim(0,hmax)
        ax[1].set_xlim(xmin,xmax)
        ax[1].grid(True)
        ax[1].invert_xaxis();
        ax[1].set_aspect('equal','box');
        display.plot_rhi("WIDTH", ax=ax[2], sweep=0, gatefilter=gatefilter,
                         norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), 
                         cmap=spw_cmap, colorbar_orient='horizontal',reverse_xaxis=True)
        ax[2].set_ylim(0,hmax)
        ax[2].set_xlim(xmin,xmax)
        ax[2].grid(True)
        ax[2].invert_xaxis();
        ax[2].set_aspect('equal','box'); 
        display.plot_rhi("LDR", ax=ax[3], sweep=0,vmin=-35, vmax=5, norm=None, gatefilter=gatefilter,
                         cmap=ldr_cmap, colorbar_orient='horizontal')
        ax[3].set_ylim(0,hmax)
        ax[3].set_xlim(xmin,xmax)
        ax[3].grid(True)
        ax[3].invert_xaxis();
        ax[3].set_aspect('equal','box');
        figname = f'ncas-radar-mobile-ka-band-1_chilbolton_ccrest-m_{dtime0_str}_l1_{product_version}.png'
        plt.savefig(os.path.join(figpath,figname),dpi=300);
        plt.close();






def make_ccrest_vpt_plot_day(datestr,inpath,figpath,blflag=False):
    
    if blflag:
        hmax = 4;
    else:
        hmax = 12;

    dbz_cmap = 'pyart_HomeyerRainbow';
    vel_cmap = 'pyart_balance';
    ldr_cmap = 'pyart_SpectralExtended';
    #ldr_cmap = 'pyart_ChaseSpectral';
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
    gatefilter.exclude_below('SNR', -20)
    display = pyart.graph.RadarDisplay(RadarSweepDS);

    import cftime
    dtime = cftime.num2pydate(RadarDS_VPT.time['data'],RadarDS_VPT.time['units'])

    dt_min = dtime[0].replace(hour=0,minute=0,second=0); #datetime.strptime(datestr, '%Y%m%d')
    dt_max = dt_min + datetime.timedelta(days=1)

    time_str = dtime[0].strftime("%Y-%m-%d");
    
   
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "DBZ")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter,
                     vmin=-40, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, colorbar_orient='horizontal')
    ax[0].set_ylim(0,hmax)
    ax[0].grid(True);
    field_name = pyart.graph.common.generate_field_name(RadarSweepDS, "VEL")
    title = f"{pyart.graph.common.generate_radar_name(RadarSweepDS)} {time_str} " + "\n" + field_name;
    display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, title=title, edges=False, gatefilter=gatefilter,
                     vmin=velmin, vmax=velmax, norm=None, cmap=vel_cmap, colorbar_orient='horizontal')
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

    for s in range(1,nsweeps):
        RadarSweepDS = RadarDS_VPT.extract_sweeps([s])
        print(f"sweep = {s}");
        gatefilter = pyart.correct.GateFilter(RadarSweepDS)
        gatefilter.exclude_below('SNR', -20)
        display = pyart.graph.RadarDisplay(RadarSweepDS);
        display.plot_vpt("DBZ", ax=ax[0], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                         vmin=-40, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, 
                         colorbar_flag=False, title_flag=False)
        display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                         vmin=velmin, vmax=velmax, norm=None, cmap=vel_cmap, colorbar_flag=False, title_flag=False)
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
                             vmin=-40, vmax=40, norm=None, filter_transitions=True, cmap=dbz_cmap, 
                             colorbar_flag=False, title_flag=False)
            display.plot_vpt("VEL", ax=ax[1], time_axis_flag=True, edges=False, gatefilter=gatefilter,
                             vmin=velmin, vmax=velmax, norm=None, cmap=vel_cmap, colorbar_flag=False, title_flag=False)
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

    figname = f'ncas-radar-mobile-ka-band-1_chilbolton_ccrest-m_{datestr}_vpt_l1_{product_version}.png';
    
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



def make_ccrest1_rhi_plots_day(datestr,inpath,figpath,blflag=False):
    inpath_date = os.path.join(inpath,datestr);
    
    os.chdir(inpath_date);
    ccrest1_rhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*ccrest1*.nc'.format(datestr))]

    print(f'ccrest-1 files = ',ccrest1_rhi_files);

    for f in ccrest1_rhi_files:
        make_ccrest_rhi_plot(f,figpath,'ccrest-1',blflag=blflag);

    return

def make_ccrest2_rhi_plots_day(datestr,inpath,figpath,blflag=False):
    inpath_date = os.path.join(inpath,datestr);
    
    os.chdir(inpath_date);
    ccrest2_rhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*ccrest2*.nc'.format(datestr))]

    print(f'ccrest-2 files = ',ccrest2_rhi_files);

    for f in ccrest2_rhi_files:
        make_ccrest_rhi_plot(f,figpath,'ccrest-2',blflag=blflag);

    return


make_ccrest_vpt_plot_day(datestr,inpath,figpath,blflag=False);
make_ccrest1_rhi_plots_day(datestr,inpath,figpath,blflag=False);
make_ccrest2_rhi_plots_day(datestr,inpath,figpath,blflag=False);
zlevels = np.arange(100, 15000, 100);  # height above radar
make_ccrest_vad_plot_day(datestr,inpath,figpath,zlevels);

