#!/usr/bin/env python
# coding: utf-8

import netCDF4 as nc4
import os

import pyart
import datetime
import numpy as np
import numpy.ma as ma
import shutil
import glob
import gzip

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cmocean
import getpass, socket

import pandas as pd

import cftime

version = 0.1

import sys
sys.path.append('/home/users/cjwalden/my-packages')

import kepler_utils

datestr = '20190509';

mmclxpath = os.path.join('/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data/campaign/danebury-hill/raw/',datestr);
#os.chdir('/gws/nopw/j04/ncas_obs');
os.chdir(mmclxpath);
files = [os.path.join(mmclxpath,f) for f in glob.glob('*.mmclx')]



scan_types = ['ppi','rhi','vert','man']

for f in files:
    for elem in scan_types:
        scan_type = None
        if elem in f.lower():
            scan_type = elem
            break;

    print("{} {}".format(f, scan_type));

import kepler_utils as kepler

dir(kepler_utils)

print(files[1]);
Radar = kepler.read_mira35_mmclx(files[5]);

outpath = '/home/users/cjwalden/dataproc/woest'

outfile = os.path.join(outpath,'ncas-mobile-radar-ka-band-1_{}_{}.nc'.format(datestr,scan_type));

print(Radar.sweep_mode['data'])

pyart.io.write_cfradial(outfile, Radar, format='NETCDF4', time_reference=None);


# In[14]:


import geopandas as gpd
import fiona
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
#gpd.io.file.fiona.drvsupport.supported_drivers['kml'] = 'rw#'
#gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
#gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
#fiona.drvsupport.supported_drivers['kml'] = 'rw'  # enable KML support which is disabled by default
##fiona.drvsupport.supported_drivers['KML'] = 'rw'  # enable KML support which is disabled by default
#fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'  # enable KML support which is disabled by default
WesconBox1 = gpd.read_file('/home/users/cjwalden/WesConBox1.kml', driver='LIBKML')
WardonHill = gpd.read_file('/home/users/cjwalden/WardonHill.kml', driver='LIBKML')


filepath = '/home/users/cjwalden/WesConGrid_Verticals.kml';
gdf_list = []
for layer in fiona.listlayers(filepath):    
    gdf = gpd.read_file(filepath, driver='LIBKML', layer=layer)
    gdf_list.append(gdf)

WesConGrid_Verticals = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))

filepath = '/home/users/cjwalden/WesConGrid_Horizontals.kml';
gdf_list = []
for layer in fiona.listlayers(filepath):    
    gdf = gpd.read_file(filepath, driver='LIBKML', layer=layer)
    gdf_list.append(gdf)

WesConGrid_Horizontals = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True))


surface_sites = gpd.read_file('/home/users/cjwalden/surface_sites.kml', driver='LIBKML')
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import OSM
from cartopy.feature import NaturalEarthFeature, LAND, COASTLINE
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt

osm_tiles = OSM()


display = pyart.graph.RadarMapDisplay(Radar);
# Setting projection and ploting the second tilt
projection = osm_tiles.crs
#ccrs.LambertConformal(
#    central_latitude=Radar.latitude["data"],
#    central_longitude=Radar.longitude["data"],
#)


fig, ax  = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=projection))
display.plot_ppi_map(
    "Zg",
    0,
    vmin=-40,
    vmax=40,
    min_lon=-3.6,
    max_lon=-1.0,
    min_lat=50.4,
    max_lat=51.6,
    #lon_lines=np.arange(-3.,-1.,0.2),
    resolution="10m",
    #lat_lines=np.arange(50.5,51.5, 0.1),
    projection=projection,
    fig=fig,
    ax = ax,
    lat_0=Radar.latitude["data"],
    lon_0=Radar.longitude["data"],
    cmap='pyart_HomeyerRainbow',
    colorbar_flag=False
)


ax.set_extent([-3.6, -1.0, 50.4, 51.6],ccrs.PlateCarree())
ax.add_image(osm_tiles, 9)

#ax.add_feature(LAND)
#ax.add_feature(COASTLINE)
gl = ax.gridlines(draw_labels=True)
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.bottom_labels = True
gl.top_labels = False
gl.right_labels = False

# Plot range rings at 10, 20, 30, 40km
display.plot_range_ring(10.0, line_style="k-",lw=0.75)
display.plot_range_ring(20.0, line_style="k--",lw=0.75)
display.plot_range_ring(30.0, line_style="k-",lw=0.75)
display.plot_range_ring(40.0, line_style="k--",lw=0.75)

# Plot cross hairs
display.plot_line_xy(
    np.array([-40000.0, 40000.0]), np.array([0.0, 0.0]), line_style="k-",
lw=0.75)
display.plot_line_xy(
    np.array([0.0, 0.0]), np.array([-40000.0, 40000.0]), line_style="k-",
lw=0.75)

# Indicate the radar location with a point
display.plot_point(Radar.longitude["data"], Radar.latitude["data"])
display.plot_colorbar(orient='horizontal',shrink=0.8);
#ax.coastlines('10m')

wesconbox1_ae = WesconBox1.to_crs(projection)
wardonhill_ae = WardonHill.to_crs(projection)
surfacesites_ae = surface_sites.to_crs(projection)
wescongrid_v_ae = WesConGrid_Verticals.to_crs(projection)
wescongrid_h_ae = WesConGrid_Horizontals.to_crs(projection)

# Here's what the plot looks like in GeoPandas
wesconbox1_ae.plot(ax=ax,color='None',edgecolor='Blue')
wardonhill_ae.plot(ax=ax,color='Magenta')
surfacesites_ae.plot(ax=ax,color='Magenta')
wescongrid_v_ae.plot(ax=ax);
wescongrid_h_ae.plot(ax=ax);

print(wescongrid_v_ae.loc[0]);
from shapely.geometry import Point
first = Point(wescongrid_v_ae.loc[0].geometry.coords)

print(first)
xgrid=first.x;
ygrid=first.y;

print(xgrid,ygrid);
#ax.annotate('A',xy=(xgrid,ygrid), xytext=(3,3),textcoords='offset_points')

xwh = wardonhill_ae.loc[0].geometry.x
ywh = wardonhill_ae.loc[0].geometry.y
ax.annotate('Wardon Hill',xy=(xwh,ywh),xytext=(3,3),textcoords='offset points');

plt.show()


# In[16]:


DS = nc4.Dataset(outfile);

print(DS);
newradar = pyart.io.read_cfradial(outfile);

display = pyart.graph.RadarMapDisplay(newradar);
# Setting projection and ploting the second tilt
projection = osm_tiles.crs

fig, ax  = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=projection))
display.plot_ppi_map(
    "Zg",
    0,
    vmin=-40,
    vmax=40,
    min_lon=-3.6,
    max_lon=-1.0,
    min_lat=50.4,
    max_lat=51.6,
    #lon_lines=np.arange(-3.,-1.,0.2),
    resolution="10m",
    #lat_lines=np.arange(50.5,51.5, 0.1),
    projection=projection,
    fig=fig,
    ax = ax,
    lat_0=newradar.latitude["data"],
    lon_0=newradar.longitude["data"],
    cmap='pyart_HomeyerRainbow',
    colorbar_flag=False
)






