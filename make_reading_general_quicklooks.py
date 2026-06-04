#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
make_reading_general_quicklooks.py

Generate quicklook plots for Reading General campaign radar data.

Produces daily VPT (Vertical Pointing) time-height plots of DBZ, VEL, WIDTH
and LDR from the Reading Mobile Ka-band Radar (Kepler) at Chilbolton.

Usage:
    python make_reading_general_quicklooks.py -d YYYYMMDD
    python make_reading_general_quicklooks.py -d YYYYMMDD -i /path/to/input -o /path/to/output

Arguments:
    -d, --date:     Date string in YYYYMMDD format (default: today)
    -i, --inpath:   Input directory containing processed CF-Radial files
                    (default: /data/processing/kepler/reading-general)
    -o, --outpath:  Output directory for quicklook images
                    (default: <inpath>/quicklooks)
    -b:             Boundary layer mode (limit height axis to 4 km)

Author: Chris Walden, Science and Technology Facilities Council (STFC)
        as part of UK Research and Innovation (UKRI)
"""

import getopt
import re
import sys
import os
import datetime
import glob

import netCDF4 as nc4
import pyart
import numpy as np
import numpy.ma as ma
import cftime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib import colors

from kepler_utils import get_valid_sweep_indices

# ── Configuration ─────────────────────────────────────────────────────────────

CAMPAIGN = 'reading-general'
INSTRUMENT_NAME = 'reading-mobile-ka-band-radar-1'
RADAR_LABEL = 'University of Reading Ka-band Cloud Radar (Kepler)'
RADAR_LOCATION = 'UKRI-STFC Chilbolton Observatory'
DEFAULT_INPATH = '/data/processing/kepler/reading-general'

COLORMAPS = {
    'dbz': 'HomeyerRainbow',
    'vel': 'balance',
    'width': 'SpectralExtended',
    'ldr': 'SpectralExtended',
}


# ── Location helpers ─────────────────────────────────────────────────────────

def _sanitize_location(value):
    """Normalize location strings for filename-safe tokens."""
    text = str(value).strip().lower()
    text = text.replace(' ', '-')
    text = re.sub(r'[^a-z0-9_-]+', '', text)
    return text if text else 'unknown'


def get_location_from_dataset(ds):
    """Extract location from NetCDF global attributes with fallbacks."""
    for attr in ('location', 'platform_location', 'site_name'):
        if hasattr(ds, attr):
            return _sanitize_location(getattr(ds, attr))
    if hasattr(ds, 'platform'):
        if 'cao' in str(getattr(ds, 'platform')).lower():
            return 'cao'
    return 'unknown'


# ── Geospatial helpers ────────────────────────────────────────────────────────

def _parse_geospatial_bounds(bounds_str):
    """Parse a geospatial_bounds string and return (lon_min, lon_max, lat_min, lat_max) or None."""
    m = re.search(
        r'([+-]?\d+\.?\d*)([NS])\s+([+-]?\d+\.?\d*)([EW])'
        r'\s*,\s*'
        r'([+-]?\d+\.?\d*)([NS])\s+([+-]?\d+\.?\d*)([EW])',
        bounds_str,
    )
    if m is None:
        return None
    lat1 = float(m.group(1)) * (-1 if m.group(2) == 'S' else 1)
    lon1 = float(m.group(3)) * (-1 if m.group(4) == 'W' else 1)
    lat2 = float(m.group(5)) * (-1 if m.group(6) == 'S' else 1)
    lon2 = float(m.group(7)) * (-1 if m.group(8) == 'W' else 1)
    return min(lon1, lon2), max(lon1, lon2), min(lat1, lat2), max(lat1, lat2)


def _compute_day_outermost_ppi_bounds(ppi_files):
    """Return the outermost lon/lat envelope across all parseable PPI files for a day."""
    day_bounds = None
    print(f'Computing day-outermost bounds from {len(ppi_files)} PPI files')
    for ncfile in ppi_files:
        try:
            with nc4.Dataset(ncfile) as ds:
                bounds = _parse_geospatial_bounds(getattr(ds, 'geospatial_bounds', ''))
                if bounds is None:
                    lo_min = getattr(ds, 'geospatial_lon_min', None)
                    lo_max = getattr(ds, 'geospatial_lon_max', None)
                    la_min = getattr(ds, 'geospatial_lat_min', None)
                    la_max = getattr(ds, 'geospatial_lat_max', None)
                    if None not in (lo_min, lo_max, la_min, la_max):
                        bounds = (float(lo_min), float(lo_max), float(la_min), float(la_max))
        except Exception as exc:
            print(f'Warning: cannot read bounds from {os.path.basename(ncfile)}: {exc}')
            continue
        if bounds is None:
            continue
        if day_bounds is None:
            day_bounds = list(bounds)
        else:
            day_bounds[0] = min(day_bounds[0], bounds[0])
            day_bounds[1] = max(day_bounds[1], bounds[1])
            day_bounds[2] = min(day_bounds[2], bounds[2])
            day_bounds[3] = max(day_bounds[3], bounds[3])
    return tuple(day_bounds) if day_bounds is not None else None


# ── PPI file collection ───────────────────────────────────────────────────────

def collect_ppi_files(datestr, inpath):
    """Return sorted list of PPI CF-Radial files for a given day."""
    inpath_date = os.path.join(inpath, datestr)
    return sorted(glob.glob(os.path.join(inpath_date, '*ppi*.nc')))


# ── Command-line parsing ───────────────────────────────────────────────────────

def parse_command_line():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], 'd:i:o:b',
            ['date=', 'inpath=', 'outpath=', 'vpt-only', 'ppi-only', 'vad-only',
             'rhi-only', 'no-rhi', 'no-vad', 'ppi-map-day-extent',
             'lat-min=', 'lat-max=', 'lon-min=', 'lon-max=', 'az-offset=', 'alpha=']
        )
    except getopt.GetoptError as err:
        print(f'Error: {err}')
        print('Usage: python make_reading_general_quicklooks.py -d YYYYMMDD '
              '[-i inpath] [-o outpath] [-b] [--vpt-only] [--ppi-only] [--vad-only] '
              '[--rhi-only] [--no-rhi] [--no-vad] [--ppi-map-day-extent] '
              '[--lat-min LAT] [--lat-max LAT] [--lon-min LON] [--lon-max LON] '
              '[--az-offset DEGREES] [--alpha 0-1]')
        sys.exit(2)

    datestr = datetime.datetime.now().strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False
    vpt_only = False
    ppi_only = False
    vad_only = False
    rhi_only = False
    no_rhi = False
    no_vad = False
    ppi_map_day_extent = False
    lat_min = None
    lat_max = None
    lon_min = None
    lon_max = None
    az_offset = 0.0
    map_alpha = 0.25

    for option, argument in opts:
        if option in ('-d', '--date'):
            datestr = argument
        elif option in ('-i', '--inpath'):
            inpath = argument
        elif option in ('-o', '--outpath'):
            outpath = argument
        elif option == '-b':
            blflag = True
        elif option == '--vpt-only':
            vpt_only = True
        elif option == '--ppi-only':
            ppi_only = True
        elif option == '--vad-only':
            vad_only = True
        elif option == '--rhi-only':
            rhi_only = True
        elif option == '--no-rhi':
            no_rhi = True
        elif option == '--no-vad':
            no_vad = True
        elif option == '--ppi-map-day-extent':
            ppi_map_day_extent = True
        elif option == '--lat-min':
            lat_min = float(argument)
        elif option == '--lat-max':
            lat_max = float(argument)
        elif option == '--lon-min':
            lon_min = float(argument)
        elif option == '--lon-max':
            lon_max = float(argument)
        elif option == '--az-offset':
            az_offset = float(argument)
        elif option == '--alpha':
            map_alpha = float(argument)
        else:
            print(f'Unhandled option: {option}')
            sys.exit(2)

    # Build explicit extent tuple if any bound was supplied
    explicit_extent = None
    if any(v is not None for v in (lon_min, lon_max, lat_min, lat_max)):
        if not all(v is not None for v in (lon_min, lon_max, lat_min, lat_max)):
            print('Error: all four of --lat-min, --lat-max, --lon-min, --lon-max must be provided together')
            sys.exit(2)
        explicit_extent = (lon_min, lon_max, lat_min, lat_max)

    return datestr, inpath, outpath, blflag, vpt_only, ppi_only, vad_only, rhi_only, no_rhi, no_vad, ppi_map_day_extent, explicit_extent, az_offset, map_alpha


# ── Path setup ─────────────────────────────────────────────────────────────────

def setup_paths(datestr, inpath=None, outpath=None):
    if inpath is None:
        inpath = DEFAULT_INPATH

    if outpath is None:
        figpath = os.path.join(inpath, 'quicklooks')
    else:
        figpath = outpath

    os.makedirs(figpath, exist_ok=True)
    return inpath, figpath


# ── VPT quicklook ─────────────────────────────────────────────────────────────

def _plot_vpt_fields(display, axes, radar_ds, gatefilter, time_str, hmax,
                     vel_min, vel_max):
    """Plot VPT fields with per-panel titles and colorbars (NCAS style)."""
    vpt_specs = [
        ('DBZ',   axes[0], dict(vmin=-40, vmax=40,   cmap=COLORMAPS['dbz'])),
        ('VEL',   axes[1], dict(vmin=vel_min, vmax=vel_max, cmap=COLORMAPS['vel'])),
        ('WIDTH', axes[2], dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1),
                                                     vmax=np.sqrt(1e1)),
                                cmap=COLORMAPS['width'])),
        ('LDR',   axes[3], dict(vmin=-35, vmax=5,    cmap=COLORMAPS['ldr'])),
    ]
    for field, ax, kwargs in vpt_specs:
        field_name = pyart.graph.common.generate_field_name(radar_ds, field)
        display.plot_vpt(field, ax=ax, time_axis_flag=True, edges=False,
                         gatefilter=gatefilter, colorbar_flag=False,
                         title=field_name, title_flag=True,
                         filter_transitions=False, **kwargs)
        ax.set_ylim(0, hmax)
        mappable = display.plots[-1]
        fld_meta = display.fields[field]
        label = (f"{fld_meta.get('long_name', field)} "
                 f"({fld_meta.get('units', '')})")
        cb = ax.get_figure().colorbar(mappable, ax=ax,
                                      orientation='horizontal',
                                      shrink=0.9, pad=0.20)
        cb.set_label(label)


def _plot_vpt_fields_overlay(display, axes, gatefilter, vel_min, vel_max):
    """Overlay additional VPT sweeps without redrawing titles or colorbars."""
    display.plot_vpt('DBZ',   ax=axes[0], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=-40, vmax=40,
                     cmap=COLORMAPS['dbz'],
                     colorbar_flag=False, title_flag=False,
                     filter_transitions=False)
    display.plot_vpt('VEL',   ax=axes[1], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=vel_min, vmax=vel_max,
                     cmap=COLORMAPS['vel'],
                     colorbar_flag=False, title_flag=False,
                     filter_transitions=False)
    display.plot_vpt('WIDTH', ax=axes[2], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter,
                     norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1),
                                         vmax=np.sqrt(1e1)),
                     cmap=COLORMAPS['width'],
                     colorbar_flag=False, title_flag=False,
                     filter_transitions=False)
    display.plot_vpt('LDR',   ax=axes[3], time_axis_flag=True, edges=False,
                     gatefilter=gatefilter, vmin=-35, vmax=5,
                     cmap=COLORMAPS['ldr'],
                     colorbar_flag=False, title_flag=False,
                     filter_transitions=False)


def make_vpt_plot_day(datestr, inpath, figpath, blflag=False):
    """Generate a daily VPT time-height plot (DBZ, VEL, WIDTH, LDR)."""

    hmax = 4 if blflag else 12

    inpath_date = os.path.join(inpath, datestr)
    vpt_files = sorted(glob.glob(os.path.join(inpath_date, f'*{datestr}*vpt*.nc')))

    if not vpt_files:
        print(f'No VPT files found for {datestr} in {inpath_date}')
        return

    print(f'Found {len(vpt_files)} VPT file(s) for {datestr}')

    # Read product version and velocity limits from the first file
    with nc4.Dataset(vpt_files[0]) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
    radar_first = pyart.io.read_cfradial(vpt_files[0])
    vel_field = radar_first.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -5.0)
    vel_max = vel_field.get('field_limit_upper',  5.0)

    # Time axis limits for the day
    dtime = cftime.num2pydate(radar_first.time['data'], radar_first.time['units'])
    dt_min = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dtime[0].strftime('%Y-%m-%d')

    fig, axes = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=False)
    fig.suptitle(f'{RADAR_LABEL} | {RADAR_LOCATION}\n{time_str}', fontsize=13, fontweight='bold', y=0.99)
    fig.subplots_adjust(left=0.10, right=0.95, top=0.94, bottom=0.10,
                        hspace=0.50)

    colorbars_done = False

    # Try previous day for midnight overspill only when today's VPT starts
    # after 00:00:00.  If it starts at midnight the spillover has already been
    # prepended/merged into the file by fix_midnight_crossing_vpt, so there is
    # nothing extra to overlay.
    if dtime[0] > dt_min:
        prev_date = (datetime.datetime.strptime(datestr, '%Y%m%d')
                     - datetime.timedelta(days=1))
        prevstr = prev_date.strftime('%Y%m%d')
        prev_files = sorted(glob.glob(
            os.path.join(inpath, prevstr, f'*{prevstr}*vpt*.nc')))

        prev_sweeps = []  # list of (radar_obj, sweep_index) in chronological order
        for pf in prev_files:
            try:
                rp = pyart.io.read_cfradial(pf)
            except Exception:
                continue
            ray_times = cftime.num2pydate(rp.time['data'], rp.time['units'])
            # Only bother with this file if any ray falls in the current day
            if not any(t >= dt_min for t in ray_times):
                continue
            for s in range(rp.nsweeps):
                s0, s1 = rp.get_start(s), rp.get_end(s)
                if any(t >= dt_min for t in ray_times[s0:s1 + 1]):
                    prev_sweeps.append((rp, s))

        if prev_sweeps:
            print(f'  Including {len(prev_sweeps)} midnight-overspill sweep(s) '
                  f'from {prevstr}')
        for rp, sidx in prev_sweeps:
            sweep = rp.extract_sweeps([sidx])
            if np.all(sweep.antenna_transition['data'] == 1):
                continue
            gf = pyart.correct.GateFilter(sweep)
            gf.exclude_below('SNR', -20)
            gf.exclude_transition()
            disp = pyart.graph.RadarDisplay(sweep)
            _plot_vpt_fields_overlay(disp, axes, gf, vel_min, vel_max)

    # All sweeps for the current day across all VPT files
    for vpt_file in vpt_files:
        try:
            radar = pyart.io.read_cfradial(vpt_file)
        except Exception as exc:
            print(f'Error reading {os.path.basename(vpt_file)}: {exc}')
            continue
        for s in range(radar.nsweeps):
            sweep = radar.extract_sweeps([s])
            if np.all(sweep.antenna_transition['data'] == 1):
                continue
            gf = pyart.correct.GateFilter(sweep)
            gf.exclude_below('SNR', -20)
            gf.exclude_transition()
            disp = pyart.graph.RadarDisplay(sweep)
            if not colorbars_done:
                _plot_vpt_fields(disp, axes, sweep, gf, time_str, hmax,
                                 vel_min, vel_max)
                colorbars_done = True
            else:
                _plot_vpt_fields_overlay(disp, axes, gf, vel_min, vel_max)

    for ax in axes:
        ax.set_xlim(dt_min, dt_max)
        ax.set_xlabel('Time (UTC)')
        ax.grid(True)

    figpath_vpt = os.path.join(figpath, 'vpt')
    os.makedirs(figpath_vpt, exist_ok=True)

    figname = (f'{INSTRUMENT_NAME}_cao_{CAMPAIGN}_{datestr}_vpt_l1'
               f'_{product_version}.png')
    if blflag:
        figname = figname.replace('.png', '_bl.png')

    plt.savefig(os.path.join(figpath_vpt, figname), dpi=300)
    print(f'Saved: {figname}')
    plt.close()


# ── PPI quicklooks ────────────────────────────────────────────────────────────

def make_ppi_plot(ncfile, datestr, figpath, blflag=False, az_offset=0.0):
    """Create a standard (Cartesian x/y) PPI quicklook for each sweep in ncfile."""
    max_range = 20.0 if blflag else 30.0

    with nc4.Dataset(ncfile) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
        location = get_location_from_dataset(ds)

    radar = pyart.io.read_cfradial(ncfile)
    if az_offset:
        radar.azimuth['data'] = (radar.azimuth['data'] + az_offset) % 360.0
    display = pyart.graph.RadarDisplay(radar)

    vel_field = radar.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10.66)
    vel_max = vel_field.get('field_limit_upper', 10.66)

    outdir = os.path.join(figpath, 'ppi', datestr)
    os.makedirs(outdir, exist_ok=True)

    valid_sweeps = get_valid_sweep_indices(radar)

    for s in valid_sweeps:
        sweep_start = radar.get_start(s)
        dtime_sweep = cftime.num2pydate(radar.time['data'][sweep_start], radar.time['units'])
        dtime_str = dtime_sweep.strftime('%Y%m%d-%H%M%S')
        ppi_el = float(np.nanmean(radar.get_elevation(s)))

        gf = pyart.correct.GateFilter(radar)
        gf.exclude_below('SNR', -20)

        fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=False)
        fig.suptitle(
            f"{CAMPAIGN.upper()} | {location.upper()}\n"
            f"{INSTRUMENT_NAME} PPI El {ppi_el:.1f}\u00b0 "
            f"{dtime_sweep.strftime('%Y-%m-%dT%H:%M:%SZ')}",
            fontsize=14, y=0.98,
        )
        fig.subplots_adjust(left=0.07, right=0.97, top=0.87, bottom=0.08,
                            wspace=0.28, hspace=0.62)

        field_specs = [
            ('DBZ', 'Equivalent Reflectivity Factor', axes[0, 0],
             dict(vmin=-40, vmax=40, cmap='HomeyerRainbow', mask_outside=True)),
            ('LDR', 'Radar Linear Depolarization Ratio', axes[0, 1],
             dict(vmin=-35, vmax=5, cmap='SpectralExtended', mask_outside=True)),
            ('VEL', 'Radial Velocity', axes[1, 0],
             dict(vmin=vel_min, vmax=vel_max, cmap='balance', mask_outside=True)),
            ('WIDTH', 'Radar Doppler Spectrum Width', axes[1, 1],
             dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
                  cmap='SpectralExtended', mask_outside=False)),
        ]

        for field, title, ax, kwargs in field_specs:
            display.plot_ppi(field, s, ax=ax, colorbar_flag=True, title_flag=False,
                             gatefilter=gf, filter_transitions=True,
                             **kwargs)
            display.plot_range_rings([10, 20, int(max_range)], ax=ax, col='0.5', lw=0.8)
            ax.set_xlim(-max_range, max_range)
            ax.set_ylim(-max_range, max_range)
            ax.set_aspect('equal')
            ax.grid(True)
            ax.set_title(title, fontsize=10, pad=6)

        figname = (
            f'{INSTRUMENT_NAME}_{location}_{CAMPAIGN}_{dtime_str}_'
            f'ppi_el{ppi_el:0.2f}_l1_{product_version}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=300)
        plt.close()
        print(f'Saved PPI plot: {figname}')


def make_ppi_map_plot(ncfile, datestr, figpath, blflag=False,
                      day_extent=None, colorbar_shrink=0.86, az_offset=0.0,
                      map_alpha=0.25):
    """Create a georeferenced (map) PPI quicklook for each sweep in ncfile."""
    try:
        import cartopy.crs as ccrs
    except ImportError:
        print('Skipping PPI map quicklook: cartopy is not available')
        return

    try:
        import contextily as ctx
    except ImportError:
        ctx = None

    max_range = 20.0 if blflag else 30.0

    with nc4.Dataset(ncfile) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
        location = get_location_from_dataset(ds)
        _bounds = _parse_geospatial_bounds(getattr(ds, 'geospatial_bounds', ''))

    radar = pyart.io.read_cfradial(ncfile)
    if az_offset:
        radar.azimuth['data'] = (radar.azimuth['data'] + az_offset) % 360.0
    display = pyart.graph.RadarMapDisplay(radar)

    vel_field = radar.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10.66)
    vel_max = vel_field.get('field_limit_upper', 10.66)

    outdir = os.path.join(figpath, 'ppi_map', datestr)
    os.makedirs(outdir, exist_ok=True)

    valid_sweeps = get_valid_sweep_indices(radar)

    field_specs = [
        ('DBZ', 'Equivalent Reflectivity Factor',
         dict(vmin=-40, vmax=40, cmap='HomeyerRainbow')),
        ('LDR', 'Radar Linear Depolarization Ratio',
         dict(vmin=-35, vmax=5, cmap='SpectralExtended')),
        ('VEL', 'Radial Velocity',
         dict(vmin=vel_min, vmax=vel_max, cmap='balance')),
        ('WIDTH', 'Radar Doppler Spectrum Width',
         dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
              cmap='SpectralExtended')),
    ]

    for s in valid_sweeps:
        sweep_start = radar.get_start(s)
        dtime_sweep = cftime.num2pydate(radar.time['data'][sweep_start], radar.time['units'])
        dtime_str = dtime_sweep.strftime('%Y%m%d-%H%M%S')
        ppi_el = float(np.nanmean(radar.get_elevation(s)))

        radar_lon = radar.longitude['data'][0]
        radar_lat = radar.latitude['data'][0]
        # AzimuthalEquidistant: same projection approach as make_picasso_quicklooks.py
        # where alpha=map_alpha on plot_ppi_map is known to work correctly.
        projection = ccrs.AzimuthalEquidistant(central_longitude=radar_lon,
                                               central_latitude=radar_lat)
        if day_extent is not None:
            lon_min, lon_max, lat_min, lat_max = day_extent
        elif _bounds is not None:
            lon_min, lon_max, lat_min, lat_max = _bounds
        else:
            lon_min, lon_max = radar_lon - 0.50, radar_lon + 0.50
            lat_min, lat_max = radar_lat - 0.35, radar_lat + 0.35

        fig, axes = plt.subplots(
            2, 2, figsize=(16, 16),
            subplot_kw=dict(projection=projection),
            constrained_layout=False,
        )
        fig.subplots_adjust(left=0.07, right=0.97, top=0.90,
                            bottom=0.05, wspace=0.25, hspace=0.25)
        fig.suptitle(
            f"{CAMPAIGN.upper()} | {location.upper()}\n"
            f"{INSTRUMENT_NAME} PPI Map | El {ppi_el:.1f}\u00b0 "
            f"{dtime_sweep.strftime('%Y-%m-%dT%H:%M:%SZ')}",
            fontsize=12, y=0.97,
        )

        for (field, short_title, kwargs), panel_ax in zip(field_specs, axes.flatten()):
            # Draw radar with alpha — basemap tiles default to zorder=0 so
            # pcolormesh (zorder=1) naturally sits on top when basemap is added after.
            display.plot_ppi_map(
                field, s,
                projection=projection, ax=panel_ax,
                lon_lines=None, lat_lines=None,
                colorbar_flag=False, title_flag=False,
                filter_transitions=True, gatefilter=None,
                alpha=map_alpha,
                **kwargs,
            )
            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max],
                                 crs=ccrs.PlateCarree())

        for (field, short_title, _kwargs), panel_ax in zip(field_specs, axes.flatten()):
            # Add basemap after all radar panels are drawn (same ordering as Picasso).
            # Pass the same CRS object used for plot_ppi_map so contextily reprojects correctly.
            if ctx is not None:
                try:
                    ctx.add_basemap(panel_ax, zoom=13,
                                    source=ctx.providers.OpenTopoMap,
                                    crs=projection, reset_extent=False)
                except Exception:
                    try:
                        ctx.add_basemap(panel_ax, zoom=13,
                                        source=ctx.providers.OpenStreetMap.Mapnik,
                                        crs=projection, reset_extent=False)
                    except Exception as e:
                        print(f'Warning: could not add basemap tiles: {e}')
            else:
                import cartopy.feature as cfeature
                panel_ax.add_feature(cfeature.LAND.with_scale('10m'),
                                     facecolor='#f5f5e6', zorder=0)
                panel_ax.add_feature(cfeature.OCEAN.with_scale('10m'),
                                     facecolor='#c8d8e8', zorder=0)
                panel_ax.add_feature(cfeature.COASTLINE.with_scale('10m'),
                                     linewidth=0.6, zorder=1)
                panel_ax.add_feature(cfeature.BORDERS.with_scale('10m'),
                                     linewidth=0.4, linestyle=':', zorder=1)
                panel_ax.add_feature(cfeature.RIVERS.with_scale('10m'),
                                     edgecolor='#7ab6d4', linewidth=0.4, zorder=1)
                panel_ax.add_feature(cfeature.LAKES.with_scale('10m'),
                                     facecolor='#c8d8e8', linewidth=0.4, zorder=1)

            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max],
                                 crs=ccrs.PlateCarree())

        # Collect the radar mappables before adding any overlays
        # display.plots is populated in order: [DBZ_mesh, LDR_mesh, VEL_mesh, WIDTH_mesh]
        radar_mappables = list(display.plots[-len(field_specs):])

        for (field, short_title, _kwargs), panel_ax, mappable in zip(
                field_specs, axes.flatten(), radar_mappables):
            panel_ax.set_title(short_title, fontsize=10, pad=6)
            display.plot_range_rings([10, 20, int(max_range)], ax=panel_ax,
                                     col='0.5', lw=0.8)
            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max],
                                 crs=ccrs.PlateCarree())

            gl = panel_ax.gridlines(draw_labels=True, linewidth=0.4,
                                    color='gray', alpha=0.6,
                                    x_inline=False, y_inline=False)
            gl.right_labels = False
            gl.top_labels = False

            field_meta = display.fields[field]
            label = f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})"
            # Build a full-opacity ScalarMappable so the colorbar is not
            # washed out by the alpha applied to the radar data on the map.
            from matplotlib.cm import ScalarMappable
            opaque_sm = ScalarMappable(norm=mappable.norm,
                                       cmap=mappable.cmap)
            opaque_sm.set_array([])
            cb = fig.colorbar(opaque_sm, ax=panel_ax, orientation='vertical',
                              shrink=colorbar_shrink, pad=0.02)
            cb.set_label(label)
            cb.ax.tick_params(labelsize=8)

        figname = (
            f'{INSTRUMENT_NAME}_{location}_{CAMPAIGN}_{dtime_str}_'
            f'ppi-map_el{ppi_el:0.2f}_l1_{product_version}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=300)
        plt.close()
        print(f'Saved PPI map plot: {figname}')




# ── VAD quicklook ─────────────────────────────────────────────────────────────

def make_vad_plot_day(datestr, inpath, figpath, blflag=False):
    """Generate a daily VAD (Velocity Azimuth Display) wind plot.

    Produces a two-panel time–height figure showing wind direction (from) and
    wind speed (with barbs) retrieved via the Browning & Wexler VAD method
    applied to each sweep in the VAD CF-Radial file for the day.

    Args:
        datestr:  Date string in YYYYMMDD format.
        inpath:   Base input directory (date sub-directory appended automatically).
        figpath:  Base output quicklooks directory (``vad/`` appended automatically).
        blflag:   If True, limit height axis to 4 km.
    """
    hmax = 4 if blflag else 12

    inpath_date = os.path.join(inpath, datestr)
    vad_files = sorted(glob.glob(os.path.join(inpath_date, f'*{datestr}*vad*.nc')))

    if not vad_files:
        print(f'No VAD files found for {datestr} in {inpath_date}')
        return

    print(f'Found {len(vad_files)} VAD file(s) for {datestr}')
    vad_file = vad_files[0]

    with nc4.Dataset(vad_file) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
        location = get_location_from_dataset(ds)

    radar = pyart.io.read_cfradial(vad_file)

    # Build per-sweep time and VAD arrays
    vad_ray_index_start = [radar.sweep_start_ray_index['data'][s] for s in range(radar.nsweeps)]
    vad_ray_index_end   = [radar.sweep_end_ray_index['data'][s]   for s in range(radar.nsweeps)]

    dt_vad_start = cftime.num2pydate(radar.time['data'][vad_ray_index_start], radar.time['units'])
    dt_vad_end   = cftime.num2pydate(radar.time['data'][vad_ray_index_end],   radar.time['units'])

    # Height levels for VAD retrieval (50 m steps from 100 m to hmax km)
    zlevels = np.arange(100, hmax * 1000 + 50, 50, dtype=float)

    u_allsweeps = []
    v_allsweeps = []

    for s in range(radar.nsweeps):
        try:
            radar_1sweep = radar.extract_sweeps([s])
            gf = pyart.correct.GateFilter(radar_1sweep)
            gf.exclude_below('SNR', 0)
            vad = pyart.retrieve.vad_browning(radar_1sweep, 'VEL', z_want=zlevels, gatefilter=gf)
            u_allsweeps.append(vad.u_wind)
            v_allsweeps.append(vad.v_wind)
        except Exception as e:
            print(f'  VAD retrieval failed for sweep {s}: {e}')
            u_allsweeps.append(np.full(len(zlevels), np.nan))
            v_allsweeps.append(np.full(len(zlevels), np.nan))

    u_vel = np.array(u_allsweeps)   # (nsweeps, nlevels)
    v_vel = np.array(v_allsweeps)

    speed       = np.sqrt(u_vel**2 + v_vel**2)
    orientation = np.rad2deg(np.arctan2(-u_vel, -v_vel)) % 360

    speed       = ma.masked_where(speed > 100., speed)
    orientation = ma.masked_where(speed > 100., orientation)

    # Fixed display duration per sweep
    vad_duration = dt_vad_end - dt_vad_start
    vad_duration[:] = datetime.timedelta(minutes=12)
    dt_vad_mid = dt_vad_start + np.array([0.5 * d for d in vad_duration])

    # Day time axis limits
    dt_min = dt_vad_start[0].replace(hour=0, minute=0, second=0, microsecond=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dt_min.strftime('%Y-%m-%d')

    myFmt = mdates.DateFormatter('%H:%M')

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    fig.suptitle(f'{RADAR_LABEL} | {RADAR_LOCATION}\n{time_str}', fontsize=12, fontweight='bold')
    fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2, wspace=0.2)

    z_km = zlevels / 1000.

    # ── Panel 0: wind direction ──────────────────────────────────────────────
    ax0 = axes[0]
    ax0.xaxis.set_major_formatter(myFmt)
    ax0.set_facecolor('white')
    ax0.grid(True)

    h0 = ax0.pcolormesh(
        [dt_vad_start[0] - 0.5 * vad_duration[0],
         dt_vad_end[0]   + 0.5 * vad_duration[0]],
        z_km, orientation[0:1, :-1].T,
        cmap='twilight_shifted', vmin=0, vmax=360,
    )
    for s in range(1, radar.nsweeps):
        ax0.pcolormesh(
            [dt_vad_start[s] - 0.5 * vad_duration[s],
             dt_vad_end[s]   + 0.5 * vad_duration[s]],
            z_km, orientation[s:s + 1, :-1].T,
            cmap='twilight_shifted', vmin=0, vmax=360,
        )
    cb0 = plt.colorbar(h0, ax=ax0, orientation='horizontal', shrink=0.8)
    cb0.ax.set_xlabel('Wind from direction (°)')
    ax0.set_ylim(0, hmax)
    ax0.set_xlim(dt_min, dt_max)
    ax0.set_ylabel('Distance above radar (km)')
    ax0.set_xlabel('Time (UTC)')
    ax0.set_title(
        f'{INSTRUMENT_NAME}  {time_str}\n'
        f'Wind direction from VAD at elevation {radar.fixed_angle["data"][0]:.1f}°'
    )

    # ── Panel 1: wind speed + barbs ──────────────────────────────────────────
    ax1 = axes[1]
    ax1.xaxis.set_major_formatter(myFmt)
    ax1.set_facecolor('gainsboro')
    ax1.grid(True)

    h1 = ax1.pcolormesh(
        [dt_vad_start[0] - 0.5 * vad_duration[0],
         dt_vad_end[0]   + 0.5 * vad_duration[0]],
        z_km, speed[0:1, :-1].T,
        cmap='viridis', vmin=0, vmax=50,
    )
    for s in range(1, radar.nsweeps):
        ax1.pcolormesh(
            [dt_vad_start[s] - 0.5 * vad_duration[s],
             dt_vad_end[s]   + 0.5 * vad_duration[s]],
            z_km, speed[s:s + 1, :-1].T,
            cmap='viridis', vmin=0, vmax=50,
        )
    cb1 = plt.colorbar(h1, ax=ax1, orientation='horizontal', shrink=0.8)
    cb1.ax.set_xlabel('Wind speed (m s\u207b\u00b9)')

    # Barbs on clock hours in time and ~500 m spacing in height.
    # Find the sweep whose mid-time is closest to each clock hour.
    nlevels = len(zlevels)
    nsweeps = radar.nsweeps
    barb_target_m = 500.0   # metres
    z_stride = max(1, round(barb_target_m / (zlevels[1] - zlevels[0]))) if nlevels > 1 else 1
    dt_mid_arr = np.array(dt_vad_mid)   # (nsweeps,) of datetime objects
    clock_hours = [dt_min + datetime.timedelta(hours=h) for h in range(25)]
    hour_indices = []
    for target in clock_hours:
        diffs = np.array([(dt - target).total_seconds() for dt in dt_mid_arr])
        idx = int(np.argmin(np.abs(diffs)))
        if idx not in hour_indices:
            hour_indices.append(idx)
    hour_indices = np.array(hour_indices)
    x_barb = np.tile(dt_vad_mid, (nlevels, 1)).T   # (nsweeps, nlevels)
    y_barb = np.tile(z_km,       (nsweeps, 1))      # (nsweeps, nlevels)
    ax1.barbs(
        x_barb[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))],
        y_barb[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))],
        u_vel[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))],
        v_vel[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))],
        np.sqrt(
            u_vel[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))]**2 +
            v_vel[np.ix_(hour_indices, np.arange(0, nlevels, z_stride))]**2
        ),
        sizes=dict(emptybarb=0.), length=6,
        cmap='Grays', clim=[0, 50],
    )
    ax1.set_ylim(0, hmax)
    ax1.set_xlim(dt_min, dt_max)
    ax1.set_ylabel('Distance above radar (km)')
    ax1.set_xlabel('Time (UTC)')
    ax1.set_title(
        f'{INSTRUMENT_NAME}  {time_str}\n'
        f'Wind speed from VAD at elevation {radar.fixed_angle["data"][0]:.1f}°'
    )

    figpath_vad = os.path.join(figpath, 'vad')
    os.makedirs(figpath_vad, exist_ok=True)

    figname = f'{INSTRUMENT_NAME}_{location}_{CAMPAIGN}_{datestr}_vad_l1_{product_version}.png'
    if blflag:
        figname = figname.replace('.png', '_bl.png')

    plt.savefig(os.path.join(figpath_vad, figname), dpi=300)
    print(f'Saved VAD plot: {figname}')
    plt.close()


# ── RHI quicklook ─────────────────────────────────────────────────────────────

def make_rhi_plot(ncfile, datestr, figpath, blflag=False):
    """Create a 4-panel RHI quicklook (DBZ, VEL, WIDTH, LDR) for one RHI file."""
    hmax = 4 if blflag else 12
    max_range = 15.0

    with nc4.Dataset(ncfile) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
        location = get_location_from_dataset(ds)

    radar = pyart.io.read_cfradial(ncfile)
    vel_field = radar.fields['VEL']
    vel_min = vel_field.get('field_limit_lower', -10.66)
    vel_max = vel_field.get('field_limit_upper',  10.66)

    valid_sweeps = get_valid_sweep_indices(radar)

    outdir = os.path.join(figpath, 'rhi', datestr)
    os.makedirs(outdir, exist_ok=True)

    for s in valid_sweeps:
        sweep_start = radar.get_start(s)
        sweep_end   = radar.get_end(s)
        dtime_sweep = cftime.num2pydate(radar.time['data'][sweep_start], radar.time['units'])
        dtime_str = dtime_sweep.strftime('%Y%m%d-%H%M%S')
        rhi_az = float(np.nanmean(radar.get_azimuth(s)))

        # Determine whether any ray in this sweep is past-zenith (el > 90°).
        # When past-zenith, pyart plots horizontal distance as negative, so
        # extend the x-axis to [-max_range, max_range].
        sweep_el = radar.elevation['data'][sweep_start:sweep_end + 1]
        past_zenith = bool(np.any(sweep_el > 90.0))
        x_min = -max_range if past_zenith else 0.0

        gf = pyart.correct.GateFilter(radar)
        gf.exclude_below('SNR', -20)
        gf.exclude_transition()
        disp = pyart.graph.RadarDisplay(radar)

        field_specs = [
            ('DBZ',   dict(vmin=-40, vmax=40,   cmap=COLORMAPS['dbz'])),
            ('VEL',   dict(vmin=vel_min, vmax=vel_max, cmap=COLORMAPS['vel'])),
            ('WIDTH', dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1),
                                               vmax=np.sqrt(1e1)),
                           cmap=COLORMAPS['width'])),
            ('LDR',   dict(vmin=-35, vmax=5,    cmap=COLORMAPS['ldr'])),
        ]

        fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=False)
        fig.suptitle(
            f'{RADAR_LABEL} | {RADAR_LOCATION}\n'
            f'RHI Az {rhi_az:.1f}\u00b0  {dtime_sweep.strftime("%Y-%m-%dT%H:%M:%SZ")}',
            fontsize=12, fontweight='bold', y=0.99,
        )
        fig.subplots_adjust(left=0.08, right=0.97, top=0.91, bottom=0.12,
                            wspace=0.35, hspace=0.45)

        for (field, kwargs), ax in zip(field_specs, axes.flatten()):
            disp.plot_rhi(field, s, ax=ax, gatefilter=gf,
                          colorbar_flag=False, title_flag=True,
                          filter_transitions=True,
                          **kwargs)
            ax.set_xlim(x_min, max_range)
            ax.set_ylim(0, hmax)
            ax.set_aspect('equal')
            fld_meta = disp.fields[field]
            label = (f"{fld_meta.get('long_name', field)} "
                     f"({fld_meta.get('units', '')})")
            cb = fig.colorbar(disp.plots[-1], ax=ax,
                              orientation='horizontal', shrink=0.9, pad=0.18)
            cb.set_label(label)

        # Numbered-file suffix: _01, _02, …  extracted from filename if present
        suffix_m = re.search(r'_(0\d+)\.nc$', os.path.basename(ncfile))
        suffix = f'_{suffix_m.group(1)}' if suffix_m else ''

        figname = (
            f'{INSTRUMENT_NAME}_{location}_{CAMPAIGN}_{dtime_str}_'
            f'rhi_az{rhi_az:06.2f}_l1_{product_version}{suffix}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=150)
        plt.close()
        print(f'  Saved RHI plot: {figname}')


def make_rhi_plots_day(datestr, inpath, figpath, blflag=False):
    """Generate RHI quicklooks for all RHI files on a given day."""
    inpath_date = os.path.join(inpath, datestr)
    rhi_files = sorted(glob.glob(os.path.join(inpath_date, f'*{datestr}*rhi*.nc')))

    if not rhi_files:
        print(f'No RHI files found for {datestr} in {inpath_date}')
        return

    print(f'Found {len(rhi_files)} RHI file(s) for {datestr}')
    for ncfile in rhi_files:
        print(f'  Processing: {os.path.basename(ncfile)}')
        try:
            make_rhi_plot(ncfile, datestr, figpath, blflag=blflag)
        except Exception as exc:
            import traceback
            print(f'  Error processing {os.path.basename(ncfile)}: {exc}')
            traceback.print_exc()


def make_ppi_plots_day(datestr, inpath, figpath, blflag=False,
                       use_day_outermost_extent=False, explicit_extent=None,
                       az_offset=0.0, map_alpha=0.65):
    """Generate standard and map PPI quicklooks for all PPI files on a given day."""
    ppi_files = collect_ppi_files(datestr, inpath)
    if not ppi_files:
        print(f'No PPI files found for {datestr}')
        return

    print(f'Found {len(ppi_files)} PPI file(s) for {datestr}')

    day_extent = explicit_extent
    if day_extent is not None:
        print(f'Using explicit map extent: lon [{day_extent[0]}, {day_extent[1]}], '
              f'lat [{day_extent[2]}, {day_extent[3]}]')
    elif use_day_outermost_extent:
        day_extent = _compute_day_outermost_ppi_bounds(ppi_files)
        if day_extent is not None:
            print(f'Day-outermost PPI extent: lon [{day_extent[0]:.4f}, {day_extent[1]:.4f}], '
                  f'lat [{day_extent[2]:.4f}, {day_extent[3]:.4f}]')
        else:
            print('Warning: no parseable bounds found; using per-file fallback extents')

    for ncfile in ppi_files:
        print(f'Processing: {os.path.basename(ncfile)}')
        make_ppi_plot(ncfile, datestr, figpath, blflag=blflag, az_offset=az_offset)
        make_ppi_map_plot(ncfile, datestr, figpath, blflag=blflag,
                          day_extent=day_extent, az_offset=az_offset,
                          map_alpha=map_alpha)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    (datestr, inpath, outpath, blflag, vpt_only, ppi_only, vad_only, rhi_only, no_rhi, no_vad,
     ppi_map_day_extent, explicit_extent, az_offset, map_alpha) = parse_command_line()
    inpath, figpath = setup_paths(datestr, inpath, outpath)

    print('Reading General Quicklooks')
    print(f'  Date:        {datestr}')
    print(f'  Input path:  {inpath}')
    print(f'  Output path: {figpath}')
    print(f'  BL mode:     {blflag}')
    if az_offset:
        print(f'  Az offset:   {az_offset}°')
    if explicit_extent:
        print(f'  Map extent:  lon [{explicit_extent[0]}, {explicit_extent[1]}] '
              f'lat [{explicit_extent[2]}, {explicit_extent[3]}]')
    if vad_only:
        print('  Mode:        VAD only')
    elif vpt_only:
        print('  Mode:        VPT only')
    elif ppi_only:
        print('  Mode:        PPI only')
    elif rhi_only:
        print('  Mode:        RHI only')
    elif no_rhi:
        print('  Mode:        no RHI (VPT + PPI + VAD)')
    print('-' * 60)

    if not ppi_only and not vad_only:
        print('Generating VPT plots...')
        try:
            make_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
        except Exception as e:
            import traceback
            print(f'Error generating VPT plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    if not vpt_only and not vad_only and not ppi_only:
        print('Generating VAD plots...')
        try:
            make_vad_plot_day(datestr, inpath, figpath, blflag=blflag)
        except Exception as e:
            import traceback
            print(f'Error generating VAD plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    if vad_only:
        print('Generating VAD plots...')
        try:
            make_vad_plot_day(datestr, inpath, figpath, blflag=blflag)
        except Exception as e:
            import traceback
            print(f'Error generating VAD plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    if not vpt_only and not vad_only and not rhi_only:
        print('Generating PPI plots...')
        try:
            make_ppi_plots_day(datestr, inpath, figpath, blflag=blflag,
                               use_day_outermost_extent=ppi_map_day_extent,
                               explicit_extent=explicit_extent,
                               az_offset=az_offset, map_alpha=map_alpha)
        except Exception as e:
            import traceback
            print(f'Error generating PPI plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    if not vpt_only and not vad_only and not ppi_only and not no_rhi:
        print('Generating RHI plots...')
        try:
            make_rhi_plots_day(datestr, inpath, figpath, blflag=blflag)
        except Exception as e:
            import traceback
            print(f'Error generating RHI plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    print('Done.')


if __name__ == '__main__':
    main()
