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
import cmocean

from kepler_utils import get_valid_sweep_indices

# ── Configuration ─────────────────────────────────────────────────────────────

CAMPAIGN = 'reading-general'
INSTRUMENT_NAME = 'reading-mobile-ka-band-radar-1'
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
            ['date=', 'inpath=', 'outpath=', 'vpt-only', 'ppi-only', 'ppi-map-day-extent']
        )
    except getopt.GetoptError as err:
        print(f'Error: {err}')
        print('Usage: python make_reading_general_quicklooks.py -d YYYYMMDD '
              '[-i inpath] [-o outpath] [-b] [--vpt-only] [--ppi-only] [--ppi-map-day-extent]')
        sys.exit(2)

    datestr = datetime.datetime.now().strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False
    vpt_only = False
    ppi_only = False
    ppi_map_day_extent = False

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
        elif option == '--ppi-map-day-extent':
            ppi_map_day_extent = True
        else:
            print(f'Unhandled option: {option}')
            sys.exit(2)

    return datestr, inpath, outpath, blflag, vpt_only, ppi_only, ppi_map_day_extent


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

def make_vpt_plot_day(datestr, inpath, figpath, blflag=False):
    """Generate a daily VPT time-height plot (DBZ, VEL, WIDTH, LDR)."""

    hmax = 4 if blflag else 12

    inpath_date = os.path.join(inpath, datestr)
    vpt_files = sorted(glob.glob(os.path.join(inpath_date, f'*{datestr}*vpt*.nc')))

    if not vpt_files:
        print(f'No VPT files found for {datestr} in {inpath_date}')
        return

    vpt_file = vpt_files[0]
    print(f'Reading VPT file: {vpt_file}')

    with nc4.Dataset(vpt_file) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')

    RadarDS_VPT = pyart.io.read_cfradial(vpt_file)
    nsweeps = RadarDS_VPT.nsweeps

    # Try previous day for continuity
    prev_date = (datetime.datetime.strptime(datestr, '%Y%m%d')
                 - datetime.timedelta(days=1))
    prevstr = prev_date.strftime('%Y%m%d')
    RadarDS_VPT_prev = None
    nsweeps_prev = 0
    try:
        prev_files = sorted(glob.glob(
            os.path.join(inpath, prevstr, f'*{prevstr}*vpt*.nc')))
        if prev_files:
            RadarDS_VPT_prev = pyart.io.read_cfradial(prev_files[0])
            nsweeps_prev = RadarDS_VPT_prev.nsweeps
    except Exception:
        pass

    # Time axis limits for the day
    dtime = cftime.num2pydate(RadarDS_VPT.time['data'], RadarDS_VPT.time['units'])
    dt_min = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    dt_max = dt_min + datetime.timedelta(days=1)
    time_str = dtime[0].strftime('%Y-%m-%d')

    fig, axes = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=1/72, h_pad=1/72,
                                    hspace=0.2, wspace=0.2)

    def _add_sweep(radar, sweep_idx, first=False):
        sweep = radar.extract_sweeps([sweep_idx])
        gf = pyart.correct.GateFilter(sweep)
        gf.exclude_below('SNR', -20)
        disp = pyart.graph.RadarDisplay(sweep)

        kw_first = dict(title=f'{pyart.graph.common.generate_radar_name(sweep)} '
                               f'{time_str}\n'
                               f'{pyart.graph.common.generate_field_name(sweep, "DBZ")}',
                        colorbar_orient='horizontal')
        kw_rest = dict(colorbar_flag=False, title_flag=False)

        disp.plot_vpt('DBZ', ax=axes[0], time_axis_flag=True, edges=False,
                      gatefilter=gf, vmin=-40, vmax=40, norm=None,
                      filter_transitions=True, cmap=COLORMAPS['dbz'],
                      **(kw_first if first else kw_rest))
        disp.plot_vpt('VEL', ax=axes[1], time_axis_flag=True, edges=False,
                      gatefilter=gf, vmin=-5, vmax=5, norm=None,
                      cmap=COLORMAPS['vel'],
                      **(kw_first if first else kw_rest))
        disp.plot_vpt('WIDTH', ax=axes[2], time_axis_flag=True, edges=False,
                      gatefilter=gf,
                      norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),
                                          vmax=np.sqrt(1e1)),
                      cmap=COLORMAPS['width'],
                      **(kw_first if first else kw_rest))
        disp.plot_vpt('LDR', ax=axes[3], time_axis_flag=True, edges=False,
                      gatefilter=gf, vmin=-35, vmax=5, norm=None,
                      cmap=COLORMAPS['ldr'],
                      **(kw_first if first else kw_rest))

    # Add previous day's last sweep for continuity
    if nsweeps_prev > 0:
        _add_sweep(RadarDS_VPT_prev, nsweeps_prev - 1, first=False)

    # Add all sweeps for the current day
    for s in range(nsweeps):
        _add_sweep(RadarDS_VPT, s, first=(s == 0 and nsweeps_prev == 0))

    for ax in axes:
        ax.set_xlim(dt_min, dt_max)
        ax.set_ylim(0, hmax)
        ax.set_xlabel('Time (UTC)')
        ax.grid(True)

    figpath_vpt = os.path.join(figpath, 'vpt')
    os.makedirs(figpath_vpt, exist_ok=True)

    figname = (f'{INSTRUMENT_NAME}_cao_{CAMPAIGN}_{datestr}_vpt_l1'
               f'_{product_version}.png')
    if blflag:
        figname = figname.replace('.png', '_bl.png')

    plt.savefig(os.path.join(figpath_vpt, figname), dpi=150)
    print(f'Saved: {figname}')
    plt.close()


# ── PPI quicklooks ────────────────────────────────────────────────────────────

def make_ppi_plot(ncfile, datestr, figpath, blflag=False):
    """Create a standard (Cartesian x/y) PPI quicklook for each sweep in ncfile."""
    max_range = 20.0 if blflag else 30.0

    with nc4.Dataset(ncfile) as ds:
        product_version = getattr(ds, 'product_version', 'v1.0.0')
        location = get_location_from_dataset(ds)

    radar = pyart.io.read_cfradial(ncfile)
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
             dict(vmin=-40, vmax=40, cmap='HomeyerRainbow')),
            ('LDR', 'Radar Linear Depolarization Ratio', axes[0, 1],
             dict(vmin=-35, vmax=5, cmap='SpectralExtended')),
            ('VEL', 'Radial Velocity', axes[1, 0],
             dict(vmin=vel_min, vmax=vel_max, cmap='balance')),
            ('WIDTH', 'Radar Doppler Spectrum Width', axes[1, 1],
             dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
                  cmap='SpectralExtended')),
        ]

        for field, title, ax, kwargs in field_specs:
            display.plot_ppi(field, s, ax=ax, colorbar_flag=True, title_flag=False,
                             gatefilter=gf, filter_transitions=True,
                             mask_outside=True, **kwargs)
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
        plt.savefig(os.path.join(outdir, figname), dpi=150)
        plt.close()
        print(f'Saved PPI plot: {figname}')


def make_ppi_map_plot(ncfile, datestr, figpath, blflag=False,
                      day_extent=None, colorbar_shrink=0.86):
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
        if day_extent is not None:
            lon_min, lon_max, lat_min, lat_max = day_extent
        elif _bounds is not None:
            lon_min, lon_max, lat_min, lat_max = _bounds
        else:
            lon_min, lon_max = radar_lon - 0.50, radar_lon + 0.50
            lat_min, lat_max = radar_lat - 0.35, radar_lat + 0.35

        fig, axes = plt.subplots(
            2, 2, figsize=(16, 16),
            subplot_kw=dict(projection=ccrs.OSGB()),
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
            display.plot_ppi_map(
                field, s,
                projection=ccrs.OSGB(), ax=panel_ax,
                lon_lines=None, lat_lines=None,
                colorbar_flag=False, title_flag=False,
                filter_transitions=True, gatefilter=None,
                **kwargs,
            )
            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max],
                                 crs=ccrs.PlateCarree())

            if ctx is not None:
                try:
                    ctx.add_basemap(panel_ax, zoom=11,
                                    source=ctx.providers.OpenTopoMap,
                                    crs='EPSG:27700', reset_extent=False)
                except Exception:
                    try:
                        ctx.add_basemap(panel_ax, zoom=11,
                                        source=ctx.providers.OpenStreetMap.Mapnik,
                                        crs='EPSG:27700', reset_extent=False)
                    except Exception as e:
                        print(f'Warning: could not add basemap tiles: {e}')

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

            mappable = display.plots[-1]
            field_meta = display.fields[field]
            label = f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})"
            cb = fig.colorbar(mappable, ax=panel_ax, orientation='vertical',
                              shrink=colorbar_shrink, pad=0.02)
            cb.set_label(label)
            cb.ax.tick_params(labelsize=8)

        figname = (
            f'{INSTRUMENT_NAME}_{location}_{CAMPAIGN}_{dtime_str}_'
            f'ppi-map_el{ppi_el:0.2f}_l1_{product_version}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=150)
        plt.close()
        print(f'Saved PPI map plot: {figname}')


def make_ppi_plots_day(datestr, inpath, figpath, blflag=False,
                       use_day_outermost_extent=False):
    """Generate standard and map PPI quicklooks for all PPI files on a given day."""
    ppi_files = collect_ppi_files(datestr, inpath)
    if not ppi_files:
        print(f'No PPI files found for {datestr}')
        return

    print(f'Found {len(ppi_files)} PPI file(s) for {datestr}')

    day_extent = None
    if use_day_outermost_extent:
        day_extent = _compute_day_outermost_ppi_bounds(ppi_files)
        if day_extent is not None:
            print(f'Day-outermost PPI extent: lon [{day_extent[0]:.4f}, {day_extent[1]:.4f}], '
                  f'lat [{day_extent[2]:.4f}, {day_extent[3]:.4f}]')
        else:
            print('Warning: no parseable bounds found; using per-file fallback extents')

    for ncfile in ppi_files:
        print(f'Processing: {os.path.basename(ncfile)}')
        make_ppi_plot(ncfile, datestr, figpath, blflag=blflag)
        make_ppi_map_plot(ncfile, datestr, figpath, blflag=blflag,
                          day_extent=day_extent)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    datestr, inpath, outpath, blflag, vpt_only, ppi_only, ppi_map_day_extent = parse_command_line()
    inpath, figpath = setup_paths(datestr, inpath, outpath)

    print('Reading General Quicklooks')
    print(f'  Date:        {datestr}')
    print(f'  Input path:  {inpath}')
    print(f'  Output path: {figpath}')
    print(f'  BL mode:     {blflag}')
    if vpt_only:
        print('  Mode:        VPT only')
    elif ppi_only:
        print('  Mode:        PPI only')
    print('-' * 60)

    if not ppi_only:
        print('Generating VPT plots...')
        try:
            make_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
        except Exception as e:
            import traceback
            print(f'Error generating VPT plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    if not vpt_only:
        print('Generating PPI plots...')
        try:
            make_ppi_plots_day(datestr, inpath, figpath, blflag=blflag,
                               use_day_outermost_extent=ppi_map_day_extent)
        except Exception as e:
            import traceback
            print(f'Error generating PPI plots: {e}')
            traceback.print_exc()
        print('-' * 60)

    print('Done.')


if __name__ == '__main__':
    main()
