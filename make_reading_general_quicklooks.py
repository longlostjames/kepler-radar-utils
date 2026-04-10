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


# ── Command-line parsing ───────────────────────────────────────────────────────

def parse_command_line():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], 'd:i:o:b',
            ['date=', 'inpath=', 'outpath=']
        )
    except getopt.GetoptError as err:
        print(f'Error: {err}')
        print('Usage: python make_reading_general_quicklooks.py -d YYYYMMDD '
              '[-i inpath] [-o outpath] [-b]')
        sys.exit(2)

    datestr = datetime.datetime.now().strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False

    for option, argument in opts:
        if option in ('-d', '--date'):
            datestr = argument
        elif option in ('-i', '--inpath'):
            inpath = argument
        elif option in ('-o', '--outpath'):
            outpath = argument
        elif option == '-b':
            blflag = True
        else:
            print(f'Unhandled option: {option}')
            sys.exit(2)

    return datestr, inpath, outpath, blflag


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


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    datestr, inpath, outpath, blflag = parse_command_line()
    inpath, figpath = setup_paths(datestr, inpath, outpath)

    print('Reading General Quicklooks')
    print(f'  Date:        {datestr}')
    print(f'  Input path:  {inpath}')
    print(f'  Output path: {figpath}')
    print(f'  BL mode:     {blflag}')
    print('-' * 60)

    print('Generating VPT plots...')
    try:
        make_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
    except Exception as e:
        import traceback
        print(f'Error generating VPT plots: {e}')
        traceback.print_exc()

    print('-' * 60)
    print('Done.')


if __name__ == '__main__':
    main()
