#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_hannington_range_ppi.py

Plot Ka-band reflectivity (DBZ) from today's sector PPIs (azimuth 10–40°) at
the range gate corresponding to the Hannington TV mast
(Hannington transmitting station, Hampshire).

The Hannington TV mast is located at approximately 51.2431°N, 1.1956°W.
Range from the Chilbolton Observatory: ~20.05 km, bearing ~57°.

This script extracts the DBZ value at the Hannington range gate for every
ray in each of today's PPI files and plots reflectivity vs. azimuth.

Usage:
    python plot_hannington_range_ppi.py [-d YYYYMMDD] [-i inpath] [-o outpath]

Arguments:
    -d / --date       Date string YYYYMMDD (default: today)
    -i / --inpath     Directory containing processed CF-Radial files
                      (default: /data/processing/kepler/reading-general)
    -o / --outpath    Output directory for the figure
                      (default: <inpath>/quicklooks)
"""

import argparse
import datetime
import glob
import os
import sys

import cftime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma
import pyart  # registers pyart colormaps (HomeyerRainbow etc.) with matplotlib

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Chilbolton Observatory (from reading-general_project.yml)
RADAR_LAT = 51.14634   # degrees N
RADAR_LON = -1.43846   # degrees E (negative = West)

# Hannington transmitting station (TV mast) on Cottington Hill, Hampshire
# Source: Wikipedia – 51.3080°N 1.2447°W
MAST_LAT = 51.3080    # degrees N
MAST_LON = -1.2447    # degrees E (negative = West)

DEFAULT_INPATH = '/data/processing/kepler/reading-general'

AZIMUTH_MIN = 10.0    # degrees, true north
AZIMUTH_MAX = 40.0    # degrees, true north

SNR_MIN = -20.0       # dB – gate filter threshold

NGATE_HALF = 5        # half-width of range-gate window either side of target


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _haversine_range_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two lat/lon points (degrees)."""
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2))
         * np.sin(dlon / 2) ** 2)
    return R * 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))


def _bearing_deg(lat1, lon1, lat2, lon2):
    """Initial bearing (true north, clockwise, 0–360°) from point 1 to point 2."""
    y = np.sin(np.radians(lon2 - lon1)) * np.cos(np.radians(lat2))
    x = (np.cos(np.radians(lat1)) * np.sin(np.radians(lat2))
         - np.sin(np.radians(lat1)) * np.cos(np.radians(lat2))
         * np.cos(np.radians(lon2 - lon1)))
    return (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0


# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------

def extract_dbz_sector(ncfile, target_range_km, ngate_half,
                       az_min, az_max, snr_min, az_offset=0.0):
    """
    Extract a 2-D DBZ sector (azimuth × range) around target_range_km.

    Parameters
    ----------
    ncfile : str
        Path to a CF-Radial PPI file.
    target_range_km : float
        Centre range in km (Hannington mast distance).
    ngate_half : int
        Number of range gates either side of the centre gate to include.
    az_min, az_max : float
        Azimuth filter limits (degrees, true north).
    snr_min : float
        Minimum SNR threshold; gates below this are masked.

    Returns
    -------
    tuple (az_sorted, rng_window, dbz_sector, scan_time, centre_gate_idx)
        az_sorted    : ndarray (N,)   – azimuths of retained rays (deg)
        rng_window   : ndarray (M,)   – range in km for each gate in window
        dbz_sector   : masked array (N, M) – DBZ values
        scan_time    : datetime
        centre_range : float – actual range of centre gate (km)
    """
    with nc4.Dataset(ncfile) as ds:
        rng_km = ds['range'][:] / 1000.0
        az_raw = ds['azimuth'][:]            # already true-north corrected
        dbz_2d = ds['DBZ'][:]               # (nrays, ngates)
        snr_2d = ds['SNR'][:] if 'SNR' in ds.variables else None
        ant_trans = ds['antenna_transition'][:]

        time_units = ds['time'].units
        time_data = ds['time'][:]
        scan_time = cftime.num2pydate(float(time_data[0]), time_units)

    # SNR gate filter
    if snr_2d is not None:
        dbz_2d = ma.masked_where(snr_2d < snr_min, dbz_2d)

    # Exclude antenna-transition rays
    trans_2d = np.broadcast_to(ant_trans[:, np.newaxis] == 1, dbz_2d.shape)
    dbz_2d = ma.masked_where(trans_2d, dbz_2d)

    # Apply azimuth offset correction (e.g. north_angle refinement)
    az_raw = (az_raw + az_offset) % 360.0

    # Centre gate and range window
    centre_idx = int(np.argmin(np.abs(rng_km - target_range_km)))
    g0 = max(centre_idx - ngate_half, 0)
    g1 = min(centre_idx + ngate_half + 1, len(rng_km))
    rng_window = rng_km[g0:g1]
    centre_range = float(rng_km[centre_idx])

    # Select rays within azimuth window
    az_mod = az_raw % 360.0
    az_mask = (az_mod >= az_min) & (az_mod <= az_max)
    az_sel = az_mod[az_mask]
    dbz_sel = dbz_2d[az_mask, g0:g1]

    # Sort by azimuth
    sort_idx = np.argsort(az_sel)
    az_sel = az_sel[sort_idx]
    dbz_sel = dbz_sel[sort_idx]

    return az_sel, rng_window, dbz_sel, scan_time, centre_range


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_hannington_range_dbz(ppi_files, target_range_km, mast_bearing_deg,
                               datestr, figpath, ngate_half=NGATE_HALF, az_offset=0.0):
    """
    Two-panel 2-D sector plot of DBZ (azimuth × range) around the Hannington
    TV mast range, one panel per PPI scan.

    Parameters
    ----------
    ppi_files : list of str
        Sorted list of CF-Radial PPI files for the day.
    target_range_km : float
        Range of Hannington TV mast (km).
    mast_bearing_deg : float
        True-north bearing to Hannington TV mast (degrees).
    datestr : str
        Date string YYYYMMDD.
    figpath : str
        Output directory.
    ngate_half : int
        Half-width of range-gate window (gates either side of centre).
    """
    all_sectors = []  # list of (az, rng, dbz_2d, label, centre_range)
    for ncfile in ppi_files:
        try:
            az, rng, dbz, scan_time, centre_range = extract_dbz_sector(
                ncfile, target_range_km, ngate_half,
                AZIMUTH_MIN, AZIMUTH_MAX, SNR_MIN, az_offset=az_offset)
        except Exception as exc:
            print(f'Warning: could not read {os.path.basename(ncfile)}: {exc}')
            continue
        if len(az) == 0:
            print(f'  No rays in {AZIMUTH_MIN}–{AZIMUTH_MAX}° in '
                  f'{os.path.basename(ncfile)}')
            continue
        label = scan_time.strftime('%H:%M:%S UTC')
        all_sectors.append((az, rng, dbz, label, centre_range))
        gate_spacing_m = float(np.diff(rng).mean() * 1000) if len(rng) > 1 else 0
        print(f'  {os.path.basename(ncfile)}: {len(az)} rays, '
              f'{len(rng)} gates ({gate_spacing_m:.0f} m spacing), '
              f'centre {centre_range:.3f} km, time {label}')

    if not all_sectors:
        print('No data found to plot.')
        return

    nscan = len(all_sectors)
    dbz_cmap = matplotlib.colormaps['HomeyerRainbow']
    dbz_norm = matplotlib.colors.Normalize(vmin=-40, vmax=40)

    fig, axes = plt.subplots(1, nscan, figsize=(8 * nscan, 6),
                             sharey=True, constrained_layout=True)
    if nscan == 1:
        axes = [axes]

    date_label = datetime.datetime.strptime(datestr, '%Y%m%d').strftime('%Y-%m-%d')
    gate_half_m = ngate_half * 60  # assumes 60 m gate spacing
    fig.suptitle(
        f'Ka-band reflectivity near Hannington TV mast | {date_label}\n'
        f'Range window: {target_range_km:.2f} km ± {gate_half_m} m  '
        f'| Azimuth {AZIMUTH_MIN}–{AZIMUTH_MAX}°',
        fontsize=12,
    )

    for ax, (az, rng, dbz, label, centre_range) in zip(axes, all_sectors):
        # pcolormesh needs cell edges; build half-cell-padded edges
        def _edges(arr):
            if len(arr) < 2:
                d = 0.03  # fallback half-width (km)
                return np.array([arr[0] - d, arr[0] + d])
            d = np.diff(arr)
            lo = np.concatenate([[arr[0] - d[0] / 2], arr[:-1] + d / 2])
            hi = np.array([arr[-1] + d[-1] / 2])
            return np.concatenate([lo, hi])

        az_edges = _edges(az)
        rng_edges = _edges(rng)

        # pcolormesh: x=azimuth, y=range, C=DBZ (dims must be [nrng, naz])
        mesh = ax.pcolormesh(az_edges, rng_edges, dbz.T,
                             cmap=dbz_cmap, norm=dbz_norm,
                             shading='flat', rasterized=True)

        # Hannington range line
        ax.axhline(centre_range, color='white', lw=1.2, ls='--',
                   label=f'Hannington range ({centre_range:.2f} km)')

        # Hannington bearing if it falls in sector
        if AZIMUTH_MIN <= mast_bearing_deg <= AZIMUTH_MAX:
            ax.axvline(mast_bearing_deg, color='white', lw=1.2, ls=':',
                       label=f'Hannington bearing ({mast_bearing_deg:.1f}°)')

        ax.set_xlim(AZIMUTH_MIN, AZIMUTH_MAX)
        ax.set_xlabel('Azimuth from true north (°)')
        ax.set_title(label, fontsize=11)
        ax.legend(fontsize=8, loc='upper right',
                  facecolor='0.2', labelcolor='white', framealpha=0.7)
        ax.grid(True, alpha=0.25, color='white', lw=0.5)

    axes[0].set_ylabel('Range (km)')

    # Shared colorbar
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=dbz_norm, cmap=dbz_cmap),
                      ax=axes, orientation='vertical', shrink=0.85, pad=0.02)
    cb.set_label('Equivalent reflectivity factor DBZ (dBZ)')

    os.makedirs(figpath, exist_ok=True)
    figname = f'kepler_reading-general_{datestr}_dbz_hannington_range.png'
    outfile = os.path.join(figpath, figname)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outfile}')


def plot_hannington_azimuth_lines(ppi_files, target_range_km, mast_bearing_deg,
                                  datestr, figpath, ngate_half=NGATE_HALF, az_offset=0.0):
    """
    Line plot of DBZ vs azimuth, one line per range gate in the window,
    one panel per PPI scan.

    Parameters
    ----------
    ppi_files : list of str
    target_range_km : float
    mast_bearing_deg : float
    datestr : str
    figpath : str
    ngate_half : int
    """
    all_sectors = []
    for ncfile in ppi_files:
        try:
            az, rng, dbz, scan_time, centre_range = extract_dbz_sector(
                ncfile, target_range_km, ngate_half,
                AZIMUTH_MIN, AZIMUTH_MAX, SNR_MIN, az_offset=az_offset)
        except Exception as exc:
            print(f'Warning: could not read {os.path.basename(ncfile)}: {exc}')
            continue
        if len(az) == 0:
            continue
        all_sectors.append((az, rng, dbz, scan_time.strftime('%H:%M:%S UTC'), centre_range))

    if not all_sectors:
        print('No data found for line plot.')
        return

    nscan = len(all_sectors)
    ngate = len(all_sectors[0][1])  # number of gates in window

    # Colour each gate line by its offset from the centre gate
    gate_cmap = matplotlib.colormaps['coolwarm'].resampled(ngate)

    fig, axes = plt.subplots(1, nscan, figsize=(9 * nscan, 5),
                             sharey=True, constrained_layout=True)
    if nscan == 1:
        axes = [axes]

    date_label = datetime.datetime.strptime(datestr, '%Y%m%d').strftime('%Y-%m-%d')
    gate_half_m = ngate_half * 60
    fig.suptitle(
        f'Ka-band reflectivity vs azimuth near Hannington TV mast | {date_label}\n'
        f'Range window: {target_range_km:.2f} km ± {gate_half_m} m  '
        f'| Azimuth {AZIMUTH_MIN}–{AZIMUTH_MAX}°',
        fontsize=12,
    )

    for ax, (az, rng, dbz, label, centre_range) in zip(axes, all_sectors):
        centre_idx = int(np.argmin(np.abs(rng - centre_range)))
        for gi in range(len(rng)):
            offset_gates = gi - centre_idx
            offset_m = int(round(offset_gates * np.diff(rng).mean() * 1000)) if len(rng) > 1 else 0
            gate_label = (f'{rng[gi]:.3f} km'
                          + (f' ({offset_m:+d} m)' if offset_gates != 0 else ' [centre]'))
            col = gate_cmap(gi / max(ngate - 1, 1))
            lw  = 1.8 if offset_gates == 0 else 1.0
            zord = 3 if offset_gates == 0 else 2
            ax.plot(az, dbz[:, gi], color=col, lw=lw, zorder=zord, label=gate_label)

        # Mast bearing marker
        ax.axvline(mast_bearing_deg, color='k', lw=1.2, ls='--',
                   label=f'Hannington bearing ({mast_bearing_deg:.1f}°)', zorder=4)

        ax.set_xlim(AZIMUTH_MIN, AZIMUTH_MAX)
        ax.set_ylim(-40, 40)
        ax.set_xlabel('Azimuth from true north (°)')
        ax.set_title(label, fontsize=11)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=8, loc='upper left', title='Range gate')

    axes[0].set_ylabel('Equivalent reflectivity factor DBZ (dBZ)')

    os.makedirs(figpath, exist_ok=True)
    figname = f'kepler_reading-general_{datestr}_dbz_hannington_azimuth_lines.png'
    outfile = os.path.join(figpath, figname)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {outfile}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot Ka-band DBZ at Hannington TV mast range from sector PPIs.')
    parser.add_argument('-d', '--date',
                        default=datetime.datetime.now().strftime('%Y%m%d'),
                        help='Date YYYYMMDD (default: today)')
    parser.add_argument('-i', '--inpath',
                        default=DEFAULT_INPATH,
                        help='Input directory of processed CF-Radial files')
    parser.add_argument('-o', '--outpath',
                        default=None,
                        help='Output directory for figure (default: <inpath>/quicklooks)')
    parser.add_argument('-n', '--ngate-half', type=int, default=NGATE_HALF,
                        help=f'Gates either side of Hannington range (default: {NGATE_HALF})')
    parser.add_argument('--az-offset', type=float, default=0.0,
                        help='Azimuth correction to add to stored azimuths (degrees, default: 0)')
    return parser.parse_args()


def main():
    args = parse_args()
    datestr = args.date
    inpath = args.inpath
    figpath = args.outpath or os.path.join(inpath, 'quicklooks')

    # Compute Hannington geometry
    target_range_km = _haversine_range_km(RADAR_LAT, RADAR_LON, MAST_LAT, MAST_LON)
    mast_bearing_deg = _bearing_deg(RADAR_LAT, RADAR_LON, MAST_LAT, MAST_LON)
    print(f'Hannington TV mast: range = {target_range_km:.3f} km, '
          f'bearing = {mast_bearing_deg:.2f}° (true north)')

    # Collect PPI files for the day
    ppi_files = sorted(glob.glob(
        os.path.join(inpath, datestr, '*ppi*.nc')))
    if not ppi_files:
        print(f'No PPI files found for {datestr} in {inpath}')
        sys.exit(1)
    print(f'Found {len(ppi_files)} PPI file(s) for {datestr}')

    plot_hannington_range_dbz(ppi_files, target_range_km, mast_bearing_deg,
                               datestr, figpath, ngate_half=args.ngate_half,
                               az_offset=args.az_offset)
    plot_hannington_azimuth_lines(ppi_files, target_range_km, mast_bearing_deg,
                                  datestr, figpath, ngate_half=args.ngate_half,
                                  az_offset=args.az_offset)


if __name__ == '__main__':
    main()
