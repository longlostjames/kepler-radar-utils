#!/usr/bin/env python

import getopt, sys, os
import re

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
import contextily as ctx

import pandas as pd

import cftime

# Import kepler utilities
from kepler_utils import get_valid_sweep_indices

VERSION = 0.1
TRACKING_TAG = 'AMOF_20230201132601'
CAMPAIGN = 'ccrest-m'
DEFAULT_DATA_VERSION = '1.0.1'

NCAS_RADAR_PATH = '/gws/nopw/j04/ncas_radar_vol1'
NCAS_OBS_PROC_PATH = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-mobile-ka-band-radar-1'
NCAS_INSTRUMENT_NAME = 'ncas-mobile-ka-band-radar-1'
CCREST_RHI_TARGETS = {
    'ccrest-1': 270.0,
    'ccrest-2': 246.0,
}
CCREST_RHI_TOLERANCE_DEG = 10.0

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')


def parse_command_line():
    """Parse command line arguments."""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "d:i:o:b",
            [
                "date=",
                "inpath=",
                "outpath=",
                "skip-all-transition",
                "vpt-only",
                "ppi-only",
                "pointing-only",
                "ppi-map-day-extent",
                "ppi-map-colorbar-shrink=",
            ],
        )
    except getopt.GetoptError as err:
        print(f"Error: {err}")
        print(
            "Usage: python make_ccrest_quicklooks.py -d YYYYMMDD -i input_path -o output_path [-b] [--skip-all-transition] [--vpt-only] [--ppi-only] [--pointing-only] [--ppi-map-day-extent] [--ppi-map-colorbar-shrink VALUE]"
        )
        sys.exit(2)

    run_date = datetime.datetime.now().strftime('%Y%m%d')
    inpath = None
    outpath = None
    blflag = False
    skip_all_transition = False
    vpt_only = False
    ppi_only = False
    pointing_only = False
    ppi_map_day_extent = False
    ppi_map_colorbar_shrink = 0.86

    for option, argument in opts:
        if option in ("-d", "--date"):
            run_date = argument
        elif option in ("-i", "--inpath"):
            inpath = argument
        elif option in ("-o", "--outpath"):
            outpath = argument
        elif option == "-b":
            blflag = True
        elif option == "--skip-all-transition":
            skip_all_transition = True
        elif option == "--vpt-only":
            vpt_only = True
        elif option == "--ppi-only":
            ppi_only = True
        elif option == "--pointing-only":
            pointing_only = True
        elif option == "--ppi-map-day-extent":
            ppi_map_day_extent = True
        elif option == "--ppi-map-colorbar-shrink":
            try:
                ppi_map_colorbar_shrink = float(argument)
            except ValueError:
                print(f"Error: invalid --ppi-map-colorbar-shrink value '{argument}'")
                sys.exit(2)
        else:
            print(f"Unhandled option: {option}")
            sys.exit(2)

    if not (0.2 <= ppi_map_colorbar_shrink <= 1.0):
        print("Error: --ppi-map-colorbar-shrink must be between 0.2 and 1.0")
        sys.exit(2)

    exclusive = sum([vpt_only, ppi_only, pointing_only])
    if exclusive > 1:
        print("Error: --vpt-only, --ppi-only and --pointing-only are mutually exclusive")
        sys.exit(2)

    return (
        run_date,
        inpath,
        outpath,
        blflag,
        skip_all_transition,
        vpt_only,
        ppi_only,
        pointing_only,
        ppi_map_day_extent,
        ppi_map_colorbar_shrink,
    )


def setup_paths(run_date):
    """Set up input and figure output paths."""
    base_campaign_path = os.path.join(NCAS_OBS_PROC_PATH, CAMPAIGN)

    # inpath = os.path.join(NCAS_RADAR_PATH, 'cjw', 'projects', CAMPAIGN, 'kepler', 'L1a')
    inpath = os.path.join(base_campaign_path, f'L1_v{DEFAULT_DATA_VERSION}')

    # If the preferred default version has no folder for the requested date,
    # fall back to the newest available L1_v* directory that does.
    if run_date:
        preferred_date_path = os.path.join(inpath, run_date)
        if not os.path.isdir(preferred_date_path):
            version_dirs = [
                d for d in glob.glob(os.path.join(base_campaign_path, 'L1_v*')) if os.path.isdir(d)
            ]

            def version_key(path):
                version_str = os.path.basename(path).replace('L1_v', '')
                parts = []
                for token in version_str.split('.'):
                    try:
                        parts.append(int(token))
                    except ValueError:
                        parts.append(0)
                return tuple(parts)

            version_dirs.sort(key=version_key, reverse=True)

            for candidate in version_dirs:
                if os.path.isdir(os.path.join(candidate, run_date)):
                    inpath = candidate
                    break

    figpath = os.path.join(inpath, 'quicklooks')

    return inpath, figpath


def _sanitize_location(value):
    """Normalize location strings for filename-safe NCAS tokens."""
    if value is None:
        return 'unknown'

    text = str(value).strip().lower()
    text = text.replace(' ', '-')
    text = re.sub(r'[^a-z0-9_-]+', '', text)
    return text if text else 'unknown'


def get_location_from_dataset(ds):
    """Extract location from NetCDF global attributes with sensible fallbacks."""
    for attr in ('location', 'platform_location', 'site_name'):
        if hasattr(ds, attr):
            return _sanitize_location(getattr(ds, attr))

    # Some files only carry a verbose platform description.
    if hasattr(ds, 'platform'):
        platform_text = str(getattr(ds, 'platform')).lower()
        if 'cao' in platform_text:
            return 'cao'

    return 'unknown'


def angular_difference_deg(angle1, angle2):
    """Return the smallest angular separation in degrees."""
    return abs((angle1 - angle2 + 180.0) % 360.0 - 180.0)


def classify_ccrest_rhi_file(ncfile):
    """Classify a CCREST RHI file as ccrest-1 or ccrest-2 using fixed_angle."""
    try:
        with nc4.Dataset(ncfile) as ds:
            if 'fixed_angle' not in ds.variables:
                return None

            fixed_angle = float(ds.variables['fixed_angle'][0])
    except Exception as exc:
        print(f"Skipping {ncfile}: unable to read fixed_angle ({exc})")
        return None

    best_label = None
    best_diff = None
    for label, target_angle in CCREST_RHI_TARGETS.items():
        diff = angular_difference_deg(fixed_angle, target_angle)
        if diff <= CCREST_RHI_TOLERANCE_DEG and (best_diff is None or diff < best_diff):
            best_label = label
            best_diff = diff

    if best_label is None:
        print(f"Skipping {ncfile}: fixed_angle={fixed_angle:.2f} does not match CCREST RHI targets")
        return None

    print(f"Classified {os.path.basename(ncfile)} as {best_label} from fixed_angle={fixed_angle:.2f}")
    return best_label


def collect_ccrest_rhi_files(datestr, inpath):
    """Group CCREST RHI files by target label using fixed_angle metadata."""
    inpath_date = os.path.join(inpath, datestr)
    os.chdir(inpath_date)

    grouped_files = {label: [] for label in CCREST_RHI_TARGETS}
    candidate_files = sorted(glob.glob('*.nc'))

    for filename in candidate_files:
        lower_name = filename.lower()
        if 'vpt' in lower_name or 'vad' in lower_name or 'pointing' in lower_name:
            continue

        filepath = os.path.join(inpath_date, filename)
        label = classify_ccrest_rhi_file(filepath)
        if label is not None:
            grouped_files[label].append(filepath)

    return grouped_files


def collect_ccrest_pointing_files(datestr, inpath):
    """Collect CCREST pointing files for a given day."""
    inpath_date = os.path.join(inpath, datestr)
    os.chdir(inpath_date)
    return sorted([os.path.join(inpath_date, f) for f in glob.glob('*pointing*.nc')])


def collect_ccrest_ppi_files(datestr, inpath):
    """Collect CCREST PPI files for a given day."""
    inpath_date = os.path.join(inpath, datestr)
    os.chdir(inpath_date)
    return sorted([os.path.join(inpath_date, f) for f in glob.glob('*ppi*.nc')])


def _plot_ccrest_rhi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max,
                            hmax, xmin, xmax):
    """Plot the standard 2x2 CCREST RHI field set."""
    from matplotlib import colors

    field_specs = [
        ("DBZ", "Equivalent Reflectivity Factor", axes[0, 0], dict(vmin=-40, vmax=40, cmap='HomeyerRainbow')),
        ("LDR", "Radar Linear Depolarization Ratio", axes[0, 1], dict(vmin=-35, vmax=5, cmap='SpectralExtended')),
        ("VEL", "Radial Velocity", axes[1, 0], dict(vmin=vel_min, vmax=vel_max, cmap='balance')),
        (
            "WIDTH",
            "Radar Doppler Spectrum Width",
            axes[1, 1],
            dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
                 cmap='SpectralExtended'),
        ),
    ]

    for field, short_title, ax, kwargs in field_specs:
        display.plot_rhi(
            field,
            ax=ax,
            sweep=sweep_idx,
            colorbar_flag=False,
            title_flag=False,
            gatefilter=gatefilter,
            filter_transitions=True,
            reverse_xaxis=True,
            **kwargs,
        )
        _setup_ccrest_rhi_axes(ax, hmax, xmin, xmax)
        ax.set_title(short_title, fontsize=11, pad=8)
        mappable = display.plots[-1]
        field_meta = display.fields[field]
        label = f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})"
        cb = ax.get_figure().colorbar(
            mappable,
            ax=ax,
            orientation='horizontal',
            shrink=0.92,
            pad=0.12,
        )
        cb.set_label(label)
        cb.ax.tick_params(labelsize=9)


def _setup_ccrest_rhi_axes(ax, hmax, xmin, xmax):
    """Apply consistent axis styling to CCREST RHI subplots."""
    ax.set_ylim(0, hmax)
    ax.set_xlim(xmin, xmax)
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_aspect('equal')


def _build_ccrest_rhi_gatefilter(radar, sweep_idx, snr_threshold=-8.0):
    """Build a gate filter that also masks likely skipped-ray artifacts."""
    gatefilter = pyart.correct.GateFilter(radar)
    gatefilter.exclude_below('SNR', snr_threshold)

    nrays = radar.nrays
    ngates = radar.ngates
    sweep_slice = radar.get_slice(sweep_idx)
    ray_bad = np.zeros(nrays, dtype=bool)
    scan_mask = np.ones(sweep_slice.stop - sweep_slice.start, dtype=bool)

    # Respect antenna transition flags if present.
    if hasattr(radar, 'antenna_transition') and radar.antenna_transition is not None:
        at_data = radar.antenna_transition.get('data')
        if at_data is not None:
            scan_mask = at_data[sweep_slice] == 0
            ray_bad[sweep_slice] |= ~scan_mask

    # Mask rays around abrupt elevation jumps (typical signature of skipped rays).
    sweep_elev = radar.get_elevation(sweep_idx)
    scan_elev = sweep_elev[scan_mask]
    if scan_elev.size > 3:
        elev_diff = np.abs(np.diff(scan_elev))
        finite_diff = elev_diff[np.isfinite(elev_diff)]
        if finite_diff.size > 0:
            median_step = np.median(finite_diff)
            jump_threshold = max(0.5, 4.0 * median_step)
            jump_idx = np.where(elev_diff > jump_threshold)[0]
            if jump_idx.size > 0:
                local_bad = np.zeros(sweep_elev.shape[0], dtype=bool)
                scan_ray_idx = np.flatnonzero(scan_mask)
                local_bad[scan_ray_idx[jump_idx]] = True
                local_bad[scan_ray_idx[np.clip(jump_idx + 1, 0, scan_ray_idx.size - 1)]] = True
                ray_bad[sweep_slice] |= local_bad

    n_bad = int(np.sum(ray_bad[sweep_slice]))
    if n_bad > 0:
        bad_gate_mask = np.zeros((nrays, ngates), dtype=bool)
        bad_gate_mask[ray_bad, :] = True
        gatefilter.exclude_gates(bad_gate_mask)
        print(f"Masked {n_bad} suspect ray(s) in sweep {sweep_idx} (transitions/skipped rays)")

    return gatefilter


def _centers_to_edges(values):
    """Convert 1D center coordinates to edge coordinates."""
    centers = np.asarray(values, dtype=float)
    if centers.size == 0:
        return np.array([0.0, 1.0], dtype=float)
    if centers.size == 1:
        delta = 1.0
        return np.array([centers[0] - 0.5 * delta, centers[0] + 0.5 * delta], dtype=float)

    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def _plot_vpt_field_manual(radar_sweep, field, ax, gatefilter=None, vmin=None, vmax=None,
                           norm=None, cmap=None, title=None, colorbar_flag=True,
                           colorbar_orient='horizontal'):
    """Plot VPT field via pcolormesh using explicit edges to avoid shape mismatch bugs."""
    field_data = ma.array(radar_sweep.fields[field]['data'], copy=True)
    if gatefilter is not None:
        field_data = ma.masked_where(gatefilter.gate_excluded, field_data)

    dtime = cftime.num2pydate(radar_sweep.time['data'], radar_sweep.time['units'])
    time_centers = mdates.date2num(dtime)
    time_edges = _centers_to_edges(time_centers)

    range_km = radar_sweep.range['data'] / 1000.0
    range_edges = _centers_to_edges(range_km)

    mesh = ax.pcolormesh(
        time_edges,
        range_edges,
        field_data.transpose(),
        shading='flat',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=norm,
    )

    ax.xaxis_date()
    if title is not None:
        ax.set_title(title)

    if colorbar_flag:
        cb = ax.figure.colorbar(mesh, ax=ax, orientation=colorbar_orient)
        field_meta = radar_sweep.fields[field]
        cb.set_label(f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})")
        cb.ax.tick_params(labelsize=8)

    return mesh


def make_ccrest_rhi_plot(ncfile,figpath,ccrest_az,blflag=False,skip_all_transition=False):

    if blflag:
        hmax = 4;
        xmin = 0;
        xmax = 20;
    else:
        hmax = 12;
        xmin = 0;
        xmax = 25;
    
    DS = nc4.Dataset(ncfile);
    product_version = DS.product_version;
    location = get_location_from_dataset(DS)
    print(f'product version = {product_version}')
    DS.close();

    RadarDS = pyart.io.read_cfradial(ncfile);

    dtime0 = cftime.num2pydate(RadarDS.time['data'][0],RadarDS.time['units']);
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S");
    nsweeps = RadarDS.nsweeps;
    
    vel_field = RadarDS.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10.66)
    vel_limit_upper = vel_field.get('field_limit_upper', 10.66)

    # Get valid sweep indices (optionally filtering out all-transition sweeps)
    valid_sweep_indices = get_valid_sweep_indices(RadarDS, skip_all_transition=skip_all_transition)
    
    if skip_all_transition and len(valid_sweep_indices) < nsweeps:
        print(f"Skipping {nsweeps - len(valid_sweep_indices)} sweep(s) that are 100% antenna transitions")

#    fig, ax = plt.subplots(nsweeps,4,figsize=(24,nsweeps*4),constrained_layout=True)
#    fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)

    display = pyart.graph.RadarDisplay(RadarDS);

    figpath = os.path.join(figpath,'rhi',datestr);
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    figpath_ccrest_az = os.path.join(figpath,ccrest_az);
    if not os.path.isdir(figpath_ccrest_az):
        os.makedirs(figpath_ccrest_az);


    if nsweeps>1:

        for s in valid_sweep_indices:
            print(f"sweep {s}/{nsweeps}");
            rhi_az = RadarDS.get_azimuth(s)[0];
            sweep_start_index = RadarDS.get_start(s);
            dtime_sweep = cftime.num2pydate(RadarDS.time['data'][sweep_start_index],RadarDS.time['units']);
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
            fig, ax = plt.subplots(2,2,figsize=(16,12), constrained_layout=False)
            fig.suptitle(
                f"{CAMPAIGN.upper()} Campaign | {location.upper()}\n"
                f"{NCAS_INSTRUMENT_NAME} Az {rhi_az:.1f} Deg. {dtime_sweep.strftime('%Y-%m-%dT%H:%M:%SZ')}",
                fontsize=14,
                y=0.98,
            )
            fig.subplots_adjust(left=0.07, right=0.97, top=0.87, bottom=0.08,
                                wspace=0.28, hspace=0.62)

            gatefilter = _build_ccrest_rhi_gatefilter(RadarDS, s, snr_threshold=-8.0)

            _plot_ccrest_rhi_fields(
                display,
                ax,
                s,
                gatefilter,
                vel_limit_lower,
                vel_limit_upper,
                hmax,
                xmin,
                xmax,
            )
            figname = f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime_sweep_str}_rhi_az{rhi_az:0.2f}_l1_{product_version}.png';
            plt.savefig(os.path.join(figpath_ccrest_az,figname),dpi=300);
            plt.close();

    else:
        fig, ax = plt.subplots(2,2,figsize=(16,12), constrained_layout=False)
        rhi_az = RadarDS.get_azimuth(0)[0]
        fig.suptitle(
            f"{CAMPAIGN.upper()} Campaign | {location.upper()}\n"
            f"{NCAS_INSTRUMENT_NAME} Az {rhi_az:.1f} Deg. {dtime0.strftime('%Y-%m-%dT%H:%M:%SZ')}",
            fontsize=14,
            y=0.98,
        )
        fig.subplots_adjust(left=0.07, right=0.97, top=0.87, bottom=0.08,
                            wspace=0.28, hspace=0.62)

        gatefilter = _build_ccrest_rhi_gatefilter(RadarDS, 0, snr_threshold=-8.0)

        _plot_ccrest_rhi_fields(
            display,
            ax,
            0,
            gatefilter,
            vel_limit_lower,
            vel_limit_upper,
            hmax,
            xmin,
            xmax,
        )
        figname = f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime0_str}_rhi_az{rhi_az:0.2f}_l1_{product_version}.png'
        plt.savefig(os.path.join(figpath_ccrest_az,figname),dpi=300);
        plt.close();






def _parse_geospatial_bounds(bounds_str):
    """Parse a geospatial_bounds string of the form
    'Bounding box: 50.95N -1.77E, 51.29N -1.44E'
    and return (lon_min, lon_max, lat_min, lat_max) or None on failure.
    """
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
    """Compute the outermost lon/lat bounds across all parseable PPI files for a day."""
    day_bounds = None

    print(f"Computing day-outermost bounds from {len(ppi_files)} PPI files")

    for ncfile in ppi_files:
        bounds_source = None
        try:
            with nc4.Dataset(ncfile) as ds:
                bounds = _parse_geospatial_bounds(getattr(ds, 'geospatial_bounds', ''))
                if bounds is not None:
                    bounds_source = 'geospatial_bounds'
                if bounds is None:
                    lon_min_attr = getattr(ds, 'geospatial_lon_min', None)
                    lon_max_attr = getattr(ds, 'geospatial_lon_max', None)
                    lat_min_attr = getattr(ds, 'geospatial_lat_min', None)
                    lat_max_attr = getattr(ds, 'geospatial_lat_max', None)
                    if None not in (lon_min_attr, lon_max_attr, lat_min_attr, lat_max_attr):
                        bounds = (
                            float(lon_min_attr),
                            float(lon_max_attr),
                            float(lat_min_attr),
                            float(lat_max_attr),
                        )
                        bounds_source = 'geospatial_lon/lat_min/max'
        except Exception as exc:
            print(f"Warning: unable to read bounds from {os.path.basename(ncfile)} ({exc})")
            continue

        if bounds is None:
            print(f"PPI bounds: {os.path.basename(ncfile)} -> unavailable")
            continue

        print(
            f"PPI bounds: {os.path.basename(ncfile)} ({bounds_source}) -> "
            f"lon [{bounds[0]:.4f}, {bounds[1]:.4f}], lat [{bounds[2]:.4f}, {bounds[3]:.4f}]"
        )

        if day_bounds is None:
            day_bounds = list(bounds)
        else:
            day_bounds[0] = min(day_bounds[0], bounds[0])
            day_bounds[1] = max(day_bounds[1], bounds[1])
            day_bounds[2] = min(day_bounds[2], bounds[2])
            day_bounds[3] = max(day_bounds[3], bounds[3])

        print(
            f"Running day union -> "
            f"lon [{day_bounds[0]:.4f}, {day_bounds[1]:.4f}], "
            f"lat [{day_bounds[2]:.4f}, {day_bounds[3]:.4f}]"
        )

    if day_bounds is None:
        return None

    return tuple(day_bounds)


def _plot_ccrest_ppi_fields(display, axes, sweep_idx, gatefilter, vel_min, vel_max, max_range):
    """Plot the standard 2x2 CCREST PPI field set."""
    from matplotlib import colors

    field_specs = [
        ("DBZ", "Equivalent Reflectivity Factor", axes[0, 0], dict(vmin=-40, vmax=40, cmap='HomeyerRainbow')),
        ("LDR", "Radar Linear Depolarization Ratio", axes[0, 1], dict(vmin=-35, vmax=5, cmap='SpectralExtended')),
        ("VEL", "Radial Velocity", axes[1, 0], dict(vmin=vel_min, vmax=vel_max, cmap='balance')),
        (
            "WIDTH",
            "Radar Doppler Spectrum Width",
            axes[1, 1],
            dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
                 cmap='SpectralExtended'),
        ),
    ]

    for field, short_title, ax, kwargs in field_specs:
        display.plot_ppi(
            field,
            ax=ax,
            sweep=sweep_idx,
            colorbar_flag=False,
            title_flag=False,
            gatefilter=gatefilter,
            filter_transitions=True,
            **kwargs,
        )
        _setup_ccrest_ppi_axes(ax, max_range)
        ax.set_title(short_title, fontsize=11, pad=8)
        mappable = display.plots[-1]
        field_meta = display.fields[field]
        label = f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})"
        cb = ax.get_figure().colorbar(
            mappable,
            ax=ax,
            orientation='vertical',
            shrink=0.9,
            pad=0.02,
        )
        cb.set_label(label)
        cb.ax.tick_params(labelsize=9)


def _setup_ccrest_ppi_axes(ax, max_range):
    """Apply consistent axis styling to CCREST PPI subplots."""
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_xlabel('East-West distance from radar (km)')
    ax.set_ylabel('North-South distance from radar (km)')
    ax.grid(True)
    ax.set_aspect('equal')


def make_ccrest_ppi_plot(ncfile, figpath, blflag=False, skip_all_transition=False):
    """Create standard CCREST PPI quicklook plots (cartesian x/y)."""
    max_range = 20.0 if blflag else 30.0

    DS = nc4.Dataset(ncfile)
    product_version = DS.product_version
    location = get_location_from_dataset(DS)
    DS.close()

    radar = pyart.io.read_cfradial(ncfile)
    display = pyart.graph.RadarDisplay(radar)

    vel_field = radar.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10.66)
    vel_limit_upper = vel_field.get('field_limit_upper', 10.66)

    outdir = os.path.join(figpath, 'ppi', datestr)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    valid_sweep_indices = get_valid_sweep_indices(radar, skip_all_transition=skip_all_transition)

    for s in valid_sweep_indices:
        sweep_start_index = radar.get_start(s)
        dtime_sweep = cftime.num2pydate(radar.time['data'][sweep_start_index], radar.time['units'])
        dtime_sweep_str = dtime_sweep.strftime('%Y%m%d-%H%M%S')
        ppi_el = float(np.nanmean(radar.get_elevation(s)))

        gatefilter = pyart.correct.GateFilter(radar)
        gatefilter.exclude_below('SNR', -20)

        fig, ax = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=False)
        fig.suptitle(
            f"{CAMPAIGN.upper()} Campaign | {location.upper()}\n"
            f"{NCAS_INSTRUMENT_NAME} PPI El {ppi_el:.1f} Deg. {dtime_sweep.strftime('%Y-%m-%dT%H:%M:%SZ')}",
            fontsize=14,
            y=0.98,
        )
        fig.subplots_adjust(left=0.07, right=0.97, top=0.87, bottom=0.08,
                            wspace=0.28, hspace=0.62)

        _plot_ccrest_ppi_fields(
            display,
            ax,
            s,
            gatefilter,
            vel_limit_lower,
            vel_limit_upper,
            max_range,
        )

        figname = (
            f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime_sweep_str}_'
            f'ppi_el{ppi_el:0.2f}_l1_{product_version}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=300)
        plt.close()


def make_ccrest_ppi_map_plot(ncfile, figpath, blflag=False, skip_all_transition=False,
                             day_extent=None, colorbar_shrink=0.86):
    """Create georeferenced CCREST PPI map quicklook plots."""
    try:
        import cartopy.crs as ccrs
    except ImportError:
        print('Skipping PPI map quicklook: cartopy is not available')
        return

    max_range = 20.0 if blflag else 30.0

    DS = nc4.Dataset(ncfile)
    product_version = DS.product_version
    location = get_location_from_dataset(DS)
    # Parse geospatial bounds from the NetCDF geospatial_bounds string attribute
    _bounds = _parse_geospatial_bounds(getattr(DS, 'geospatial_bounds', ''))
    DS.close()

    radar = pyart.io.read_cfradial(ncfile)
    display = pyart.graph.RadarMapDisplay(radar)

    from matplotlib import colors

    vel_field = radar.fields['VEL']
    vel_limit_lower = vel_field.get('field_limit_lower', -10.66)
    vel_limit_upper = vel_field.get('field_limit_upper', 10.66)

    field_specs = [
        ('DBZ', 'Equivalent Reflectivity Factor', dict(vmin=-40, vmax=40, cmap='HomeyerRainbow')),
        ('LDR', 'Radar Linear Depolarization Ratio', dict(vmin=-35, vmax=5, cmap='SpectralExtended')),
        ('VEL', 'Radial Velocity', dict(vmin=vel_limit_lower, vmax=vel_limit_upper, cmap='balance')),
        (
            'WIDTH',
            'Radar Doppler Spectrum Width',
            dict(norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)), cmap='SpectralExtended'),
        ),
    ]

    outdir = os.path.join(figpath, 'ppi_map', datestr)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    valid_sweep_indices = get_valid_sweep_indices(radar, skip_all_transition=skip_all_transition)

    for s in valid_sweep_indices:
        sweep_start_index = radar.get_start(s)
        dtime_sweep = cftime.num2pydate(radar.time['data'][sweep_start_index], radar.time['units'])
        dtime_sweep_str = dtime_sweep.strftime('%Y%m%d-%H%M%S')
        ppi_el = float(np.nanmean(radar.get_elevation(s)))

        # Define bounding box from day-wide extent (optional), then file metadata,
        # then fall back to computed margins.
        radar_lon = radar.longitude['data'][0]
        radar_lat = radar.latitude['data'][0]
        if day_extent is not None:
            lon_min, lon_max, lat_min, lat_max = day_extent
            extent_source = 'day-outermost'
        elif _bounds is not None:
            lon_min, lon_max, lat_min, lat_max = _bounds
            extent_source = 'file'
        else:
            lon_min, lon_max = radar_lon - 0.50, radar_lon + 0.50
            lat_min, lat_max = radar_lat - 0.35, radar_lat + 0.35
            extent_source = 'fallback'
        print(
            f"PPI map extent ({extent_source}): "
            f"lon [{lon_min:.4f}, {lon_max:.4f}], lat [{lat_min:.4f}, {lat_max:.4f}]"
        )

        fig, axes = plt.subplots(
            2, 2, figsize=(16, 16),
            subplot_kw=dict(projection=ccrs.OSGB()),
            constrained_layout=False,
        )
        fig.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.05, wspace=0.25, hspace=0.25)
        fig.suptitle(
            f"{CAMPAIGN.upper()} Campaign | {location.upper()}\n"
            f"{NCAS_INSTRUMENT_NAME} PPI Map | El {ppi_el:.1f} Deg. "
            f"{dtime_sweep.strftime('%Y-%m-%dT%H:%M:%SZ')}",
            fontsize=12,
            y=0.97,
        )

        for (field, short_title, kwargs), panel_ax in zip(field_specs, axes.flatten()):
            display.plot_ppi_map(
                field,
                s,
                projection=ccrs.OSGB(),
                ax=panel_ax,
                lon_lines=None,
                lat_lines=None,
                colorbar_flag=False,
                title_flag=False,
                filter_transitions=True,
                gatefilter=None,
                **kwargs,
            )

            # Set extent to the bounding box before adding basemap tiles so that
            # contextily fetches tiles for the correct area.
            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

            # Add base map tiles in OSGB (EPSG:27700) with OSM fallback.
            # reset_extent=False prevents contextily from overriding the explicit extent above.
            try:
                ctx.add_basemap(
                    panel_ax, zoom=11, source=ctx.providers.OpenTopoMap,
                    crs='EPSG:27700', reset_extent=False,
                )
            except Exception as e_topo:
                print(f"Warning: OpenTopoMap basemap unavailable ({e_topo}); trying OSM Mapnik")
                try:
                    ctx.add_basemap(
                        panel_ax, zoom=11, source=ctx.providers.OpenStreetMap.Mapnik,
                        crs='EPSG:27700', reset_extent=False,
                    )
                except Exception as e_osm:
                    print(f"Warning: Could not add basemap tiles: {e_osm}")
            panel_ax.set_title(short_title, fontsize=10, pad=6)

            display.plot_range_rings([10, 20, int(max_range)], ax=panel_ax, col='0.5', lw=0.8)

            # Re-enforce extent after range rings: pyart's plot_range_rings can draw
            # circles that extend beyond the specified bounds, causing Cartopy to
            # auto-expand the view.
            panel_ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

            gl = panel_ax.gridlines(
                draw_labels=True,
                linewidth=0.4,
                color='gray',
                alpha=0.6,
                x_inline=False,
                y_inline=False,
            )
            gl.right_labels = False
            gl.top_labels = False
            gl.left_labels = True
            gl.bottom_labels = True

            mappable = display.plots[-1]
            field_meta = display.fields[field]
            label = f"{field_meta.get('long_name', field)} ({field_meta.get('units', '')})"
            cb = fig.colorbar(
                mappable,
                ax=panel_ax,
                orientation='vertical',
                shrink=colorbar_shrink,
                pad=0.02,
            )
            cb.set_label(label)
            cb.ax.tick_params(labelsize=8)

        figname = (
            f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime_sweep_str}_'
            f'ppi-map_el{ppi_el:0.2f}_l1_{product_version}.png'
        )
        plt.savefig(os.path.join(outdir, figname), dpi=300)
        plt.close()


def make_ccrest_ppi_plots_day(datestr, inpath, figpath, blflag=False,
                              skip_all_transition=False, use_day_outermost_extent=False,
                              ppi_map_colorbar_shrink=0.86):
    """Generate standard and map PPI quicklooks for all PPI files on a given day."""
    ppi_files = collect_ccrest_ppi_files(datestr, inpath)
    print(f'ppi files = {ppi_files}')

    day_extent = None
    if use_day_outermost_extent:
        day_extent = _compute_day_outermost_ppi_bounds(ppi_files)
        if day_extent is not None:
            print(
                f"Using day-outermost PPI-map extent: "
                f"lon [{day_extent[0]:.4f}, {day_extent[1]:.4f}], "
                f"lat [{day_extent[2]:.4f}, {day_extent[3]:.4f}]"
            )
        else:
            print("Warning: no parseable day bounds found; falling back to per-file extents")

    for ncfile in ppi_files:
        make_ccrest_ppi_plot(
            ncfile,
            figpath,
            blflag=blflag,
            skip_all_transition=skip_all_transition,
        )
        make_ccrest_ppi_map_plot(
            ncfile,
            figpath,
            blflag=blflag,
            skip_all_transition=skip_all_transition,
            day_extent=day_extent,
            colorbar_shrink=ppi_map_colorbar_shrink,
        )


def make_ccrest_vpt_plot_day(datestr,inpath,figpath,blflag=False):
    
    if blflag:
        hmax = 4;
    else:
        hmax = 12;

    dbz_cmap = 'HomeyerRainbow';
    vel_cmap = 'balance';
    ldr_cmap = 'SpectralExtended';
    #ldr_cmap = 'pyart_ChaseSpectral';
    spw_cmap = 'SpectralExtended';
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
    vpt_files = sorted([os.path.join(inpath_date,f) for f in glob.glob('*{}*vpt*.nc'.format(datestr))])
    if not vpt_files:
        print(f'No VPT files found for {datestr}')
        return

    DS = nc4.Dataset(vpt_files[0]);
    product_version = DS.product_version;
    location = get_location_from_dataset(DS)
    print(f'product version = {product_version}')
    DS.close();

    # Load all VPT files for the day as a list of radar objects
    vpt_radars = [pyart.io.read_cfradial(f) for f in vpt_files]
    RadarDS_VPT = vpt_radars[0]
    nsweeps = sum(r.nsweeps for r in vpt_radars)
    

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
    fig.set_constrained_layout_pads(h_pad=0.16, w_pad=0.02, hspace=0.04, wspace=0.02)

    
    #pyart.correct.despeckle_field(RadarDS, "SNR", size=3, threshold=-20, gatefilter=gatefilter, delta=5.0)

    #display = pyart.graph.RadarDisplay(RadarDS_VPT);
    RadarSweepDS = RadarDS_VPT.extract_sweeps([0])
    print("sweep = 0");
    gatefilter = pyart.correct.GateFilter(RadarSweepDS)
    gatefilter.exclude_below('SNR', -20)
    if RadarSweepDS.antenna_transition is not None:
        gatefilter.gate_excluded[RadarSweepDS.antenna_transition['data'] == 1, :] = True
    import cftime
    dtime = cftime.num2pydate(RadarDS_VPT.time['data'],RadarDS_VPT.time['units'])

    dt_min = dtime[0].replace(hour=0,minute=0,second=0); #datetime.strptime(datestr, '%Y%m%d')
    dt_max = dt_min + datetime.timedelta(days=1)

    time_str = dtime[0].strftime("%Y-%m-%d");
    fig.suptitle(
        f"{CAMPAIGN.upper()} campaign | {location.upper()} | {NCAS_INSTRUMENT_NAME}\n{time_str}",
        fontsize=12,
    )
    
   
    _plot_vpt_field_manual(
        RadarSweepDS,
        "DBZ",
        ax[0],
        gatefilter=gatefilter,
        vmin=-40,
        vmax=40,
        norm=None,
        cmap=dbz_cmap,
        colorbar_flag=True,
        colorbar_orient='horizontal',
    )
    ax[0].set_ylim(0,hmax)
    ax[0].grid(True);
    _plot_vpt_field_manual(
        RadarSweepDS,
        "VEL",
        ax[1],
        gatefilter=gatefilter,
        vmin=velmin,
        vmax=velmax,
        norm=None,
        cmap=vel_cmap,
        colorbar_flag=True,
        colorbar_orient='horizontal',
    )
    ax[1].set_ylim(0,hmax)
    ax[1].grid(True)
    _plot_vpt_field_manual(
        RadarSweepDS,
        "WIDTH",
        ax[2],
        gatefilter=gatefilter,
        norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)),
        cmap=spw_cmap,
        colorbar_flag=True,
        colorbar_orient='horizontal',
    )
    ax[2].set_ylim(0,hmax)
    ax[2].grid(True)
    _plot_vpt_field_manual(
        RadarSweepDS,
        "LDR",
        ax[3],
        gatefilter=gatefilter,
        vmin=-35,
        vmax=5,
        norm=None,
        cmap=ldr_cmap,
        colorbar_flag=True,
        colorbar_orient='horizontal',
    )
    ax[3].set_ylim(0,hmax)
    ax[3].grid(True)

    ax[0].set_xlim(dt_min,dt_max);
    ax[1].set_xlim(dt_min,dt_max);
    ax[2].set_xlim(dt_min,dt_max);
    ax[3].set_xlim(dt_min,dt_max);

    time_locator = mdates.HourLocator(interval=3)
    time_fmt = mdates.DateFormatter('%H:%M')
    for axis in ax:
        axis.xaxis.set_major_locator(time_locator)
        axis.xaxis.set_major_formatter(time_fmt)
        axis.set_ylabel('Height above radar (km)')

    ax[0].set_xlabel('Time (UTC)');
    ax[1].set_xlabel('Time (UTC)');
    ax[2].set_xlabel('Time (UTC)');
    ax[3].set_xlabel('Time (UTC)');

    # Plot remaining sweeps from all VPT files for the day
    for file_idx, RadarDS in enumerate(vpt_radars):
        start_s = 1 if file_idx == 0 else 0
        for s in range(start_s, RadarDS.nsweeps):
            RadarSweepDS = RadarDS.extract_sweeps([s])
            print(f"file = {file_idx}, sweep = {s}");
            gatefilter = pyart.correct.GateFilter(RadarSweepDS)
            gatefilter.exclude_below('SNR', -20)
            if RadarSweepDS.antenna_transition is not None:
                gatefilter.gate_excluded[RadarSweepDS.antenna_transition['data'] == 1, :] = True
            _plot_vpt_field_manual(
                RadarSweepDS,
                "DBZ",
                ax[0],
                gatefilter=gatefilter,
                vmin=-40,
                vmax=40,
                norm=None,
                cmap=dbz_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "VEL",
                ax[1],
                gatefilter=gatefilter,
                vmin=velmin,
                vmax=velmax,
                norm=None,
                cmap=vel_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "WIDTH",
                ax[2],
                gatefilter=gatefilter,
                norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)),
                cmap=spw_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "LDR",
                ax[3],
                gatefilter=gatefilter,
                vmin=-35,
                vmax=5,
                norm=None,
                cmap=ldr_cmap,
                colorbar_flag=False,
            )

    if (nsweeps_prev>0):
        for s in range(nsweeps_prev-1,nsweeps_prev):
            RadarSweepDS = RadarDS_VPT_prev.extract_sweeps([s])
            print(f"sweep = {s}");
            gatefilter = pyart.correct.GateFilter(RadarSweepDS)
            gatefilter.exclude_below('SNR', -20)
            if RadarSweepDS.antenna_transition is not None:
                gatefilter.gate_excluded[RadarSweepDS.antenna_transition['data'] == 1, :] = True
            _plot_vpt_field_manual(
                RadarSweepDS,
                "DBZ",
                ax[0],
                gatefilter=gatefilter,
                vmin=-40,
                vmax=40,
                norm=None,
                cmap=dbz_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "VEL",
                ax[1],
                gatefilter=gatefilter,
                vmin=velmin,
                vmax=velmax,
                norm=None,
                cmap=vel_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "WIDTH",
                ax[2],
                gatefilter=gatefilter,
                norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)),
                cmap=spw_cmap,
                colorbar_flag=False,
            )
            _plot_vpt_field_manual(
                RadarSweepDS,
                "LDR",
                ax[3],
                gatefilter=gatefilter,
                vmin=-35,
                vmax=5,
                norm=None,
                cmap=ldr_cmap,
                colorbar_flag=False,
            )
     
    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)
    ax[3].grid(True)

    ax[0].set_xlabel('Time (UTC)');
    ax[1].set_xlabel('Time (UTC)');
    ax[2].set_xlabel('Time (UTC)');
    ax[3].set_xlabel('Time (UTC)');

    figname = f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{datestr}_vpt_l1_{product_version}.png';
    
    if blflag:
        figname = figname.replace('.png','_bl.png');

    figpath = os.path.join(figpath,'vpt');
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    plt.savefig(os.path.join(figpath,figname),dpi=300);

    plt.close();


def make_ccrest_pointing_plot(ncfile, figpath, blflag=False):
    """Make quicklook plots for pointing sweeps (time vs range)."""
    if blflag:
        hmax = 4
    else:
        hmax = 12

    dbz_cmap = 'HomeyerRainbow'
    vel_cmap = 'balance'
    ldr_cmap = 'SpectralExtended'
    spw_cmap = 'SpectralExtended'

    from matplotlib import colors

    DS = nc4.Dataset(ncfile)
    product_version = DS.product_version
    location = get_location_from_dataset(DS)
    DS.close()

    radar = pyart.io.read_cfradial(ncfile)
    vel_field = radar.fields['VEL']
    velmin = vel_field.get('field_limit_lower', -5.0)
    velmax = vel_field.get('field_limit_upper', 5.0)

    outdir = os.path.join(figpath, 'pointing', datestr)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for s in range(radar.nsweeps):
        radar_sweep = radar.extract_sweeps([s])

        gatefilter = pyart.correct.GateFilter(radar_sweep)
        gatefilter.exclude_below('SNR', -20)
        if radar_sweep.antenna_transition is not None:
            gatefilter.gate_excluded[radar_sweep.antenna_transition['data'] == 1, :] = True

        dtime = cftime.num2pydate(radar_sweep.time['data'], radar_sweep.time['units'])
        dt_min = dtime[0]
        dt_max = dtime[-1]
        if dt_max <= dt_min:
            dt_max = dt_min + datetime.timedelta(minutes=1)

        az = float(np.nanmean(radar_sweep.azimuth['data']))
        el = float(np.nanmean(radar_sweep.elevation['data']))

        fig, ax = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=True)
        fig.set_constrained_layout_pads(h_pad=0.16, w_pad=0.02, hspace=0.04, wspace=0.02)
        fig.suptitle(
            f"{CAMPAIGN.upper()} | {location.upper()} | {NCAS_INSTRUMENT_NAME}\n"
            f"Pointing Az {az:.1f} deg El {el:.1f} deg | {dtime[0].strftime('%Y-%m-%d')} UTC",
            fontsize=12,
        )

        _plot_vpt_field_manual(
            radar_sweep,
            'DBZ',
            ax[0],
            gatefilter=gatefilter,
            vmin=-40,
            vmax=40,
            norm=None,
            cmap=dbz_cmap,
            colorbar_flag=True,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep,
            'VEL',
            ax[1],
            gatefilter=gatefilter,
            vmin=velmin,
            vmax=velmax,
            norm=None,
            cmap=vel_cmap,
            colorbar_flag=True,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep,
            'WIDTH',
            ax[2],
            gatefilter=gatefilter,
            norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
            cmap=spw_cmap,
            colorbar_flag=True,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep,
            'LDR',
            ax[3],
            gatefilter=gatefilter,
            vmin=-35,
            vmax=5,
            norm=None,
            cmap=ldr_cmap,
            colorbar_flag=True,
            colorbar_orient='horizontal',
        )

        locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
        formatter = mdates.DateFormatter('%H:%M')
        for axis in ax:
            axis.set_xlim(dt_min, dt_max)
            axis.set_ylim(0, hmax)
            axis.grid(True)
            axis.set_ylabel('Height above radar (km)')
            axis.set_xlabel('Time (UTC)')
            axis.xaxis.set_major_locator(locator)
            axis.xaxis.set_major_formatter(formatter)

        dtime0_str = dtime[0].strftime('%Y%m%d-%H%M%S')
        figname = (
            f"{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime0_str}_pointing_"
            f"s{s:02d}_l1_{product_version}.png"
        )
        plt.savefig(os.path.join(outdir, figname), dpi=300)
        plt.close()


def make_ccrest_pointing_multisweep_plot(ncfile, figpath, blflag=False):
    """Make a single quicklook displaying all sweeps of a multi-sweep pointing file.

    All sweeps are plotted on shared time axes (one column of 4 panels: DBZ, VEL,
    WIDTH, LDR), colorbars drawn once from the first sweep.
    """
    if blflag:
        hmax = 4
    else:
        hmax = 12

    dbz_cmap = 'HomeyerRainbow'
    vel_cmap = 'balance'
    ldr_cmap = 'SpectralExtended'
    spw_cmap = 'SpectralExtended'

    from matplotlib import colors

    DS = nc4.Dataset(ncfile)
    product_version = DS.product_version
    location = get_location_from_dataset(DS)
    DS.close()

    radar = pyart.io.read_cfradial(ncfile)
    vel_field = radar.fields['VEL']
    velmin = vel_field.get('field_limit_lower', -5.0)
    velmax = vel_field.get('field_limit_upper', 5.0)

    outdir = os.path.join(figpath, 'pointing', datestr)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Time bounds for xlim: use full extent across all sweeps
    dtime_all = cftime.num2pydate(radar.time['data'], radar.time['units'])
    dt_min = dtime_all[0]
    dt_max = dtime_all[-1]
    if dt_max <= dt_min:
        dt_max = dt_min + datetime.timedelta(minutes=1)

    az_mean = float(np.nanmean(radar.azimuth['data']))
    el_mean = float(np.nanmean(radar.elevation['data']))
    time_str = dtime_all[0].strftime('%Y-%m-%d')

    fig, ax = plt.subplots(4, 1, figsize=(12, 18), constrained_layout=True)
    fig.set_constrained_layout_pads(h_pad=0.16, w_pad=0.02, hspace=0.04, wspace=0.02)
    fig.suptitle(
        f"{CAMPAIGN.upper()} | {location.upper()} | {NCAS_INSTRUMENT_NAME}\n"
        f"Pointing Az {az_mean:.1f}\u00b0 El {el_mean:.1f}\u00b0 | {radar.nsweeps} sweeps | {time_str} UTC",
        fontsize=12,
    )

    for s in range(radar.nsweeps):
        radar_sweep = radar.extract_sweeps([s])

        gatefilter = pyart.correct.GateFilter(radar_sweep)
        gatefilter.exclude_below('SNR', -20)
        if radar_sweep.antenna_transition is not None:
            gatefilter.gate_excluded[radar_sweep.antenna_transition['data'] == 1, :] = True

        colorbar_flag = (s == 0)

        _plot_vpt_field_manual(
            radar_sweep, 'DBZ', ax[0],
            gatefilter=gatefilter, vmin=-40, vmax=40, norm=None,
            cmap=dbz_cmap, colorbar_flag=colorbar_flag,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep, 'VEL', ax[1],
            gatefilter=gatefilter, vmin=velmin, vmax=velmax, norm=None,
            cmap=vel_cmap, colorbar_flag=colorbar_flag,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep, 'WIDTH', ax[2],
            gatefilter=gatefilter,
            norm=colors.LogNorm(vmin=1e-1 * np.sqrt(1e-1), vmax=np.sqrt(1e1)),
            cmap=spw_cmap, colorbar_flag=colorbar_flag,
            colorbar_orient='horizontal',
        )
        _plot_vpt_field_manual(
            radar_sweep, 'LDR', ax[3],
            gatefilter=gatefilter, vmin=-35, vmax=5, norm=None,
            cmap=ldr_cmap, colorbar_flag=colorbar_flag,
            colorbar_orient='horizontal',
        )

    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    formatter = mdates.DateFormatter('%H:%M:%S')
    for axis in ax:
        axis.set_xlim(dt_min, dt_max)
        axis.set_ylim(0, hmax)
        axis.grid(True)
        axis.set_ylabel('Height above radar (km)')
        axis.set_xlabel('Time (UTC)')
        axis.xaxis.set_major_locator(locator)
        axis.xaxis.set_major_formatter(formatter)

    dtime0_str = dtime_all[0].strftime('%Y%m%d-%H%M%S')
    figname = (
        f"{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{dtime0_str}_pointing_"
        f"l1_{product_version}.png"
    )
    plt.savefig(os.path.join(outdir, figname), dpi=300)
    plt.close()


def make_ccrest_pointing_plots_day(datestr, inpath, figpath, blflag=False):
    """Generate quicklooks for all pointing files on a given day."""
    pointing_files = collect_ccrest_pointing_files(datestr, inpath)
    print(f'pointing files = {pointing_files}')
    for ncfile in pointing_files:
        radar = pyart.io.read_cfradial(ncfile)
        if radar.nsweeps > 1:
            make_ccrest_pointing_multisweep_plot(ncfile, figpath, blflag=blflag)
        else:
            make_ccrest_pointing_plot(ncfile, figpath, blflag=blflag)

def make_ccrest_vad_plot_day(datestr,inpath,figpath,zlevels,blflag=False):

    if blflag:
        hmax = 4;
    else:
        hmax = 12;

    inpath_date = os.path.join(inpath,datestr);

    current_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    prev_date = current_date - datetime.timedelta(days=1)
    prevstr = prev_date.strftime('%Y%m%d')

    # Search current day's directory for VAD files named with datestr
    vad_files = sorted(glob.glob(os.path.join(inpath_date, '*{}*vad*.nc'.format(datestr))))

    # Also search previous day's directory — a sweep starting just after midnight
    # may have been written there during campaign processing
    inpath_prev = os.path.join(inpath, prevstr)
    if os.path.isdir(inpath_prev):
        vad_files += sorted(glob.glob(os.path.join(inpath_prev, '*{}*vad*.nc'.format(datestr))))

    vad_files = sorted(set(vad_files))

    if not vad_files:
        print(f'No VAD files found for {datestr}')
        return
    vad_file = vad_files[0];

    DS = nc4.Dataset(vad_file);
    product_version = DS.product_version;
    location = get_location_from_dataset(DS)
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
    fig.suptitle(f"{CAMPAIGN.upper()} | {location.upper()}")
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

    figname = f'{NCAS_INSTRUMENT_NAME}_{location}_ccrest-m_{datestr}_vad_l1_{product_version}.png'
    plt.savefig(os.path.join(figpath, figname),dpi=300)

    plt.close();



def make_ccrest1_rhi_plots_day(datestr,inpath,figpath,grouped_rhi_files,blflag=False,skip_all_transition=False):
    ccrest1_rhi_files = grouped_rhi_files['ccrest-1']

    print(f'ccrest-1 files = ',ccrest1_rhi_files);

    for f in ccrest1_rhi_files:
        make_ccrest_rhi_plot(f,figpath,'ccrest-1',blflag=blflag,skip_all_transition=skip_all_transition);

    return

def make_ccrest2_rhi_plots_day(datestr,inpath,figpath,grouped_rhi_files,blflag=False,skip_all_transition=False):
    ccrest2_rhi_files = grouped_rhi_files['ccrest-2']

    print(f'ccrest-2 files = ',ccrest2_rhi_files);

    for f in ccrest2_rhi_files:
        make_ccrest_rhi_plot(f,figpath,'ccrest-2',blflag=blflag,skip_all_transition=skip_all_transition);

    return

def main():
    """Main function."""
    global datestr

    (
        datestr,
        inpath_arg,
        outpath_arg,
        blflag,
        skip_all_transition,
        vpt_only,
        ppi_only,
        pointing_only,
        ppi_map_day_extent,
        ppi_map_colorbar_shrink,
    ) = parse_command_line()

    inpath, figpath = setup_paths(datestr)

    if inpath_arg:
        inpath = inpath_arg
    if outpath_arg:
        figpath = outpath_arg

    print(f"Processing date: {datestr}")
    print(f"Input path: {inpath}")
    print(f"Output path: {figpath}")
    print(f"Boundary layer mode: {blflag}")
    print(f"Skip all-transition sweeps: {'enabled' if skip_all_transition else 'disabled'}")
    print(f"PPI-map day-outermost extent: {'enabled' if ppi_map_day_extent else 'disabled'}")
    print(f"PPI-map colorbar shrink: {ppi_map_colorbar_shrink:.2f}")
    if vpt_only:
        plot_mode = 'VPT only'
    elif ppi_only:
        plot_mode = 'PPI + PPI-MAP only'
    elif pointing_only:
        plot_mode = 'Pointing only'
    else:
        plot_mode = 'VPT + PPI + PPI-MAP + RHI + VAD'
    print(f"Plot mode: {plot_mode}")

    inpath_date = os.path.join(inpath, datestr)
    if not os.path.isdir(inpath_date):
        print(f"No input directory for date {datestr}: {inpath_date}")
        print("Skipping this date gracefully.")
        return

    if vpt_only:
        make_ccrest_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
        return

    if pointing_only:
        make_ccrest_pointing_plots_day(datestr, inpath, figpath, blflag=blflag)
        return

    if ppi_only:
        make_ccrest_ppi_plots_day(
            datestr,
            inpath,
            figpath,
            blflag=blflag,
            skip_all_transition=skip_all_transition,
            use_day_outermost_extent=ppi_map_day_extent,
            ppi_map_colorbar_shrink=ppi_map_colorbar_shrink,
        )
        return

    make_ccrest_vpt_plot_day(datestr, inpath, figpath, blflag=blflag)
    make_ccrest_ppi_plots_day(
        datestr,
        inpath,
        figpath,
        blflag=blflag,
        skip_all_transition=skip_all_transition,
        use_day_outermost_extent=ppi_map_day_extent,
        ppi_map_colorbar_shrink=ppi_map_colorbar_shrink,
    )
    grouped_rhi_files = collect_ccrest_rhi_files(datestr, inpath)
    make_ccrest1_rhi_plots_day(
        datestr,
        inpath,
        figpath,
        grouped_rhi_files,
        blflag=blflag,
        skip_all_transition=skip_all_transition,
    )
    make_ccrest2_rhi_plots_day(
        datestr,
        inpath,
        figpath,
        grouped_rhi_files,
        blflag=blflag,
        skip_all_transition=skip_all_transition,
    )
    make_ccrest_pointing_plots_day(datestr, inpath, figpath, blflag=blflag)
    zlevels = np.arange(100, 15000, 100)
    make_ccrest_vad_plot_day(datestr, inpath, figpath, zlevels)


if __name__ == '__main__':
    main()

