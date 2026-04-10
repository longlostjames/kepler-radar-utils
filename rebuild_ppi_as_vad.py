#!/usr/bin/env python3
"""rebuild_ppi_as_vad.py

Rebuild PPI scan files that have sentinel elevation values as a single
multisweep VAD file, recovering elevation values from the axis log.

Each input PPI file becomes one sweep in the output.  The sweep_mode is
set to "vad", fixed_angle to the median elevation recovered from the log,
and scan_rate / target_scan_rate variables are dropped so the output matches
the format produced by the normal CCREST-M processing pipeline.

Usage
-----
    python rebuild_ppi_as_vad.py YYYYMMDD [L1_PATH] [options]

Examples
--------
    python rebuild_ppi_as_vad.py 20240307
    python rebuild_ppi_as_vad.py 20240307 /path/to/L1 \\
        --logp-path /path/to/logp --output-dir /tmp/out
"""

import argparse
import datetime
import getpass
import glob
import gzip
import os
import re
import socket
import sys

import netCDF4 as nc
import numpy as np

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
L1_BASE = (
    "/gws/pw/j07/ncas_obs_vol2/cao/processing/"
    "ncas-mobile-ka-band-radar-1/ccrest-m/L1_v1.0.1"
)
LOGP_BASE = (
    "/gws/pw/j07/ncas_obs_vol2/cao/raw_data/"
    "ncas-mobile-ka-band-radar-1/data/campaign/logp"
)
SENTINEL_THRESHOLD = -900.0

# Moment variable names (time, range)
MOMENT_VARS = ["DBZ", "VEL", "WIDTH", "LDR", "SNR", "RHOHX", "PHIHX"]
# Per-ray 1-D variables (besides time, azimuth, elevation, antenna_transition)
SCALAR_RAY_VARS = ["prt", "pulse_width", "radar_measured_transmit_power_h"]
# Variables that exist in PPI files but should NOT be in a VAD file
DROP_VARS = {"scan_rate", "target_scan_rate"}


# ---------------------------------------------------------------------------
# Axis log helpers
# ---------------------------------------------------------------------------

def _load_axis_log(logp_dir, datestr):
    """Return (unix_times, az_deg, el_deg) arrays for *datestr* from the log.

    Load log files whose filename falls within ±1 h of the day, but keep
    only samples inside the calendar day (UTC).
    """
    dt_day = datetime.datetime.strptime(datestr, "%Y%m%d").replace(
        tzinfo=datetime.timezone.utc
    )
    day_start_unix = dt_day.timestamp()
    day_end_unix = (dt_day + datetime.timedelta(days=1)).timestamp()

    scan_start = dt_day - datetime.timedelta(hours=1)
    scan_end = dt_day + datetime.timedelta(days=1, hours=1)

    times_list, az_list, el_list = [], [], []

    for lf in sorted(glob.glob(os.path.join(logp_dir, "*.axis.gz"))):
        m = re.match(r"(\d{12})\.axis\.gz", os.path.basename(lf))
        if not m:
            continue
        fdt = datetime.datetime.strptime(m.group(1), "%Y%m%d%H%M").replace(
            tzinfo=datetime.timezone.utc
        )
        if fdt < scan_start or fdt > scan_end:
            continue
        try:
            with gzip.open(lf, "rt") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        t = float(parts[0])
                        a = float(parts[1])
                        e = float(parts[2])
                    except ValueError:
                        continue
                    if day_start_unix <= t < day_end_unix:
                        times_list.append(t)
                        az_list.append(a)
                        el_list.append(e)
        except Exception as exc:
            print(f"  WARNING: could not read {lf}: {exc}")

    if not times_list:
        return None, None, None

    times = np.array(times_list)
    az = np.array(az_list)
    el = np.array(el_list)
    idx = np.argsort(times)
    return times[idx], az[idx], el[idx]


def _ray_unix_times(nc_ds):
    """Return ray times as Unix timestamps from an open netCDF4 Dataset."""
    time_var = nc_ds.variables["time"]
    units = time_var.units
    m = re.match(
        r"seconds since (\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", units
    )
    if not m:
        raise ValueError(f"Cannot parse time units: {units!r}")
    ref = datetime.datetime.strptime(
        m.group(1), "%Y-%m-%dT%H:%M:%SZ"
    ).replace(tzinfo=datetime.timezone.utc)
    return ref.timestamp() + np.asarray(time_var[:])


def _interp_log_el(ray_unix, log_times, log_el):
    """Nearest-neighbour interpolation of log elevation at ray Unix times."""
    idx = np.searchsorted(log_times, ray_unix)
    idx_next = np.clip(idx, 0, len(log_times) - 1)
    idx_prev = np.clip(idx - 1, 0, len(log_times) - 1)
    diff_next = np.abs(log_times[idx_next] - ray_unix)
    diff_prev = np.abs(log_times[idx_prev] - ray_unix)
    use_prev = diff_prev < diff_next
    return np.where(use_prev, log_el[idx_prev], log_el[idx_next]).astype(
        np.float32
    )


# ---------------------------------------------------------------------------
# Main rebuild routine
# ---------------------------------------------------------------------------

def rebuild_vad(
    datestr,
    l1_path,
    logp_path,
    output_dir=None,
    sentinel_threshold=SENTINEL_THRESHOLD,
):
    day_dir = os.path.join(l1_path, datestr)
    ppi_files = sorted(
        glob.glob(os.path.join(day_dir, f"*_{datestr}*_ppi_l1_*.nc"))
    )
    if not ppi_files:
        sys.exit(f"No PPI files found for {datestr} in {day_dir}")

    print(f"Found {len(ppi_files)} PPI files")

    # ------------------------------------------------------------------
    # Load axis log
    # ------------------------------------------------------------------
    print("Loading axis log...")
    log_t, log_az, log_el = _load_axis_log(logp_path, datestr)
    if log_t is None:
        print(
            "  WARNING: no axis log data found — "
            "using fixed elevation of 82.0°"
        )
        use_log = False
    else:
        print(f"  Loaded {len(log_t):,} log samples")
        use_log = True

    # ------------------------------------------------------------------
    # Open first file to build template info
    # ------------------------------------------------------------------
    with nc.Dataset(ppi_files[0]) as f0:
        range_data = f0.variables["range"][:].data.copy()
        range_attrs = {
            k: f0.variables["range"].getncattr(k)
            for k in f0.variables["range"].ncattrs()
        }
        latitude = float(f0.variables["latitude"][:])
        longitude = float(f0.variables["longitude"][:])
        altitude = float(f0.variables["altitude"][:])
        frequency = float(f0.variables["frequency"][:])
        az_correction = float(f0.variables["azimuth_correction"][:])

        global_attrs = {k: f0.getncattr(k) for k in f0.ncattrs()}

        # Per-variable attributes and dtypes
        all_var_names = (
            ["time", "azimuth", "elevation", "antenna_transition"]
            + SCALAR_RAY_VARS
            + MOMENT_VARS
        )
        var_attrs = {}
        for vname in all_var_names:
            if vname in f0.variables:
                var_attrs[vname] = {
                    k: f0.variables[vname].getncattr(k)
                    for k in f0.variables[vname].ncattrs()
                    if k != "_FillValue"
                }

    # Time reference: midnight of the day
    time_ref_str = (
        f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:8]}T00:00:00Z"
    )
    time_units = f"seconds since {time_ref_str}"
    ref_dt = datetime.datetime.strptime(
        time_ref_str, "%Y-%m-%dT%H:%M:%SZ"
    ).replace(tzinfo=datetime.timezone.utc)
    ref_unix = ref_dt.timestamp()

    # ------------------------------------------------------------------
    # Read all sweeps
    # ------------------------------------------------------------------
    print("Reading sweeps...")
    all_times = []
    all_az = []
    all_el = []
    all_at = []
    all_prt = []
    all_pw = []
    all_txpow = []
    all_moments = {v: [] for v in MOMENT_VARS}

    sweep_start = []
    sweep_end = []
    n_rays_total = 0
    fixed_el_values = []

    for ppi_file in ppi_files:
        fname_short = os.path.basename(ppi_file)
        with nc.Dataset(ppi_file) as f:
            n_rays = len(f.dimensions["time"])

            # Ray Unix times → offset from day midnight
            ray_unix = _ray_unix_times(f)
            ray_times = (ray_unix - ref_unix).astype(np.float64)

            # Azimuth
            az = np.ma.getdata(f.variables["azimuth"][:]).astype(np.float32)

            # Elevation: replace sentinels from axis log
            el_raw = np.ma.getdata(f.variables["elevation"][:]).astype(
                np.float32
            )
            is_sentinel = el_raw < sentinel_threshold
            n_sent = int(is_sentinel.sum())
            if use_log:
                el_log = _interp_log_el(ray_unix, log_t, log_el)
                el = np.where(is_sentinel, el_log, el_raw)
            else:
                el = np.where(is_sentinel, 82.0, el_raw)

            median_el = float(np.median(el))
            fixed_el_values.append(round(median_el, 1))
            print(
                f"  {fname_short}: {n_rays} rays, "
                f"{n_sent} sentinel el replaced, "
                f"median el = {median_el:.2f}°"
            )

            # 1-D instrument parameter arrays
            at = np.ma.getdata(f.variables["antenna_transition"][:])
            prt = np.ma.getdata(f.variables["prt"][:]).astype(np.float32)
            pw = np.ma.getdata(f.variables["pulse_width"][:]).astype(
                np.float32
            )
            txpow = np.ma.getdata(
                f.variables["radar_measured_transmit_power_h"][:]
            ).astype(np.float32)

            # Moments — keep as masked arrays so fill is preserved
            sweep_start.append(n_rays_total)
            sweep_end.append(n_rays_total + n_rays - 1)
            n_rays_total += n_rays

            all_times.append(ray_times)
            all_az.append(az)
            all_el.append(el.astype(np.float32))
            all_at.append(at)
            all_prt.append(prt)
            all_pw.append(pw)
            all_txpow.append(txpow)

            for mv in MOMENT_VARS:
                data = f.variables[mv][:]  # masked array
                all_moments[mv].append(data)

    # Concatenate
    all_times = np.concatenate(all_times)
    all_az = np.concatenate(all_az)
    all_el = np.concatenate(all_el)
    all_at = np.concatenate(all_at)
    all_prt = np.concatenate(all_prt)
    all_pw = np.concatenate(all_pw)
    all_txpow = np.concatenate(all_txpow)
    for mv in MOMENT_VARS:
        all_moments[mv] = np.ma.concatenate(all_moments[mv], axis=0)

    n_sweeps = len(ppi_files)
    n_rays = n_rays_total
    n_range = len(range_data)

    # Fixed angle per sweep: use median elevation recovered from log
    fixed_angles = np.array(fixed_el_values, dtype=np.float32)

    # Time coverage strings
    tcs_str = (
        ref_dt + datetime.timedelta(seconds=float(all_times[0]))
    ).strftime("%Y-%m-%dT%H:%M:%SZ")
    tce_str = (
        ref_dt + datetime.timedelta(seconds=float(all_times[-1]))
    ).strftime("%Y-%m-%dT%H:%M:%SZ")

    # Output filename: use timestamp from first file, change scan type to vad
    m = re.search(r"_(\d{8}-\d{6})_", os.path.basename(ppi_files[0]))
    time_part = m.group(1) if m else datestr + "-000000"
    # Extract version suffix from filename (e.g. v1.0.1)
    m2 = re.search(r"_(v\d+\.\d+\.\d+)\.nc$", os.path.basename(ppi_files[0]))
    ver_suffix = m2.group(1) if m2 else "v1.0.1"
    out_fname = (
        f"ncas-mobile-ka-band-radar-1_cao_{time_part}_vad_l1_{ver_suffix}.nc"
    )

    if output_dir is None:
        output_dir = day_dir
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, out_fname)

    print(f"\nWriting {out_path}")
    print(f"  {n_sweeps} sweeps, {n_rays:,} rays, {n_range} range gates")

    # ------------------------------------------------------------------
    # Write output NetCDF
    # ------------------------------------------------------------------
    def _str_to_char(s, length=32):
        """Pack a string into a numpy S1 array of fixed length."""
        arr = np.zeros(length, dtype="S1")
        for i, c in enumerate(s[:length]):
            arr[i] = c.encode("ascii")
        return arr

    with nc.Dataset(out_path, "w", format="NETCDF4") as out:

        # Dimensions
        out.createDimension("time", None)  # unlimited
        out.createDimension("range", n_range)
        out.createDimension("sweep", n_sweeps)
        out.createDimension("string_length", 32)
        out.createDimension("frequency", 1)

        # Compression / chunking settings to match reference VAD files
        _zlib = dict(zlib=True, complevel=4, shuffle=True)

        # ---- time ----
        v = out.createVariable("time", "f8", ("time",),
                               chunksizes=(512,), **_zlib)
        v.units = time_units
        v.standard_name = "time"
        v.calendar = "gregorian"
        v.long_name = "time_since_time_reference"
        v.comment = " "
        v[:] = all_times

        # ---- range ----
        v = out.createVariable("range", "f4", ("range",), **_zlib)
        for k, val in range_attrs.items():
            setattr(v, k, val)
        v[:] = range_data

        # ---- azimuth ----
        v = out.createVariable("azimuth", "f4", ("time",),
                               chunksizes=(1024,), **_zlib)
        for k, val in var_attrs.get("azimuth", {}).items():
            setattr(v, k, val)
        v[:] = all_az

        # ---- elevation (corrected) ----
        v = out.createVariable("elevation", "f4", ("time",),
                               chunksizes=(1024,), **_zlib)
        for k, val in var_attrs.get("elevation", {}).items():
            setattr(v, k, val)
        v[:] = all_el

        # ---- antenna_transition ----
        v = out.createVariable("antenna_transition", "i1", ("time",),
                               chunksizes=(1024,), **_zlib)
        for k, val in var_attrs.get("antenna_transition", {}).items():
            setattr(v, k, val)
        v[:] = all_at

        # ---- moment variables ----
        for mv in MOMENT_VARS:
            v = out.createVariable(
                mv, "f4", ("time", "range"), fill_value=-9999.0,
                chunksizes=(1, n_range), **_zlib
            )
            for k, val in var_attrs.get(mv, {}).items():
                setattr(v, k, val)
            v[:] = all_moments[mv]

        # ---- sweep metadata ----
        v = out.createVariable("sweep_number", "i4", ("sweep",))
        v[:] = np.arange(n_sweeps, dtype=np.int32)

        v = out.createVariable("fixed_angle", "f4", ("sweep",))
        v[:] = fixed_angles

        v = out.createVariable("sweep_start_ray_index", "i4", ("sweep",))
        v[:] = np.array(sweep_start, dtype=np.int32)

        v = out.createVariable("sweep_end_ray_index", "i4", ("sweep",))
        v[:] = np.array(sweep_end, dtype=np.int32)

        v = out.createVariable("sweep_mode", "S1", ("sweep", "string_length"))
        v.comment = (
            'Options are: "sector", "coplane", "rhi", "vertical_pointing", '
            '"idle", "azimuth_surveillance", "elevation_surveillance", '
            '"sunscan", "pointing", "manual_ppi", "manual_rhi"'
        )
        v.units = ""
        mode_arr = np.zeros((n_sweeps, 32), dtype="S1")
        for i in range(n_sweeps):
            for j, c in enumerate("vad"):
                mode_arr[i, j] = c.encode("ascii")
        v[:] = mode_arr

        # ---- per-ray instrument parameters ----
        v = out.createVariable("prt", "f4", ("time",),
                               chunksizes=(512,), **_zlib)
        for k, val in var_attrs.get("prt", {}).items():
            setattr(v, k, val)
        v[:] = all_prt

        v = out.createVariable("pulse_width", "f4", ("time",),
                               chunksizes=(512,), **_zlib)
        for k, val in var_attrs.get("pulse_width", {}).items():
            setattr(v, k, val)
        v[:] = all_pw

        v = out.createVariable(
            "radar_measured_transmit_power_h", "f4", ("time",),
            chunksizes=(512,), **_zlib
        )
        for k, val in var_attrs.get(
            "radar_measured_transmit_power_h", {}
        ).items():
            setattr(v, k, val)
        v[:] = all_txpow

        # ---- static scalars ----
        v = out.createVariable("latitude", "f8")
        v.long_name = "Latitude"
        v.units = "degrees_north"
        v.standard_name = "Latitude"
        v[:] = latitude

        v = out.createVariable("longitude", "f8")
        v.long_name = "Longitude"
        v.units = "degrees_east"
        v.standard_name = "Longitude"
        v[:] = longitude

        v = out.createVariable("altitude", "f8")
        v.long_name = "Altitude"
        v.units = "metres"
        v.standard_name = "Altitude"
        v.positive = "up"
        v.comment = (
            "Altitude of the centre of rotation of the antenna above the "
            "geoid using the WGS84 ellipsoid and EGM2008 geoid model"
        )
        v[:] = altitude

        # ---- string variables ----
        for vname, val_str in [
            ("time_coverage_start", tcs_str),
            ("time_coverage_end", tce_str),
            ("time_reference", time_ref_str),
        ]:
            v = out.createVariable(vname, "S1", ("string_length",))
            if vname == "time_coverage_start":
                v.long_name = "UTC time of first ray in the file"
            elif vname == "time_coverage_end":
                v.long_name = "UTC time of last ray in the file"
            else:
                v.long_name = "UTC time reference"
            v.units = "unitless"
            v[:] = _str_to_char(val_str)

        # ---- volume number ----
        v = out.createVariable("volume_number", "i4")
        v.long_name = "Volume number"
        v.units = "unitless"
        v[:] = 0

        # ---- frequency ----
        v = out.createVariable("frequency", "f4", ("frequency",))
        v.standard_name = "radiation_frequency"
        v.long_name = "frequency of transmitted radiation"
        v.units = "s-1"
        v.meta_group = "instrument_parameters"
        v[:] = frequency

        # ---- azimuth_correction ----
        v = out.createVariable("azimuth_correction", "f4")
        v.long_name = "azimuth correction applied"
        v.units = "degrees"
        v.meta_group = "geometry_correction"
        v.comment = (
            "Azimuth correction applied. "
            "North angle relative to instrument home azimuth."
        )
        v[:] = az_correction

        # ---- global attributes ----
        skip_attrs = {
            "time_coverage_start",
            "time_coverage_end",
            "history",
            "last_revised_date",
            "date_created",
        }
        for k, val in global_attrs.items():
            if k not in skip_attrs:
                setattr(out, k, val)

        out.time_coverage_start = tcs_str
        out.time_coverage_end = tce_str

        now_str = datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S %Y")
        now_iso = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        user = getpass.getuser()
        host = socket.gethostname()
        history_new = (
            f"{now_str} - user:{user} machine:{host} "
            f"CCREST-M: Multi-sweep vad conversion, "
            f"revised_northangle={az_correction:.1f}deg, "
            f"{n_sweeps} sweeps"
        )
        existing_hist = global_attrs.get("history", "").strip()
        out.history = (
            f"{existing_hist}\n{history_new}" if existing_hist else history_new
        )
        out.last_revised_date = now_iso
        out.date_created = now_iso

    print(f"Done.")
    print(f"Output: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild PPI files with sentinel elevation as a multisweep VAD "
            "file, recovering elevation values from the axis log."
        )
    )
    parser.add_argument("datestr", help="Date string YYYYMMDD")
    parser.add_argument(
        "l1_path",
        nargs="?",
        default=L1_BASE,
        help=f"L1 root directory (default: {L1_BASE})",
    )
    parser.add_argument(
        "--logp-path",
        default=LOGP_BASE,
        help=f"Axis log directory (default: {LOGP_BASE})",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: same as input day directory)",
    )
    parser.add_argument(
        "--sentinel-threshold",
        type=float,
        default=SENTINEL_THRESHOLD,
        help=(
            "Elevation values below this are treated as sentinels "
            f"(default: {SENTINEL_THRESHOLD})"
        ),
    )
    args = parser.parse_args()

    rebuild_vad(
        args.datestr,
        args.l1_path,
        args.logp_path,
        output_dir=args.output_dir,
        sentinel_threshold=args.sentinel_threshold,
    )


if __name__ == "__main__":
    main()
