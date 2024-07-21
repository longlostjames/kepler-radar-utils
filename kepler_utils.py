#!/usr/bin/env python
# coding: utf-8

# ==========================================================================
# Module for processing mmclx radar files from Kepler (MIRA-35) radar
# Author: Chris Walden, UK Research & Innovation and
#                       National Centre for Atmospheric Science
# Last modified: 02-10-2023
# ==========================================================================

"""Module for processing mmclx radar data from Kepler (MIRA-35) radar."""

module_version = 0.2;

import datetime, cftime

import netCDF4 as nc4
import numpy as np

import getpass, socket
import pyart

from pyart.config import FileMetadata, get_fillvalue, get_metadata
from pyart.core.radar import Radar
from pyart.io.common import _test_arguments, make_time_unit_str

import yaml

import os, fnmatch
import re

import gzip

from io import StringIO

def read_mira35_mmclx_hsrhi(mmclxfiles, **kwargs):
    """
    Read a set of single-sweep netCDF mmclx files from MIRA-35 radar recorded as part of HSRHI scan strategy.

    Parameters
    ----------
    mmclxfiles : list(str)
        List of names of mmclx netCDF files to read data from.

    Returns
    -------
    radar : Radar
        Radar object.
    """

    mmclxfiles.sort();

    nsweep = len(mmclxfiles);

def read_mira35_mmclx_vpt_multi(mmclxfiles, **kwargs):
    """
    Read a set of netCDF mmclx files from MIRA-35 radar recorded as separate vertical pointing sweeps.

    Parameters
    ----------
    mmclxfiles : list(str)
        List of names of mmclx netCDF files to read data from.

    Returns
    -------
    radar : Radar
        Radar object.
    """

    mmclxfiles.sort();

    nsweep = len(mmclxfiles);



def read_mira35_mmclx(filename, gzip_flag=False, revised_northangle=55.7, **kwargs):
    """
    Read a netCDF mmclx file from MIRA-35 radar.

    Parameters
    ----------
    filename : str
        Name of mmclx netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.
    """
    # List 
    # time, range, fields, metadata, scan_name, latitude, longitude, altitude, altitude_agl,
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index, sweep_end_ray_index, rays_per_sweep,
    # target_scan_rate, rays_are_indexed, ray_angle_res,
    # azimuth, elevation, gate_x, gate_y, gate_z, gate_longitude, gate_latitude, projection, gate_altitude,
    # scan_rate, antenna_transition, 
    # instrument_parameters
    # radar_calibration
    # OK ngates
    # OK nrays
    # OK nsweeps
    
    # This routine only applies to fixed platforms
    # The following are not required for a fixed platform
    rotation = None;
    tilt = None;
    roll = None;
    drift = None;
    heading = None;
    pitch = None;
    georefs_applied = None;
     
    # -------------------------
    # test for non empty kwargs
    # -------------------------
    _test_arguments(kwargs)

    # --------------------------------
    # create metadata retrieval object
    # --------------------------------
    filemetadata = FileMetadata('mmclx')

    # -----------------
    # Open netCDF4 file
    # -----------------
    print(f"gzip_flag={gzip_flag}")
    if gzip_flag:
        gz = gzip.open(filename)
        ncobj = nc4.Dataset('dummy', mode='r', memory=gz.read())
    else:
        ncobj = nc4.Dataset(filename)
    
    nrays = len(ncobj.dimensions["time"]);
    ngates = len(ncobj.dimensions["range"]);
    nsweeps = 1; # We only have single sweep files 
    
    ncvars = ncobj.variables;

    # --------------------------------
    # latitude, longitude and altitude
    # --------------------------------
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    z = StringIO(ncobj.getncattr('Latitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lat','zs2'])
    if z1['zs2']==b'S' and z1['lat']>0: 
        latitude['data'] = np.array([-z1['lat']],dtype='f4');
    else:
        latitude['data'] = np.array([z1['lat']],dtype='f4');  
        
    z = StringIO(ncobj.getncattr('Longitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lon','zs2'])
    
    if z1['zs2']==b'W' and z1['lon']>0: 
        longitude['data'] = np.array([-z1['lon']],dtype='f4');
    else:
        longitude['data'] = np.array([z1['lon']],dtype='f4');  
    
    z = StringIO(ncobj.getncattr('Altitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['alt','zs2'])
    altitude['data'] = np.array([z1['alt']],dtype='f4')

    # Original mmclx metadata
    # -----------------------
    # convention: "CF-1.0"
    # location: "Chilbolton"
    # Altitude:
    # Latitude:
    # Longitude:
    # system: "B3XC"
    # title: "MIRA Cloud Radar Data"
    # institution: "National Centre for Atmospheric Science (NCAS)"
    # source: "220520_101223.pds"
    # reference: "Ka Band Cloud Radar MIRA 35, METEK GmbH www.metek.de"

    # Most important information from the header structure
    # ----------------------------------------------------
    # nfft: "number_of_fft_points"
    # NyquistVelocity: 
    # nave: "number_of_spectral_averages"
    # ovl: "overlapping_ffts_flag"
    # zrg: "number_of_range_gates"
    # rg0: "number_of_lowest_range_gate"
    # drg: "range_resolution"
    # lambda: "wavelength"
    # range: "range_from_antenna_to_centre_of_each_range_gate"

    # Most important information from the SRVI struture given at each dwell time
    # ---------------------------------------------------------------------------
    # time: "seconds since 1970-01-01T00:00:00Z"
    # # microsec: 
    # tpow: "average_transmit_power"
    # npw1: "noise_power_co-channel"
    # npw2: "noise_power_cross-channel"
    # cpw1: "snr_of_calibration_signal_co-channel"
    # cpw2: "snr_of_calibration_signal_cross-channel"
    # grst: "general_radar_state"
    # azi: "azimuth"
    # elv: "elevation"
    # aziv: "azimuth_angle_velocity"
    # northangle: "north_angle"
    # elvv: "elevation_angle_velocity"
    # LO_Frequency: "transmit_frequency"
    # DetuneFine: "detuning_of_local_scillator_detected_from_Txn_footprint"

    metadata_keymap = {
        "location": "platform",
        "Longitude": "longitude",
        "Latitude": "latitude",
        "Altitude": "altitude",
        "system": "instrument_serial_nnmber",
        "title": "title",
        "institution": "institution",
        "reference": "reference",
    }

 # time, range, fields, metadata, scan_name, latitude, longitude, altitude, altitude_agl,
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index, sweep_end_ray_index, rays_per_sweep,
    # target_scan_rate, rays_are_indexed, ray_angle_res,
    # azimuth, elevation, gate_x, gate_y, gate_z, gate_longitude, gate_latitude, projection, gate_altitude,
    # scan_rate, antenna_transition, 
    # None rotation, tilt, roll, drift, heading, pitch
    # ?? georefs_applied
    # instrument_parameters
    # radar_calibration



    variables_keymap = {
        "elv": "elevation",
        "elvv": "elevation_angle_velocity",
        "azi": "azimuth_angle",
        "aziv": "azimuth_angle_velocity",
        "nfft": "nfft",
        "nave": "nave",
        "prf": "prf",
        "rg0": "rg0",
    }

    fields_keymap = {
        "Zg": "DBZ",
        "VELg": "VEL",
        "RMSg": "WIDTH",
        "LDRg": "LDR",
        "SNRg": "SNR",
        "RHO": "RHOHX",
        "DPS": "PHIHX",
    }

    variables = list(variables_keymap.keys())

    # SNRg, VELg, RMSg, LDRg, NPKg, SNRg
    #RHO DPS PHOwav LDRnormal
    # HSDco HSDcx Zg ISDRco ISDRcx 
    # MRMco MRMcx RadarConst SNRCorFaCo SNRCorFaCx SKWg

    # Variables
    # nfft
    # prf
    # NyquistVelocity
    # nave
    
    # metadata
    #metadata = filemetadata("metadata")
    #metadata_mapping = {
    #    "vcp-value": "vcp",
    #    "radarName-value": "radar_name",
    #    "ConversionPlugin": "conversion_software",
    #}
    #for netcdf_attr, metadata_key in metadata_mapping.items():
    #    if netcdf_attr in dset.ncattrs():
    #        metadata[metadata_key] = dset.getncattr(netcdf_attr)

    # --------
    # metadata
    # --------
    metadata = filemetadata('metadata')
    for k in ['institution', 'title', 'used_algorithms']:
        if k in ncobj.ncattrs(): 
            metadata[k] = ncobj.getncattr(k)

    metadata['instrument_name']='ncas-radar-mobile-ka-band-1'

    # ------------------------------------------
    # sweep_start_ray_index, sweep_end_ray_index
    # ------------------------------------------
    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_end_ray_index = filemetadata("sweep_end_ray_index")
    sweep_start_ray_index["data"] = np.array([0], dtype="int32")
    sweep_end_ray_index["data"] = np.array([nrays - 1], dtype="int32")

    # ------------
    # sweep number
    # ------------
    sweep_number = filemetadata("sweep_number")
    sweep_number["data"] = np.array([0], dtype="int32")

    # -----------------------
    # sweep_mode, fixed_angle
    # -----------------------
    sweep_modes = {'ppi' : 'ppi', 'rhi' : 'rhi', 'vert' : 'vertical_pointing','man' : 'manual_rhi'}

    sweep_mode = filemetadata("sweep_mode")

    print(filename.lower());

    scan_name = None;
    sweep_mode["data"] = np.array(1 * [None]);

    for key, value in sweep_modes.items():
        if key in filename.lower(): 
            scan_name = value;
            sweep_mode["data"] = np.array(1 * [value]);
            sweep_mode["data"][0] = value;
            break;
    
    print(scan_name);

    #    #fixed_angles = {'ppi' : ncvars['elv'][0], 'rhi' : ncvars['azi'][0]+ncvars['northangle'][0], 'vertical_pointing' : ncvars['elv'][0], "manual_rhi" : ncvars['azi'][0]}
    #fixed_angles = {'ppi' : ncvars['elv'][10], 'rhi' : ncvars['azi'][10]+revised_northangle, 'vertical_pointing' : ncvars['elv'][10], "manual_rhi" : ncvars['azi'][10]}

    fixed_angle = filemetadata("fixed_angle")

    if scan_name in ['rhi','manual_rhi']:
        fixed_angle_value = np.round((ncvars['azi'][5]+revised_northangle)%360,2);
        fixed_angle["data"] = np.array(1 * [fixed_angle_value], dtype='f'); 
    elif scan_name in ['ppi','vertical_pointing']:
        fixed_angle_value = np.round(ncvars['elv'][5],2);
        fixed_angle["data"] = np.array(1 * [fixed_angle_value], dtype='f');
    else:
        fixed_angle["data"] = np.array(1 * [None]);
    
    print(fixed_angle["data"])


    # time
    # interpolate between the first and last timestamps in the Time variable
    time = filemetadata('time')
    
    dtime = cftime.num2pydate(ncvars['time'][:],'seconds since 1970-01-01 00:00:00')

    for idx, x in np.ndenumerate(dtime):
        dtime[idx]=x.replace(microsecond=ncvars['microsec'][idx])
            
    base_time = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    
    time['units'] = make_time_unit_str(base_time);  
    time['data']  = cftime.date2num(dtime,time['units']);

    # range
    # -----
    _range = filemetadata('range');
    _range['data'] = ncvars['range'][:];
    _range['units'] = 'metres';
    #_range['metres_to_centre_of_first_gate'] = _range['data'][0];
    _range['proposed_standard_name'] = "projection_range_coordinate";
    _range['long_name'] = "distance to centre of each range gate";
    # assuming the distance between all gates is constant, may not
    # always be true.
    #_range['metres_between_gates'] = (_range['data'][1] - _range['data'][0])

    # azimuth, elevation
    # ------------------
    #ray_duration0 = time['data'][1]-time['data'][0]; 
    #ray_duration = np.insert(time['data'][1:]-time['data'][:-1],0,ray_duration0);
    ray_duration = time['data'][1:]-time['data'][:-1];

    # scan rate
    # ---------
    scan_rates = {'ppi' : ncvars['aziv'][:], 'rhi' : ncvars['elvv'][:]}

    scan_rate = filemetadata("scan_rate")
    antenna_transition = filemetadata("antenna_transition")
    target_scan_rate = filemetadata("target_scan_rate")

    if scan_name in  ['ppi','rhi']:
        scan_rate['data'] = scan_rates[scan_name];
        scan_rate['units']='degrees_per_second';
        scan_rate['_FillValue'] = -9999.
        scan_rate['long_name'] = 'antenna angle scan rate';
        antenna_transition['data'] = np.where(abs(scan_rate['data'])<0.01,1,0).astype('int8');
        antenna_transition['long_name'] = "antenna is in transition between sweeps" ;
        antenna_transition['units'] = "" ;
        antenna_transition['_FillValue'] = -128 ;
        antenna_transition['comment'] = "1 if antenna is in transition, 0 otherwise" ;

        target_scan_rate = filemetadata("target_scan_rate")
        target_scan_rate["data"] = np.array([4.0], dtype="f4")
        scanning_indices = np.where(antenna_transition['data']==0)[0];
        target_scan_rate['data'] = np.mean(scan_rate['data'][scanning_indices]);
        at_speed = np.where(abs(scan_rate['data']-target_scan_rate['data'])<0.06)[0];
        target_scan_rate['data'] = np.round(np.mean(scan_rate['data'][at_speed]),2);
        target_scan_rate['long_name'] = 'target scan rate for sweep';
        target_ray_duration = np.round(np.mean(ray_duration),3);
        ok_duration = np.where(ray_duration/target_ray_duration<1.5);
        target_ray_duration = np.round(np.mean(ray_duration[ok_duration]),3);

        print(target_ray_duration);

        ray_duration = np.insert(ray_duration,0,target_ray_duration);

        long_duration  = np.where(ray_duration/target_ray_duration>1.5)[0];
    
        print(f'long-duration indices:{long_duration}')

    else:
        scan_rate  = None;
        antenna_transition = None;
        target_scan_rate = None;
    
    
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = (ncvars['azi'][:]+revised_northangle) % 360;
    elevation['data'] = ncvars['elv'][:];

    print("i am here")

    if scan_name in  ['ppi','rhi']:
        #azimuth['data'] = (ncvars['azi'][:]+ncvars['northangle'][:]) % 360;
        azimuth['data'] -= 0.5*ray_duration * ncvars['aziv'][:] ;
        azimuth['units'] = "degrees";
        azimuth['proposed_standard_name'] = "sensor_to_target_azimuth_angle";
        azimuth['long_name'] = "sensor to target azimuth angle";

        # Special case for long-duration glitches
        azimuth['data'][long_duration] = (ncvars['azi'][long_duration]+revised_northangle) % 360;
        azimuth['data'][long_duration] -= 0.5*target_ray_duration * ncvars['aziv'][long_duration];
    

        elevation['data'] -= 0.5*ray_duration * ncvars['elvv'][:];
        elevation['units'] = "degrees";
        elevation['proposed_standard_name'] = "sensor_to_target_elevation_angle";
        elevation['long_name'] = "sensor to target elevation angle";

        # Special case for long-duration glitches
        elevation['data'][long_duration] = ncvars['elv'][long_duration];
        elevation['data'][long_duration] -= 0.5*target_ray_duration * ncvars['elvv'][long_duration];


    metadata['time_coverage_start'] = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
    metadata['time_coverage_end'] = datetime.datetime.strftime(dtime[-1],'%Y-%m-%dT%H:%M:%SZ');

    print("now i am here")

    # ------
    # fields
    # ------
    # mmclx files contain the following: 
    # where x="g" (global),"" (hydrometeors),"plank" (plankton),"rain","cl" 
    # SNRx - USE THIS
    # VELx - USE THIS
    # RMSx - USE THIS
    # LDRx - USE THIS

    # LWC (liquid water content of peaks classified as rain)
    # RR (rain rate of peaks classified as rain)
    # Ze (equivalent radar reflectivity of all hydrometeors)

    # Zg (equivalent radar relectivity of all targets (global)) - USE THIS
    
    # Z (similar to Ze but with Mie correction applied, and a pressure and temperature dependent value of |Kw|^2
    # MeltHeiDet (melting layer height detected from LDR if detected, else -1)
    # MeltHeiDb (melting layer height as deduced from external sources)
    # MeltHei (melting layer height from both sources)
    # TEMP (temperature profile faked form external sources)

    # ISDRco (ratio between power integrated over the largest peak and the noise level - to indicate saturation) - MAYBE USE
    # ISDRcx (cross-polar channel version of ISDRco) - MAYBE USE
    # RadarConst (radar constant related to 5km range, Zx = RadarConst*SNRx*(range/5 km)^2), with range at the centre of each range gate) - MAYBE USE
    # SNRCorFaCo () - MAYBE USE
    # SNRCorFaCx () - MAYBE USE
    

    fields = {}


    #linear_to_log_fields = {"Zg"}
    #for key, value in fields_keymap.items():
    #    print(key)
    #    if key in ncvars:
    #        field_name = value;
    #        field_dic = filemetadata(field_name);
    #        field_dic['_FillValue'] = get_fillvalue();
    #        field_dic['units'] = ncvars[key].units;
    #        if key in linear_to_log_fields:
    #            field_dic['data'] = 10.0*np.log10(ncvars[key][:]);
    #        else:
    #            field_dic['data'] = ncvars[key][:];
    #        fields[field_name] = field_dic;

    if "Zg" in ncvars:
        field_name = fields_keymap['Zg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dBZ'
        field_dic['data'] = 10.0*np.log10(ncvars['Zg'][:]);
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "radar equivalent reflectivity factor";
        field_dic['standard_name'] = "equivalent_reflectivity_factor";
        field_dic['proposed_standard_name'] =  "radar_equivalent_reflectivity_factor";   
        fields[field_name] = field_dic
    else:
        print("Zg does not exist")

    print("Done Zg")

    if "VELg" in ncvars:
        field_name = fields_keymap['VELg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['VELg'][:];
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "radial velocity of scatterers away from instrument";
        field_dic['standard_name'] = "radial_velocity_of_scatterers_away_from_instrument";
        field_dic['field_folds'] = 'true'
        vfold = ncvars['prf'][:] * ncvars['lambda'][:] / 4.0;
        field_dic['field_limit_lower'] = -vfold;
        field_dic['field_limit_upper'] = vfold;
        fields[field_name] = field_dic
    else:
        print("VELg does not exist")

    print("Done VELg")


    if "RMSg" in ncvars:
        field_name = fields_keymap['RMSg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['RMSg'][:];
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "radar doppler spectrum width";
        field_dic['proposed_standard_name'] = "radar_doppler_spectrum_width";
        fields[field_name] = field_dic
    else:
        print("RMSg does not exist")

    print("Done RMSg");

    if "LDRg" in ncvars:
        field_name = fields_keymap['LDRg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = 10.0*np.log10(ncvars['LDRg'][:]);
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "radar linear depolarization ratio";
        field_dic['proposed_standard_name'] = "radar_linear_depolarization_ratio";
        fields[field_name] = field_dic
    else:
        print("LDRg does not exist")

    print("Done LDRg");

    if "SNRg" in ncvars:
        field_name = fields_keymap['SNRg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = 10.0*np.log10(ncvars['SNRg'][:]);
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "radar signal to noise ratio";
        field_dic['proposed_standard_name'] = "radar_signal_to_noise_ratio";
        fields[field_name] = field_dic
    else:
        print("SNRg does not exist")
    
    print("Done SNRg");


    if "RHO" in ncvars:
        field_name = fields_keymap['RHO']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = ''
        field_dic['data'] = ncvars['RHO'][:];
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "co- to cross-polar correlation ratio for horizontal transmitted polarization";
        field_dic['proposed_standard_name'] = "co_to_cross_polar_correlation_ratio_h"
        fields[field_name] = field_dic
    else:
        print("RHO does not exist")

    print("Done RHO")

    if "DPS" in ncvars:
        field_name = fields_keymap['DPS']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degrees'
        field_dic['data'] = ncvars['DPS'][:];
        isnan = np.isnan(field_dic['data'][:]);
        field_dic['data'][isnan] = field_dic['_FillValue'];
        if scan_name in  ['ppi','rhi']:
            field_dic['data'][long_duration,:] = field_dic['_FillValue'];
            field_dic['data'][long_duration-1,:] = field_dic['_FillValue'];
        field_dic['long_name'] =  "cross-polar differential phase";
        field_dic['proposed_standard_name'] = "cross_polar_differential_phase"
        fields[field_name] = field_dic
    else:
        print("DPS does not exist")

    print("Done DPS")

    # instrument_parameters
    instrument_parameters = {}
    radar_parameters = {}

    radar_calibration = {}

    hrd_variables = {};

    if 'hrd' in ncobj.ncattrs():
        hrd_attribute_value = ncobj.getncattr('hrd')

        print(hrd_attribute_value);

        if type(hrd_attribute_value) == str:
            # Parse the attribute value while ignoring lines starting with "DESCR"
            for line in hrd_attribute_value.split('\n'):
                if not line.startswith("DESCR") and ':' in line:
                    key, value = line.split(':', 1)
                    print(key, value)
                    hrd_variables[key.strip()] = convert_to_numerical(value.strip())
        else:
            print("Value of 'hrd' attribute is not a string.")
    else:
        print("'hrd' attribute not found in the NetCDF file.")


    print(hrd_variables)


    if "prf" in ncvars:
        dic = filemetadata("prt")
        prt = 1.0/ncvars["prf"][:]
        dic["data"] = np.ones((nrays,), dtype="float32") * prt
        instrument_parameters["prt"] = dic

    if "PULSE_WIDTH" in hrd_variables:
        dic = filemetadata("pulse_width")
        pulse_width = hrd_variables["PULSE_WIDTH"]
        dic["data"] = np.ones((nrays,), dtype="float32") * pulse_width
        instrument_parameters["pulse_width"] = dic

    if "tpow" in ncvars:
        dic = filemetadata("radar_measured_transmit_power_h")
        txpower = 10.0*np.log10(ncvars["tpow"][:])+30;
        dic["data"] = txpower
        instrument_parameters["radar_measured_transmit_power_h"] = dic
        
    #if "NyquistVelocity-value" in dset.ncattrs():
    #    dic = filemetadata("nyquist_velocity")
    #    nyquist_velocity = float(dset.getncattr("NyquistVelocity-value"))
    #    dic["data"] = np.ones((nrays,), dtype="float32") * nyquist_velocity
    #    instrument_parameters["nyquist_velocity"] = dic

    #if "Beamwidth" in dset.variables:
    #    dic = filemetadata("radar_beam_width_h")
    #    dic["data"] = dset.variables["Beamwidth"][:]
    #    instrument_parameters["radar_beam_width_h"] = dic

    
    ncobj.close()
    if gzip_flag:
        gz.close();

    radar = Radar(
        time,
        _range,
        fields,
        metadata,
        scan_name,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        target_scan_rate=target_scan_rate,
        scan_rate=scan_rate,
        antenna_transition=antenna_transition,
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration
    )

    print("last line")
    return radar

def convert_to_numerical(value):
    try:
        # Try converting the value to float
        return float(value)
    except ValueError:
        try:
            # If conversion to float fails, try converting to int
            return int(value)
        except ValueError:
            return value


# ===================
# CONVERSION ROUTINES
# ===================

def convert_kepler_mmclx2l1(infile,outpath,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):

    """This routine converts mmclx data from the NCAS Mobile Ka-band Radar (Kepler) to Level 1 (cfradial) data, compliant with the 
    NCAS Radar Data Standard v1.0.0.

    Metadata are added using information in two YAML files the yaml_project_file, and yaml_instrument_file.

    :param infile: Full path of NetCDF Level 0b mmclx data file, e.g. `<path-to-file>/20220907_071502.ppi.mmclx`
    :type infile: str

    :param outpath: Path where NetCDF Level 1 output file will be written
    :type outfile: str

    :param yaml_project_file: Full path of YAML file containing project-specific metadata
    :type yaml_project_file: str

    :param yaml_instrument_file: Full path of YAML file containing instrument-specific metadata
    :type yaml_instrument_file: str

    :param tracking_tag: AMOF tracking tag for the project
    "type tracking_tag: str
    
    :param data_version: Version of data product in the format `n.m.p`, where n (major version), m (minor revision) and p (patch) are integers.
    :type data_version: str
    """

    instrument_tagname = "ncas-radar-mobile-ka-band-1"

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

    radar_name = instrument["instrument_name"].lower();

    print(radar_name);

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name];

    print(project_instrument);

    location = project_instrument['platform']['location'].lower();
    
    RadarDataset = read_mira35_mmclx(infile);

    scan_name = RadarDataset.scan_name

    file_timestamp = datetime.datetime.strptime(RadarDataset["time_coverage_start"][:],'%Y-%m-%dT%H:%M:%SZ');

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    outfile = os.path.join(outpath,'{}_{}_{}_{}_l1_v{}.nc'.format(radar_name,location,dtstr,scan_name.replace('_','-',1),data_version));

    # ---------------------------------
    # Use PyART to create CfRadial file
    # ---------------------------------
    pyart.io.write_cfradial(outfile, RadarDataset, format='NETCDF4', time_reference=True)

    cfradial_add_ncas_metadata(outfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version)


    # -----------------------
    # Update history metadata
    # -----------------------
    DS = nc4.Dataset(outfile,'r+');

    user = getpass.getuser()

    updttime = datetime.datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: kepler_utils.convert_kepler_mmclx2l1"
    + " version:" + str(module_version));

    DS.history = history + "\n" + DS.history;

    DS.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DS.close();

    return

def cfradial_get_bbox(cfradfile):
    print(cfradfile)
    Radar = pyart.io.read_cfradial(cfradfile);
    latmin = np.min(Radar.gate_latitude['data']);
    lonmin = np.min(Radar.gate_longitude['data']);
    latmax = np.max(Radar.gate_latitude['data']);
    lonmax = np.max(Radar.gate_longitude['data']); 
    print(latmin,latmax,lonmin,lonmax)
    boundingbox = f"Bounding box: {latmin:.2f}N {lonmin:.2f}E, {latmax:.2f}N {lonmax:.2f}E"
    return boundingbox


def cfradial_add_ncas_metadata_v0(cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    # -------------------------------------------------------
    # Read cfradial file to add NCAS metadata
    # -------------------------------------------------------
    print("in ncas metadata")
    DS = nc4.Dataset(cfradfile,'r+');

    DS.product_version = f"v{data_version}";
    DS.processing_level = "1" ;

    DS.licence = project_instrument["data_licence"];
    DS.acknowledgement = project_instrument["acknowledgement"];

    DS.platform = project_instrument["platform"]["location"];
    DS.platform_type = project_instrument["platform"]["type"];
    DS.location_keywords = project_instrument["platform"]["location_keywords"];

    DS.deployment_mode = project_instrument["platform"]["deployment_mode"];

    DS.title = project_instrument["title"];

    DS.creator_name = project_instrument["data_creator"]["name"];
    DS.creator_email = project_instrument["data_creator"]["email"];
    DS.creator_url = project_instrument["data_creator"]["pid"];
    DS.institution = project_instrument["data_creator"]["institution"];
    DS.instrument_name = instrument["instrument_name"];
    DS.instrument_software = project_instrument["instrument_software"]["name"];
    DS.instrument_software_version = project_instrument["instrument_software"]["version"];
    DS.instrument_manufacturer = instrument['instrument_manufacturer'];
    DS.instrument_model = instrument['instrument_model'];
    DS.instrument_serial_number = instrument['instrument_serial_number'];
    DS.instrument_pid = instrument['instrument_pid']

    DS.references = instrument['references'];
    #DS.source = "NCAS Mobile Ka-band Radar (Kepler)";
    #DS.comment = "";
    DS.project = project["project_name"];
    DS.project_principal_investigator = project["principal_investigator"]["name"];
    DS.project_principal_investigator_email = project["principal_investigator"]["email"];
    DS.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DS.processing_software_url = "";
    DS.processing_software_version = "";

    #DS.time_coverage_start = datetime.datetime.strftime(dt_start,'%Y-%m-%dT%H:%M:%SZ');
    #DS.time_coverage_end = datetime.datetime.strftime(dt_end,'%Y-%m-%dT%H:%M:%SZ');

    print('getting bbox')
    print(cfradfile)
    str = cfradial_get_bbox(cfradfile)
    print(str)
    DS.geospatial_bounds = str;


    #DS.geospatial_bounds = "51.1450N -1.4384E";

    # -------------------------------------------------------
    # Now clean up some variable attributes
    # -------------------------------------------------------


    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['latitude'];
    #varout = DSout.createVariable('latitude',varin.datatype);
    #varout.standard_name = 'latitude';
    #varout.long_name = 'latitude of the antenna';
    #varout.units = 'degree_north';
    #varout[:]=51.1450;

    #varin = DSin['longitude'];
    #varout = DSout.createVariable('longitude',varin.datatype);
    #varout.standard_name = 'longitude';
    #varout.long_name = 'longitude of the antenna';
    #varout.units = 'degree_east';
    #varout[:]=-1.4384;

    #varin = DSin['height'];
    #varout = DSout.createVariable('altitude',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    #varout.units = 'm';
    #varout[:]=146.7;

    #varout = DSout.createVariable('altitude_agl',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above ground';
    #varout.units = 'm';
    #varout[:]=16.0;

    #varin = DSin['frequency'];
    #varout = DSout.createVariable('frequency',varin.datatype);
    #varout.standard_name = 'radiation_frequency';
    #varout.long_name = 'frequency of transmitted radiation';
    #varout.units = 'GHz';
    #varout[:]=varin[:];

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];


    # ---------------
    # Field variables
    # ---------------
    # These are:
    # SNR, VEL, RMS, LDR, NPK, SNRg, VELg, NPKg, RHO, DPS, RHOwav, DPSwav, 
    # HSDco, HSDcx, Ze, Zg, ISDRco, ISDRcx
    # The ones to use are:
    # SNR, VEL, RMS, LDR, NPK, RHO, DPS, HSDco, HSDcx, Ze, ISDRco, ISDRcx

    DS.close();

    return

 


def cfradial_add_instrument_parameters(mmclxfile,cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    # -------------------------------------------------------
    # Read cfradial file to add instrument parameters
    # -------------------------------------------------------
    DS = nc4.Dataset(cfradfile,'r+');

    print(cfradfile)

    if mmclxfile.endswith('.mmclx.gz'):
        gz = gzip.open(mmclxfile)
        DSin = nc4.Dataset('dummy', mode='r', memory=gz.read())
    else:
        DSin = nc4.Dataset(mmclxfile,'r');

    print('Creating frequency dimension')
    frequency = DS.createDimension("frequency", 1);

    varin = DSin['lambda'];
    lightspeed = 299792458
    tx_freq = lightspeed/varin[:];
    print(tx_freq)
    varout = DS.createVariable('frequency',varin.datatype,("frequency"));
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 's-1';
    varout[:]=tx_freq;
    print('Creating meta_group attribute')
    varout.meta_group = "instrument_parameters";

    varout = DS['radar_measured_transmit_power_h'];
    varout.long_name = "radar_measured_transmit_power_h";
    varout.units = 'dBm'
    varout.meta_group = "radar_parameters";


    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];


    DSin.close();
    DS.close();

    return


def cfradial_add_geometry_correction(cfradfile,revised_northangle):
    # -------------------------------------------------------
    # Read cfradial file to add instrument parameters
    # -------------------------------------------------------
    DS = nc4.Dataset(cfradfile,'r+');

    print(cfradfile)
    print(f'Revised northangle = {revised_northangle}')

    print('creating azimuth_correction')
    varout = DS.createVariable('azimuth_correction','float32');
    print('created')
    varout.long_name = "azimuth correction applied";
    varout.units = 'degrees'
    varout.meta_group = "geometry_correction";
    varout.comment = "Azimuth correction applied.  North angle relative to instrument home azimuth."
    varout[:] = revised_northangle;


    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];

    DS.close();

    return

def convert_kepler_cfradial2l1(infile,outpath,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):

    """This routine converts multi-sweep cfradial data from the NCAS Mobile Ka-band Radar (Kepler) to 
    Level 1 (cfradial) data, compliant with the NCAS Radar Data Standard v1.0.0.

    Metadata are added using information in two YAML files the yaml_project_file, and yaml_instrument_file.

    :param infile: Full path of NetCDF Level 1a cfradial data file, e.g. `<path-to-file>/20220907_071502_hsrhi.nc`
    :type infile: str

    :param outpath: Path where NetCDF Level 1 output file will be written
    :type outfile: str

    :param yaml_project_file: Full path of YAML file containing project-specific metadata
    :type yaml_project_file: str

    :param yaml_instrument_file: Full path of YAML file containing instrument-specific metadata
    :type yaml_instrument_file: str

    :param tracking_tag: AMOF tracking tag for the project
    "type tracking_tag: str
    
    :param data_version: Version of data product in the format `n.m.p`, where n (major version), m (minor revision) and p (patch) are integers.
    :type data_version: str
    """

    instrument_tagname = "ncas-radar-mobile-ka-band-1"

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

    radar_name = instrument["instrument_name"].lower();

    print(radar_name);

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name];

    print(project_instrument);

    location = project_instrument['platform']['location'].lower();
    
    RadarDataset = nc4.Dataset(infile);

    scan_name = RadarDataset.scan_name.lower();

    time_coverage_start = nc4.chartostring(RadarDataset['time_coverage_start']['data'][0])[()];

    file_timestamp = datetime.datetime.strptime(time_coverage_start,'%Y-%m-%dT%H:%M:%SZ');

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    outfile = os.path.join(outpath,'{}_{}_{}_{}_l1_v{}.nc'.format(radar_name,location,dtstr,scan_name.replace('_','-',1),data_version));

    if os.path.isfile(outfile):
        print("The file already exists")
    else:
        # Rename the file
        os.rename(infile,outfile);
    
    cfradial_add_ncas_metadata(outfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version)

    # -----------------------
    # Update history metadata
    # -----------------------
    updatestr = "Add NCAS metadata"
    update_history_attribute(outfile,updatestr)


    return

def cfradial_add_ncas_metadata(cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    
    instrument_tagname = "ncas-radar-mobile-ka-band-1"

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

    radar_name = instrument["instrument_name"].lower();

    print(radar_name);

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name];

    print(project_instrument);

    location = project_instrument['platform']['location'].lower();
    
    RadarDataset = nc4.Dataset(cfradfile);

    scan_name = RadarDataset.scan_name;

    print(scan_name);
    print('HERE')

    str_start = RadarDataset.variables['time_coverage_start'][:].tobytes().decode('utf-8')

    time_coverage_start = str_start[0:20];

    print('and HERE')

    str_end = RadarDataset.variables['time_coverage_end'][:].tobytes().decode('utf-8').strip()
    time_coverage_end = str_end[0:20];

    print('and now HERE')
    print(time_coverage_start)
    print(type(time_coverage_start))   
    print(len(time_coverage_start))   
    file_timestamp = datetime.datetime.strptime(time_coverage_start,'%Y-%m-%dT%H:%M:%SZ');

    print('and then HERE')

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    
    import pathlib
    outpath = pathlib.Path(cfradfile).parent.resolve();

    outfile = os.path.join(outpath,'{}_{}_{}_{}_l1_v{}.nc'.format(radar_name,location,dtstr,scan_name.replace('_','-',1),data_version));

    if os.path.isfile(outfile):
        print("The file already exists")
    else:
        # Rename the file
        os.rename(cfradfile,outfile);
    
    RadarDataset.close();
    
    # -------------------------------------------------------
    # Read cfradial file to add NCAS metadata
    # -------------------------------------------------------
    DS = nc4.Dataset(outfile,'r+');

    DS.Conventions = "NCAS-Radar-1.0 CfRadial-1.4 instrument_parameters radar_parameters geometry_correction"

    if 'version' in DS.ncattrs():
        DS.delncattr('version');

    DS.product_version = f"v{data_version}";
    DS.processing_level = "1" ;

    DS.licence = project_instrument["data_licence"];
    DS.acknowledgement = project_instrument["acknowledgement"];

    DS.platform = project_instrument["platform"]["location"];
    DS.platform_type = project_instrument["platform"]["type"];
    DS.location_keywords = project_instrument["platform"]["location_keywords"];
    if (project_instrument["platform"]["type"]=="stationary_platform"):
        DS.platform_is_mobile = "false";
    else:
        DS.platform_is_mobile = "true";


    DS.deployment_mode = project_instrument["platform"]["deployment_mode"];

    DS.title = project_instrument["title"];
    DS.source = project_instrument["source"];

    DS.creator_name = project_instrument["data_creator"]["name"];
    DS.creator_email = project_instrument["data_creator"]["email"];
    DS.creator_url = project_instrument["data_creator"]["pid"];
    DS.institution = project_instrument["data_creator"]["institution"];
    DS.instrument_name = instrument["instrument_name"];
    DS.instrument_software = project_instrument["instrument_software"]["name"];
    DS.instrument_software_version = project_instrument["instrument_software"]["version"];
    DS.instrument_manufacturer = instrument['instrument_manufacturer'];
    DS.instrument_model = instrument['instrument_model'];
    DS.instrument_serial_number = instrument['instrument_serial_number'];
    DS.instrument_pid = instrument['instrument_pid']

    DS.references = instrument['references'];
    #DS.source = "NCAS Mobile Ka-band Radar (Kepler)";
    #DS.comment = "";
    DS.project = project["project_name"];
    DS.project_principal_investigator = project["principal_investigator"]["name"];
    DS.project_principal_investigator_email = project["principal_investigator"]["email"];
    DS.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DS.processing_software_url = "";
    DS.processing_software_version = "";

    DS.time_coverage_start = time_coverage_start;
    DS.time_coverage_end = time_coverage_end;

    print('getting bbox')
    print(cfradfile)
    str = cfradial_get_bbox(outfile)
    print(str)
    DS.geospatial_bounds = str;


    if "vpt" in DS.scan_name or "VPT" in DS.scan_name or "vertical_pointing" in DS.scan_name:
        DS.featureType = 'timeSeriesProfile';

    # -------------------------------------------------------
    # Now clean up some variable attributes
    # -------------------------------------------------------
    DS['time'].comment = "";
    try:
        DS['range'].delncattr('standard_name');
    except: 
        pass
    
    DS['range'].comment = 'Range to centre of each bin';
    DS['range'].meters_to_center_of_first_gate = DS['range'][0];
    try:
        DS['azimuth'].delncattr('standard_name');
    except:
        pass
    try:
        DS['elevation'].delncattr('standard_name');
    except:
        pass
    DS['DBZ'].standard_name = 'equivalent_reflectivity_factor';
    try:
        DS['sweep_number'].delncattr('standard_name');
    except:
        pass
    try:
        DS['sweep_mode'].delncattr('standard_name');
    except:
        pass
    try:
        DS['fixed_angle'].delncattr('standard_name');
    except:
        pass
    DS['latitude'].long_name = 'latitude';
    DS['latitude'].standard_name = 'latitude';
    DS['longitude'].long_name = 'longitude';
    DS['longitude'].standard_name = 'longitude'; 
    DS['altitude'].standard_name = 'altitude';
    DS['altitude'].comment = 'Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model' 
    DS['altitude'].long_name = 'altitude';
    DS['altitude'].units = 'metres';
    try:
        DS['altitude'].delncattr('positive'); 
    except:
        pass
    DS['volume_number'].long_name = 'data volume index number';
    DS['volume_number'].units = "" ;
    #DS['volume_number']._FillValue = -9999 ;


    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['latitude'];
    #varout = DSout.createVariable('latitude',varin.datatype);
    #varout.standard_name = 'latitude';
    #varout.long_name = 'latitude of the antenna';
    #varout.units = 'degree_north';
    #varout[:]=51.1450;

    #varin = DSin['longitude'];
    #varout = DSout.createVariable('longitude',varin.datatype);
    #varout.standard_name = 'longitude';
    #varout.long_name = 'longitude of the antenna';
    #varout.units = 'degree_east';
    #varout[:]=-1.4384;

    #varin = DSin['height'];
    #varout = DSout.createVariable('altitude',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    #varout.units = 'm';
    #varout[:]=146.7;

    #varout = DSout.createVariable('altitude_agl',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above ground';
    #varout.units = 'm';
    #varout[:]=16.0;

    #varin = DSin['frequency'];
    #varout = DSout.createVariable('frequency',varin.datatype);
    #varout.standard_name = 'radiation_frequency';
    #varout.long_name = 'frequency of transmitted radiation';
    #varout.units = 'GHz';
    #varout[:]=varin[:];

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];


    # ---------------
    # Field variables
    # ---------------
    # These are:
    # SNR, VEL, RMS, LDR, NPK, SNRg, VELg, NPKg, RHO, DPS, RHOwav, DPSwav, 
    # HSDco, HSDcx, Ze, Zg, ISDRco, ISDRcx
    # The ones to use are:
    # SNR, VEL, RMS, LDR, NPK, RHO, DPS, HSDco, HSDcx, Ze, ISDRco, ISDRcx

    DS.close();

    # -----------------------
    # Update history metadata
    # -----------------------
    updatestr = "Add NCAS metadata"
    update_history_attribute(outfile,updatestr)

    return outfile

def anglicise_cfradial(cfradfile):

    # -------------------------------------------------------
    # Read cfradial file and change language localisation
    # -------------------------------------------------------

    DS = nc4.Dataset(cfradfile,'r+');

    variable = nc_file.variables['my_variable']
    variable.setncattr('attribute_name', 'new_attribute_value')

    DS.close();

    return


def amend_unitless(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'units') and var.units == 'unitless':
            var.units = ""
        if hasattr(var, 'units') and var.units == 'count':
            var.units = ""

    DS.close()

def lowercase_long_names(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'long_name'):
            var.long_name = var.long_name.lower();
            if "utc" in var.long_name:
                var.long_name = var.long_name.replace("utc", "UTC")


    DS.close()

def time_long_name(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');
    time_var = DS.variables['time'];
    if 'time_reference' in DS.variables:
        time_var.long_name = "time in seconds since time_reference"

    DS.close()

def process_kepler(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,tracking_tag):

    pattern = '*{}*.mmclx'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    mmclxfiles = [];
    mmclxdirs = [];

    vertfiles = [];

    for root,dirs,files in os.walk(datepath):
        mmclxfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        mmclxdirs += dirs;
        

    data_version = "1.0.0";

    l1path = os.path.join(outpath,'L1',datestr);

    os.makedirs(l1path,exist_ok=True);

    for dir in mmclxdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l1path,dir),exist_ok=True);

    for f in mmclxfiles:
        convert_kepler_mmclx2l1(f,l1path,yaml_project_file,yaml_instrument_file,tracking_tag);

    return

def find_mmclxfiles(start_time, end_time, sweep_type,inpath,gzip_flag=False):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')

    # Define a list to store the found files
    matching_files = []

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range
            if gzip_flag:
                if sweep_type in file and file.endswith('.mmclx.gz'):
                    fullfile = os.path.join(root,file)
                    with gzip.open(fullfile) as gz:
                        with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                            file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                            if start_datetime <= file_time <= end_datetime:
                                matching_files.append(os.path.join(root, file))
            else:
                if sweep_type in file and file.endswith('.mmclx'):
                    nc_file = nc4.Dataset(os.path.join(root, file))
                    file_time = cftime.num2pydate(nc_file['time'][0],'seconds since 1970-01-01 00:00:00')
                    nc_file.close()      
                    if start_datetime <= file_time <= end_datetime:
                        matching_files.append(os.path.join(root, file))
        print(sorted(matching_files));

    return sorted(matching_files)

def convert_angle(angle,offset):
    print(angle,offset)
    if angle >= 360+round(offset):    
        angle -= 360
    return angle
    
def find_mmclx_rhi_files(start_time, end_time,azim_min,azim_max,inpath,gzip_flag=False,azimuth_offset=-6.85,revised_northangle=302.15):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S') 
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    print(start_datetime);
    print(end_datetime);
    hrstr = start_datetime.strftime("%H");
    # Define a list to store the found files
    matching_files = []

    #az_search_offset = -38.0; # July ? onwards
    az_search_offset = -8.0; # Before ? July


    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range
            if gzip_flag:
                #if "rhi" in file and hrstr in file and file.endswith('.mmclx.gz'):
                if "rhi" in file and file.endswith('.mmclx.gz'):
                    fullfile = os.path.join(root,file)
                    with gzip.open(fullfile) as gz:
                        with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                            if len(nc.dimensions['time'])>0:
                                file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                                #azim = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) % 360;
                                azim = (nc['azi'][5]+revised_northangle) % 360;
                                #print(f'{file} ST:{start_datetime} {file_time} FI:{end_datetime}')
                                if start_datetime <= file_time <= end_datetime:
                                    #if azim_min <= convert_angle(azim,azimuth_offset) < azim_max:
                                    if azim_min <= convert_angle(azim,az_search_offset) < azim_max:
                                        print(f'{file}: {file_time} {azim_min} {convert_angle(azim,az_search_offset)} {azim_max}');
                                        matching_files.append(os.path.join(root, file))
            else:
                #if "rhi" in file and hrstr in file and file.endswith('.mmclx'):
                if "rhi" in file and file.endswith('.mmclx'):
                    nc = nc4.Dataset(os.path.join(root, file))
                    if len(nc.dimensions['time'])>0:
                        file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                        #azim = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) % 360;
                        azim = (nc['azi'][5]+revised_northangle) % 360;

                        if start_datetime <= file_time <= end_datetime:
                            if azim_min <= convert_angle(azim,az_search_offset) < azim_max:
                                print(f'{file}: {file_time} {azim_min} {convert_angle(azim,az_search_offset)} {azim_max}');
                                #print(f'{file_time} {convert_angle(azim)}');
                                matching_files.append(os.path.join(root, file))
                    nc.close()         
    return sorted(matching_files)
    

def find_mmclx_ppi_files(start_time, end_time,elev_min,elev_max,inpath,gzip_flag=False):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    hrstr = start_datetime.strftime("%H");

    # Define a list to store the found files
    matching_files = []

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            if gzip_flag:
                #if "ppi" in file and hrstr in file and file.endswith('.mmclx.gz'):
                if "ppi" in file and file.endswith('.mmclx.gz'):
                    fullfile = os.path.join(root,file)
                    with gzip.open(fullfile) as gz:
                        with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                            file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                            elev = nc['elv'][0];
                            if start_datetime <= file_time <= end_datetime:
                                print(f'{file_time} {elev_min} {elev} {elev_max}');
                                if elev_min <= elev <= elev_max:
                                    print('Match')
                                    matching_files.append(fullfile);
            else:
                #if "ppi" in file and hrstr in file and file.endswith('.mmclx'):
                if "ppi" in file and file.endswith('.mmclx'):
                    fullfile = os.path.join(root, file);
                    nc = nc4.Dataset(fullfile)
                    file_time = cftime.num2pydate(nc_file['time'][0],'seconds since 1970-01-01 00:00:00')
                    elev = nc['elv'][0];
                    if start_datetime <= file_time <= end_datetime:
                        print(f'{file_time} {elev_min} {elev} {elev_max}');
                        if elev_min <= elev <= elev_max:
                            matching_files.append(fullfile);
                    nc.close()   
            print(sorted(matching_files));
    return sorted(matching_files)

def find_mmclx_vad_files(start_time, end_time,elev_min,elev_max,inpath,gzip_flag=False):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')

    # Define a list to store the found files
    matching_files = []

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            if gzip_flag:
                if "ppi" in file and file.endswith('.mmclx.gz'):
                    fullfile = os.path.join(root,file)
                    with gzip.open(fullfile) as gz:
                        with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                            file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                            elev = nc['elv'][0];
                            if start_datetime <= file_time <= end_datetime:
                                print(file_time);
                                if elev_min <= elev <= elev_max:
                                    matching_files.append(os.path.join(root, file))
            else:
                if "ppi" in file and file.endswith('.mmclx'):
                    nc = nc4.Dataset(os.path.join(root, file))
                    file_time = cftime.num2pydate(nc_file['time'][0],'seconds since 1970-01-01 00:00:00')
                    elev = nc['elv'][0];
                    if start_datetime <= file_time <= end_datetime:
                        print(file_time);
                        if elev_min <= elev <= elev_max:
                            matching_files.append(os.path.join(root, file))
                    nc.close()         
    return sorted(matching_files)

def multi_mmclx2cfrad(
    mmclxfiles,
    output_dir,
    scan_name="HSRHI",
    gzip_flag=False,
    azimuth_offset=0.0,
    #revised_northangle=-56.85,
    revised_northangle=55.7,
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version="1.0.0",
):
    """
    Aggregates single-sweep mmclx data to a cfradial1 data.
    output_dir(str): Enter the path for output data,
    scan_name(str): "HSRHI"
    """
    #pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    #from kepler_utils import read_mira35_mmclx
    
    from pathlib import Path
    homepath = Path.home()

    yaml_project_file = os.path.join(homepath,'amof_campaigns','{}_project.yml'.format(campaign))
    yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

    out_dir = output_dir
    files = sorted(mmclxfiles)
    print(files)
    print("Number of files: ", len(files))
    print(f"gzip_flag={gzip_flag}");

    print('Start to read mmclx')

    RadarDS = read_mira35_mmclx(files[0],gzip_flag=gzip_flag,revised_northangle=revised_northangle);

    print('Done reading mmclx')
    # Read time and microsec directly from mmclx file
    if gzip_flag:
        with gzip.open(files[0]) as gz:
            with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                dtsec = cftime.num2pydate(nc['time'][:],'seconds since 1970-01-01 00:00:00')
                usec = nc['microsec'][:];           
    else:
        nc = nc4.Dataset(files[0])
        dtsec = cftime.num2pydate(nc['time'][:],'seconds since 1970-01-01 00:00:00')
        usec = nc['microsec'][:]; 
        nc.close();

    dt_ref = dtsec[0].replace(hour=0,minute=0,second=0,microsecond=0);
    time_reference = dt_ref.strftime('%Y-%m-%dT%H:%M:%SZ');
    time_units = f'seconds since {time_reference}'; 

    tsec = cftime.date2num(dtsec,time_units);
    print(time_units);
    
    print(RadarDS.latitude['data']);
    print("Merging all scans into one Volume")
    for i in range(1, len(files)):

        newRadarDS = read_mira35_mmclx(files[i],gzip_flag=gzip_flag)

        if gzip_flag:
            with gzip.open(files[i]) as gz:
                with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                    dtsec_new = cftime.num2pydate(nc['time'][:],'seconds since 1970-01-01 00:00:00')
                    usec_new = nc['microsec'][:];           
        else:
            nc = nc4.Dataset(files[i])
            dtsec_new = cftime.num2pydate(nc['time'][:],'seconds since 1970-01-01 00:00:00')
            usec_new = nc['microsec'][:]; 
            nc.close();

        if 'RHI' in scan_name or 'rhi' in scan_name:
            if np.max(newRadarDS.elevation['data'])-np.min(newRadarDS.elevation['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
                tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
                usec = np.append(usec,usec_new);

                RadarDS.scan_rate['data'] = np.append(RadarDS.scan_rate['data'],newRadarDS.scan_rate['data']);
                RadarDS.antenna_transition['data'] = np.append(RadarDS.antenna_transition['data'],newRadarDS.antenna_transition['data']);
        elif 'PPI' in scan_name or 'ppi' in scan_name or 'VAD' in scan_name or 'vad' in scan_name:
            if np.max(newRadarDS.azimuth['data'])-np.min(newRadarDS.azimuth['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
                tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
                usec = np.append(usec,usec_new);

                RadarDS.scan_rate['data'] = np.append(RadarDS.scan_rate['data'],newRadarDS.scan_rate['data']);
                RadarDS.antenna_transition['data'] = np.append(RadarDS.antenna_transition['data'],newRadarDS.antenna_transition['data']);
        elif 'VPT' in scan_name or 'vpt' in scan_name or 'VERT' in scan_name or 'vert' in scan_name:        
            RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
            tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
            usec = np.append(usec,usec_new);

            #RadarDS.scan_rate['data'] = None;
            #RadarDS.antenna_transition['data'] = None;
        
 
        
    RadarDS.time['units'] = time_units;
    RadarDS.time['data'][:] = tsec+usec*1e-6;



    #RadarDS.azimuth['data'] += azimuth_offset;


    #if 'RHI' in scan_name or 'rhi' in scan_name:
    #    RadarDS.fixed_angle['data'] += azimuth_offset;
    
    fname = os.path.basename(files[0]).split(".")[0]

    out_file = f"{fname}_{scan_name.lower()}.nc"

    print(out_file);
    out_path = os.path.join(out_dir, out_file)
    print(out_path);
    pyart.io.write_cfradial(out_path, RadarDS, format="NETCDF4")

    DS = nc4.Dataset(out_path,'r+');
    DS.scan_name = scan_name.lower();

    DS.close();

    amend_unitless(out_path)
    lowercase_long_names(out_path)
    time_long_name(out_path)

    # Update history
    update_string = 'Merge single sweep files into cfradial file'
    update_history_attribute(out_path,update_string)

    print(mmclxfiles[0]);
    cfradial_add_instrument_parameters(mmclxfiles[0],out_path,yaml_project_file,yaml_instrument_file,tracking_tag, data_version)

    cfradial_add_geometry_correction(out_path,revised_northangle);


    print(out_path);
    print('Going to add NCAS metadata')
    outfile = cfradial_add_ncas_metadata(out_path,yaml_project_file,yaml_instrument_file,tracking_tag,data_version);

 
    return outfile

def ppistack_mmclx2cfrad(
    mmclxfiles,
    output_dir,
    scan_name="BLPPI",
):
    """
    Aggregates single-sweep PPI data to a cfradial1 data.
    input_dir(str): Enter path of single sweep data directory,
    output_dir(str): Enter the path for output data,
    scan_name(str): "BLPPI"
    """
    #pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    in_dir = input_dir
    out_dir = output_dir
    files = sorted(mmclxfiles, key=natural_sort_key)
    print("Number of files: ", len(files))
 
    RadarDS = read_mira35_mmclx(files[0]);
    
    print("Merging all scans into one Volume")
    for i in range(1, len(files)):

        newRadarDS = read_mira35_mmclx(files[i])
        RadarDS = pyart.util.join(RadarDS,newRadarDS)

        fname = os.path.basename(files[0]).split(".nc")[0]

        out_file = f"cfrad_{fname}.nc"
        out_path = os.path.join(out_dir, out_file)
        pyart.io.write_cfradial(out_path, RadarDS, format="NETCDF4")

def split_monotonic_sequence(sequence):
    subsequences = [];
    increasing_subsequence = []
    for i in range(len(sequence)):
        if i == 0:
            increasing_subsequence.append(sequence[i])
            continue
        if np.round(10*sequence[i],decimals=1) > np.round(10*sequence[i-1],decimals=1):
            increasing_subsequence.append(sequence[i])
        else:
            if increasing_subsequence:
                subsequences.append(increasing_subsequence);
                increasing_subsequence = [sequence[i]]

    if increasing_subsequence:
        subsequences.append(increasing_subsequence);
    endindices = np.cumsum([len(s) for s in subsequences])-1;
    startindices =  [1+endindices[i] - len(subsequences[i]) for i in range(len(endindices))];
    return list(zip(startindices,endindices))

def read_woest_cmd_file(file_path):
    all_data = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.rstrip().replace(': ', ':')
            pairs = line.split()
            line_data = {}

            for pair in pairs:
                key, value = pair.split(':')
                line_data[key.strip()] = value.strip()

            all_data.append(line_data)

    return all_data

def process_kepler_woest_iop_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,azimuth_offset=-6.85,gzip_flag=True,revised_northangle=-56.85,data_version="1.0.0"):

    ioppath = os.path.join(indir,'iop');

    cmd_path = "/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data/campaign/woest/woest-command/";
    cmd_file = os.path.join(cmd_path,f'combined_logs_{datestr}.txt');
    cmd_dicts = read_woest_cmd_file(cmd_file);

    for d in cmd_dicts:
        print(d);
        storm_az = float(d["azimuth"]);
        storm_range = float(d["range"]);
        woest_cmd = d["command"];

        print(storm_az,storm_range,woest_cmd);

    from itertools import groupby

    # Sort the data by the 'region' key
    #cmd_dicts.sort(key=lambda x: x['region'])
    cmd_dicts.sort(key=lambda x: x['Timestamp'])

    # Group the data by the 'region' key
    grouped_data = {key: list(group) for key, group in groupby(cmd_dicts, key=lambda x: x['region'])}

    # Print the grouped data
    for region, values in grouped_data.items():
        print(f"Region: {region}")
        for value in values:
            print(f"  {value}")

    all_regions = [];
    for region, values in grouped_data.items():
        region_data = {};
        time_start = grouped_data[region][0]["Timestamp"];
        time_end = grouped_data[region][-1]["Timestamp"];
        print(region,time_start,time_end);
        region_data["region"] = region;
        region_data["time_start"] = time_start;
        region_data["time_end"] = time_end;
        all_regions.append(region_data);


    print(all_regions);
    time_start = cmd_dicts[0]["Timestamp"];
    time_end = cmd_dicts[-1]["Timestamp"];

    print(time_start);
    print(time_end);    

    # Specify the directory containing the gzipped NetCDF files
    #directory = f'/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-mobile-ka-band-1/data/campaign/woest/mom/{datestr}/iop'

    # Get all files in the specified directory
    gzipped_files = [os.path.join(ioppath, file) for file in os.listdir(ioppath) if file.endswith('mmclx.gz')]

    files_by_region = {region: [] for region,value in grouped_data.items()}
    #rhi1_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), -10, 169, hsrhipath,gzip_flag=True, azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)

    for gzipped_file in gzipped_files:
        with gzip.open(gzipped_file, 'rb') as f:
            with nc4.Dataset('dummy',mode='r',memory=f.read()) as DS:
                time_var = DS.variables['time'];  

                region_id = None;
                for region in all_regions:
                    region_start = datetime.datetime.strptime(region['time_start'], '%y%m%d-%H%M%S')
                    region_end = datetime.datetime.strptime(region['time_end'], '%y%m%d-%H%M%S')

                    time_values = time_var[:]
                    time_units = time_var.units
             
                    dtime = cftime.num2pydate(time_values,'seconds since 1970-01-01 00:00:00');
                
                    if region_start <= dtime[0] <= region_end and region_start <= dtime[-1] <= region_end:
                        print(region_start,region_end)
                        print(dtime[0],dtime[-1]);
                        region_id = region["region"]
                        print(region_id)

                    if region_id:
                        files_by_region[region_id].append(gzipped_file)
                        #break

    files_by_region = {region: files for region, files in files_by_region.items() if files}

    for region, file_list in files_by_region.items():
        file_set = set(file_list);
        iop_files = list(file_set)
        print(f"Region {region}: {file_list}")

        RadarDS_IOP = multi_mmclx2cfrad(iop_files,outdir,scan_name='SRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,data_version=data_version);

        DS = nc4.Dataset(RadarDS_IOP,'r+');
        DS.comment = f"WOEST IOP tracking storm region {region}";
        DS.close();

def process_kepler_woest_day_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,azimuth_offset=-6.85,gzip_flag=True,revised_northangle=-56.85,data_version="1.0.0"):
#def process_kepler_woest_day_step1(datestr,indir,outdir,azimuth_offset,gzip_flag=True,revised_northangle=303.15):
    # Define the start and end times for the loop


    start_date = datetime.datetime.strptime(datestr, '%Y%m%d'); # + datetime.timedelta(hours=1);
    end_date = start_date  + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);

    hsrhipath = os.path.join(indir,'hsrhi');
    blppipath = os.path.join(indir,'blppi');
    vadpath = os.path.join(indir,'vad');
    vptpath = os.path.join(indir,'vpt');
    ioppath = os.path.join(indir,'iop');
    ppipath = os.path.join(indir,'ppi');

    # Correct azimuth offset based on use of revised northangle
    

    # Iterate through each half hour for HSRHI and BLPPI files
    current_date = start_date
    while current_date <= end_date:
        print(current_date);
        next_halfhour = current_date + datetime.timedelta(minutes=30);
        
        try:
            print("searching for hsrhi1 files")
            hsrhi1_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), -10, 169, hsrhipath,gzip_flag=True, azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)
            #hsrhi1_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), -40, 139, indir,gzip_flag=True, azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)

            print(hsrhi1_files);
            list_length = len(hsrhi1_files);
            for i in range(0, list_length, 6):
                subset = hsrhi1_files[i:i+6]
                RadarDS_HSRHI1 = multi_mmclx2cfrad(subset,outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,data_version=data_version);

            remaining_files = hsrhi1_files[i+6:]                
            RadarDS_HSRHI1 = multi_mmclx2cfrad(remaining_files,outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,data_version=data_version);

            #azims = [];
            #if (len(hsrhi1_files)>0):
                #hsrhi1_files.sort();
            #    if gzip_flag:
            #        for f in hsrhi1_files:
            #            print(f);
            #            with gzip.open(f) as gz:
            #                with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
            #                    nc.set_auto_mask(False);
                                #az = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) %360;
            #                    az = (nc['azi'][0]+revised_northangle) %360;
            #                    azims.append(az);
            #    else:
            #        for f in hsrhi1_files:
            #            nc = nc4.Dataset(f);
            #            nc.set_auto_mask(False);
                        #az = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) %360;
            #            az = (nc['azi'][0]+revised_northangle) %360;
            #            azims.append(az);
            #            nc.close();
            #    print(azims);
            #    idx = split_monotonic_sequence(azims);
            #    print(idx);
            #    for l in idx:
                    #if len(hsrhi1_files)==1:
                    #    RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files,outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
            #        if l[0]==l[1]:
            #            RadarDS_HSRHI1 = multi_mmclx2cfrad([hsrhi1_files[l[0]]],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);
            #        elif l[1]==len(hsrhi1_files)-1:
            #            RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[l[0]:],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);
            #        else:
            #            RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[l[0]:l[1]+1],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);

               # if (len(hsrhi1_files)>6):
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[:6],outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[6:],outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
               # else:
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files,outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
        except:
            print("Problem");
            pass

        try:        
            print("searching for hsrhi2 files")
            hsrhi2_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), 170, 349, hsrhipath,gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)
            #hsrhi2_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), 140, 319, indir,gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)

            print(hsrhi2_files);
            list_length = len(hsrhi2_files);
            for i in range(0, list_length, 6):
                subset = hsrhi2_files[i:i+6]
                RadarDS_HSRHI2 = multi_mmclx2cfrad(subset,outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,data_version=data_version);

            remaining_files = hsrhi2_files[i+6:]                
            RadarDS_HSRHI2 = multi_mmclx2cfrad(remaining_files,outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,data_version=data_version);

            #azims = [];
            #if (len(hsrhi2_files)>0):
            #    hsrhi2_files.sort();
            #    if gzip_flag:
            #        for f in hsrhi2_files:
            #            with gzip.open(f) as gz:
            #                with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
            #                    nc.set_auto_mask(False);
                                #az = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) %360;
            #                    az = (nc['azi'][0]+revised_northangle) %360;
            #                    azims.append(az);
            #    else:
            #        for f in hsrhi2_files:
            #            nc = nc4.Dataset(f);
            #            nc.set_auto_mask(False);
                        #az = (nc['azi'][0]+nc['northangle'][0]+azimuth_offset) %360;
            #            az = (nc['azi'][0]+revised_northangle) %360;
            #            azims.append(az);
            #            nc.close();
            #    print(azims);
            #    idx = split_monotonic_sequence(azims);
            #    print(idx);
            #    for l in idx:
                    #if len(hsrhi2_files)==1:
                    #    RadarDS_HSRHI2 = multi_mmclx2cfrad(hsrhi2_files,outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
            #        if l[0]==l[1]:
            #            RadarDS_HSRHI2 = multi_mmclx2cfrad([hsrhi2_files[l[0]]],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
            #        elif l[1]==len(hsrhi2_files)-1:
            #            RadarDS_HSRHI2 = multi_mmclx2cfrad(hsrhi2_files[l[0]:],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
            #        else:
            #            RadarDS_HSRHI2 = multi_mmclx2cfrad(hsrhi2_files[l[0]:l[1]+1],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
        except:
            print('Problem');
            pass
        
        try:
            print("searching for blppi files")
            blppi_files = find_mmclx_ppi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), 0, 80, blppipath,gzip_flag=gzip_flag)
            print('FOUND BLPPI FILES')
            print(blppi_files);
            RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files,outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle,data_version=data_version);

            #elevs=[];
            #if (len(blppi_files)>0):
            #    blppi_files.sort();
            #    if gzip_flag:
            #        for f in blppi_files:
            #            with gzip.open(f) as gz:
            #                with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
            #                    nc.set_auto_mask(False);
            #                    el = nc['elv'][0];
            #                    elevs.append(el);
            #    else:
            #        for f in blppi_files:
            #            nc = nc4.Dataset(f,'r');
            #            nc.set_auto_mask(False);
            #            el = nc['elv'][0];
            #            elevs.append(el);
            #            nc.close();

            #    print(elevs);
            #    idx = split_monotonic_sequence(elevs);
            #    print(idx);
            #    for l in idx:
                    #if len(blppi_files)==1:
                    #    RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files,outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset);
            #        if l[0]==l[1]:
            #            RadarDS_BLPPI = multi_mmclx2cfrad([blppi_files[l[0]]],outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
            #        elif l[1]==len(blppi_files)-1:   
            #            RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files[l[0]:],outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
            #        else:
            #            RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files[l[0]:l[1]+1],outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle);
        except:
            print('Problem with BLPPI processing')
            pass

        current_date = next_halfhour


    # Individual PPI files for whole day
    try:
        ppi_files = find_mmclx_ppi_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'), 0, 80, ppipath,gzip_flag=gzip_flag)
        print(ppi_files);
        if (len(ppi_files)>0):
            ppi_files.sort();
            for f in ppi_files:
                RadarDS_PPI = multi_mmclx2cfrad([f],outdir,scan_name='PPI',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle,data_version=data_version);
    except:
        pass
    
    # Vertically pointing files for whole day
    try:
        vpt_files = find_mmclxfiles(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', vptpath,gzip_flag=True);
        print(vpt_files);
        if (len(vpt_files)>0):
            vpt_files.sort();
            RadarDS_VPT = multi_mmclx2cfrad(vpt_files,outdir,scan_name='VPT',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle,data_version=data_version);
    except:
        pass
    # VAD files for whole day
    vad_dt_start = datetime.datetime.strptime(datestr,"%Y%m%d");
    vad_dt_end = vad_dt_start+datetime.timedelta(days=1);
    try:
        vad_files = find_mmclx_vad_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),80,90, vadpath,gzip_flag=True);
        if (len(vad_files)>0):
            vad_files.sort();
            RadarDS_VAD = multi_mmclx2cfrad(vad_files,outdir,scan_name='VAD',gzip_flag=True,azimuth_offset=azimuth_offset,revised_northangle=revised_northangle,data_version=data_version);
    except:
        pass



def process_kepler_ccrest_day_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,azimuth_offset=0.0,gzip_flag=True,revised_northangle=55.7,data_version="1.0.0"):

    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    end_date = start_date + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);

    # Correct azimuth offset based on use of revised northangle
    azimuth_offset = 0.0;

    rhi_files_270 = find_mmclx_rhi_files(start_date.strftime('%Y-%m-%d %H:%M:%S'), end_date.strftime('%Y-%m-%d %H:%M:%S'), 260, 280, indir,gzip_flag=True, azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)
    rhi_files_246 = find_mmclx_rhi_files(start_date.strftime('%Y-%m-%d %H:%M:%S'), end_date.strftime('%Y-%m-%d %H:%M:%S'), 236, 256, indir,gzip_flag=True, azimuth_offset=azimuth_offset,revised_northangle=revised_northangle)

    
    print(rhi_files_270);
    print(rhi_files_246);
    print(len(rhi_files_270));
    print(len(rhi_files_246));

    if (len(rhi_files_270)>0):
        if gzip_flag:
            RadarDS_RHI = multi_mmclx2cfrad(rhi_files_270,outdir,scan_name='RHI-CCREST1',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
        else:
            RadarDS_RHI = multi_mmclx2cfrad(rhi_files_270,outdir,scan_name='RHI-CCREST1',gzip_flag=False,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);


    if (len(rhi_files_246)>0):
        if gzip_flag:
            RadarDS_RHI = multi_mmclx2cfrad(rhi_files_246,outdir,scan_name='RHI-CCREST2',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
        else:
            RadarDS_RHI = multi_mmclx2cfrad(rhi_files_246,outdir,scan_name='RHI-CCREST2',gzip_flag=False,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);

    
    # Vertically pointing files for whole day
    try:
        vpt_files = find_mmclxfiles(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', indir,gzip_flag=True);
        print(vpt_files)
        #vpt_files_unzipped = find_mmclxfiles(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', indir,gzip_flag=False);
        if (len(vpt_files)>0):
            RadarDS_VPT = multi_mmclx2cfrad(vpt_files,outdir,scan_name='VPT',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
        #elif (len(vpt_files_unzipped)>0):
        #    RadarDS_VPT = multi_mmclx2cfrad(vpt_files_unzipped,outdir,scan_name='VPT',gzip_flag=False,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=55.7);
    except:
        print("VPT problem")
        pass
    # VAD files for whole day
    vad_dt_start = datetime.datetime.strptime(datestr,"%Y%m%d");
    vad_dt_end = vad_dt_start+datetime.timedelta(days=1);
    try:
        vad_files = find_mmclx_vad_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),80,90, indir,gzip_flag=True);
        if (len(vad_files)>0):
            RadarDS_VAD = multi_mmclx2cfrad(vad_files,outdir,scan_name='VAD',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
    except:
        pass

    
def update_history_attribute(ncfile,update):

    dataset = nc4.Dataset(ncfile,'r+');

    user = getpass.getuser();

    updttime = datetime.datetime.utcnow();
    updttimestr = updttime.ctime();

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " " + update);

    dataset.history = history + "\n" + dataset.history;

    dataset.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    dataset.close();



