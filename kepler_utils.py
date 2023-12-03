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



def read_mira35_mmclx(filename, gzip_flag=False, **kwargs):
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

    antenna_transition = None;  #  NEED TO CHECK THIS
     
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
        latitude['data'] = np.array([-z1['lat']]);
    else:
        latitude['data'] = np.array([z1['lat']]);  
        
    z = StringIO(ncobj.getncattr('Longitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lon','zs2'])
    
    if z1['zs2']==b'W' and z1['lon']>0: 
        longitude['data'] = np.array([-z1['lon']]);
    else:
        longitude['data'] = np.array([z1['lon']]);  
    
    z = StringIO(ncobj.getncattr('Altitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['alt','zs2'])
    altitude['data'] = np.array([z1['alt']])

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
        print(key)
        if key in filename.lower(): 
            scan_name = value;
            sweep_mode["data"] = np.array(1 * [value]);
            sweep_mode["data"][0] = value;
            break;

    fixed_angles = {'ppi' : ncvars['elv'][0], 'rhi' : ncvars['azi'][0]+ncvars['northangle'][0], 'vertical_pointing' : ncvars['elv'][0], "manual_rhi" : ncvars['azi'][0]}

    fixed_angle = filemetadata("fixed_angle")

    if scan_name is not None:
        fixed_angle["data"] = np.array(1 * [fixed_angles[scan_name] % 360]) 
    else:
        fixed_angle["data"] = np.array(1 * [None]) 


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
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = (ncvars['azi'][:]+ncvars['northangle'][:]) % 360;
    azimuth['units'] = "degrees";
    azimuth['proposed_standard_name'] = "sensor_to_target_azimuth_angle";
    azimuth['long_name'] = "sensor to target azimuth angle";

    elevation['data'] = ncvars['elv'][:];
    elevation['units'] = "degrees";
    elevation['proposed_standard_name'] = "sensor_to_target_elevation_angle";
    elevation['long_name'] = "sensor to target elevation angle";


    metadata['time_coverage_start'] = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
    metadata['time_coverage_end'] = datetime.datetime.strftime(dtime[-1],'%Y-%m-%dT%H:%M:%SZ');

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
        field_dic['long_name'] =  "radar equivalent reflectivity factor";
        field_dic['standard_name'] = "equivalent_reflectivity_factor";
        field_dic['proposed_standard_name'] =  "radar_equivalent_reflectivity_factor";   
        fields[field_name] = field_dic
    else:
        print("Zg does not exist")

    if "VELg" in ncvars:
        field_name = fields_keymap['VELg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['VELg'][:];
        field_dic['long_name'] =  "radial velocity of scatterers away from instrument";
        field_dic['standard_name'] = "radial_velocity_of_scatterers_away_from_instrument";
        fields[field_name] = field_dic
    else:
        print("VELg does not exist")

    if "RMSg" in ncvars:
        field_name = fields_keymap['RMSg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['RMSg'][:];
        field_dic['long_name'] =  "radar doppler spectrum width";
        field_dic['proposed_standard_name'] = "radar_doppler_spectrum_width";
        fields[field_name] = field_dic
    else:
        print("RMSg does not exist")

    if "LDRg" in ncvars:
        field_name = fields_keymap['LDRg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = 10.0*np.log10(ncvars['LDRg'][:]);
        field_dic['long_name'] =  "radar linear depolarization ratio";
        field_dic['proposed_standard_name'] = "radar_linear_depolarization_ratio";
        fields[field_name] = field_dic
    else:
        print("LDRg does not exist")

    if "SNRg" in ncvars:
        field_name = fields_keymap['SNRg']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = 10.0*np.log10(ncvars['SNRg'][:]);
        field_dic['long_name'] =  "radar signal to noise ratio";
        field_dic['proposed_standard_name'] = "radar_signal_to_noise_ratio";
        fields[field_name] = field_dic
    else:
        print("SNRg does not exist")
    

    # instrument_parameters
    instrument_parameters = {}

    radar_calibration = {}


    #if "PRF-value" in dset.ncattrs():
    #    dic = filemetadata("prt")
    #    prt = 1.0 / float(dset.getncattr("PRF-value"))
    #    dic["data"] = np.ones((nrays,), dtype="float32") * prt
    #    instrument_parameters["prt"] = dic

    #if "PulseWidth-value" in dset.ncattrs():
    #    dic = filemetadata("pulse_width")
    #    pulse_width = dset.getncattr("PulseWidth-value") * 1.0e-6
    #    dic["data"] = np.ones((nrays,), dtype="float32") * pulse_width
    #    instrument_parameters["pulse_width"] = dic

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

    return Radar(
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
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration
    )



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

def cfradial_add_ncas_metadata(cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    # -------------------------------------------------------
    # Read cfradial file to add NCAS metadata
    # -------------------------------------------------------
    DS = nc4.Dataset(cfradfile,'r+');

    DS.product_version = "v{}".format(data_version) ;
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

    DSin = nc4.Dataset(mmclxfile,'r');

    frequency = DS.createDimension("frequency", 1);

    varin = DSin['frequency'];
    varout = DSout.createVariable('frequency',varin.datatype,("frequency"));
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=varin[:];
    varout.meta_group = "instrument_parameters";

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

    time_coverage_start = str(nc4.chartostring(RadarDataset.variables['time_coverage_start'][:]));
    time_coverage_end = str(nc4.chartostring(RadarDataset.variables['time_coverage_end'][:]));


    file_timestamp = datetime.datetime.strptime(time_coverage_start,'%Y-%m-%dT%H:%M:%SZ');

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

    DS.Conventions = "NCAS-Radar-1.0 CfRadial-1.4"
    DS.delncattr('version');
    DS.product_version = "v{}".format(data_version) ;
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
    DS.geospatial_bounds = "";

    if "vpt" in DS.scan_name or "VPT" in DS.scan_name or "vertical_pointing" in DS.scan_name:
        DS.featureType = 'timeSeriesProfile';

    # -------------------------------------------------------
    # Now clean up some variable attributes
    # -------------------------------------------------------
    DS['time'].comment = "";
    DS['range'].delncattr('standard_name');
    DS['range'].comment = 'Range to centre of each bin';
    DS['range'].meters_to_center_of_first_gate = DS['range'][0];
    DS['azimuth'].delncattr('standard_name');
    DS['elevation'].delncattr('standard_name');
    DS['DBZ'].standard_name = 'equivalent_reflectivity_factor';
    DS['sweep_number'].delncattr('standard_name');
    DS['sweep_mode'].delncattr('standard_name');
    DS['fixed_angle'].delncattr('standard_name');
    DS['latitude'].long_name = 'latitude';
    DS['latitude'].standard_name = 'latitude';
    DS['longitude'].long_name = 'longitude';
    DS['longitude'].standard_name = 'longitude'; 
    DS['altitude'].standard_name = 'altitude';
    DS['altitude'].comment = 'Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model' 
    DS['altitude'].long_name = 'altitude';
    DS['altitude'].units = 'metres';
    DS['altitude'].delncattr('positive'); 
    DS['volume_number'].long_name = 'volume number';


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

    return

def anglicise_cfradial(cfradfile):

    # -------------------------------------------------------
    # Read cfradial file and change language localisation
    # -------------------------------------------------------

    DS = nc4.Dataset(cfradfile,'r+');

    variable = nc_file.variables['my_variable']
    variable.setncattr('attribute_name', 'new_attribute_value')

    DS.close();

    return


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
            #print(file);
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

    return sorted(matching_files)

def find_mmclx_rhi_files(start_time, end_time,azim_min,azim_max,inpath,gzip_flag=False):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    hrstr = start_datetime.strftime("%H");
    # Define a list to store the found files
    matching_files = []

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range
            if gzip_flag:
                if "rhi" in file and hrstr in file and file.endswith('.mmclx.gz'):
                    fullfile = os.path.join(root,file)
                    with gzip.open(fullfile) as gz:
                        with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                            file_time = cftime.num2pydate(nc['time'][0],'seconds since 1970-01-01 00:00:00')
                            azim = (nc['azi'][0]+nc['northangle'][0]) % 360;
                            if start_datetime <= file_time <= end_datetime:
                                print(file_time);
                                if azim_min < azim <= azim_max:
                                    matching_files.append(os.path.join(root, file))
            else:
                if "rhi" in file and hrstr in file and file.endswith('.mmclx'):
                    nc = nc4.Dataset(os.path.join(root, file))
                    file_time = cftime.num2pydate(nc_file['time'][0],'seconds since 1970-01-01 00:00:00')
                    azim = (nc['azi'][0]+nc['northangle'][0]) % 360;
                    if start_datetime <= file_time <= end_datetime:
                        print(file_time);
                        if azim_min < azim <= azim_max:
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
                if "ppi" in file and hrstr in file and file.endswith('.mmclx.gz'):
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
                if "ppi" in file and hrstr in file and file.endswith('.mmclx'):
                    nc = nc4.Dataset(os.path.join(root, file))
                    file_time = cftime.num2pydate(nc_file['time'][0],'seconds since 1970-01-01 00:00:00')
                    elev = nc['elv'][0];
                    if start_datetime <= file_time <= end_datetime:
                        print(file_time);
                        if elev_min <= elev <= elev_max:
                            matching_files.append(os.path.join(root, file))
                    nc.close()         
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
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version=0.1,
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
    print("Number of files: ", len(files))
    print(f"gzip_flag={gzip_flag}");
    RadarDS = read_mira35_mmclx(files[0],gzip_flag=gzip_flag);

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
        print(len(newRadarDS.latitude["data"]))

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

        tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
        usec = np.append(usec,usec_new);

        RadarDS = pyart.util.join_radar(RadarDS,newRadarDS)

    RadarDS.time['units'] = time_units;
    RadarDS.time['data'][:] = tsec+usec*1e-6;

    RadarDS.azimuth['data'] += azimuth_offset;


    if 'RHI' in scan_name or 'rhi' in scan_name:
        RadarDS.fixed_angle['data'] += azimuth_offset;
    
    fname = os.path.basename(files[0]).split(".")[0]

    out_file = f"{fname}_{scan_name.lower()}.nc"

    print(out_file);
    out_path = os.path.join(out_dir, out_file)
    print(out_path);
    pyart.io.write_cfradial(out_path, RadarDS, format="NETCDF4")

    DS = nc4.Dataset(out_path,'r+');
    DS.scan_name = scan_name.lower();

    DS.close();

    # Update history
    update_string = 'Merge single sweep files into cfradial file'
    update_history_attribute(out_path,update_string)

    cfradial_add_ncas_metadata(out_path,yaml_project_file,yaml_instrument_file,tracking_tag,data_version);
    
    return

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
    endindices = np.cumsum([len(s) for s in subseq])-1;
    print(endindices);
    startindices =  endindex+1 - [len(s) for s in subseq];
    print(startindices);
    return list(zip(startindices,endindices))

def process_kepler_woest_day_step1(datestr,indir,outdir,azimuth_offset,gzip_flag=True):
    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    end_date = start_date + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);

    # Iterate through each half hour for HSRHI and BLPPI files
    current_date = start_date
    while current_date <= end_date:
        print(current_date);
        next_halfhour = current_date + datetime.timedelta(minutes=30);
        
        try:
            hsrhi1_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), -15, 165, indir,gzip_flag=True)
            
            if (len(hsrhi1_files)>0):
                azims = [];
                if gzip_flag:
                    for f in hsrhi1_files:
                        with gzip.open(f) as gz:
                            with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                                nc.set_auto_mask(False);
                                az = nc['azi'][0];
                                azims.append(az);
                else:
                    for f in hsrhi1_files:
                        nc = nc4.Dataset(f);
                        nc.set_auto_mask(False);
                        az = nc['azi'][0];
                        azims.append(az);
                        nc.close();
                print(azims);
                idx = split_monotonic_sequence(azims);
                print(idx);
                for l in idx:
                    if l[0]==l[1]:
                        RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[l[0]],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);
                    else:   
                        RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[l[0]:l[1]],outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);

               # if (len(hsrhi1_files)>6):
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[:6],outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files[6:],outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
               # else:
               #     RadarDS_HSRHI1 = multi_mmclx2cfrad(hsrhi1_files,outdir,scan_name='HSRHI',gzip_flag=True,azimuth_offset=azimuth_offset);
        except:
            pass

        try:
            hsrhi2_files = find_mmclx_rhi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), 165, 345, indir, gzip_flag=gzip_flag)
            if (len(hsrhi2_files)>0):
                RadarDS_HSRHI2 = multi_mmclx2cfrad(hsrhi2_files,outdir,scan_name='HSRHI',gzip_flag=gzip_flag,azimuth_offset=azimuth_offset);
        except:
            pass
        
        try:
            blppi_files = find_mmclx_ppi_files(current_date.strftime('%Y-%m-%d %H:%M:%S'), next_halfhour.strftime('%Y-%m-%d %H:%M:%S'), 0, 80, indir,gzip_flag=gzip_flag)
            elevs=[];
            if (len(blppi_files)>0):
                if gzip_flag:
                    for f in blppi_files:
                        with gzip.open(f) as gz:
                            with nc4.Dataset('dummy', mode='r', memory=gz.read()) as nc:
                                nc.set_auto_mask(False);
                                el = nc['elv'][0];
                                elevs.append(el);
                else:
                    for f in blppi_files:
                        nc = nc4.Dataset(f,'r');
                        nc.set_auto_mask(False);
                        el = nc['elv'][0];
                        elevs.append(el);
                        nc.close();

                print(elevs);
                idx = split_monotonic_sequence(elevs);
                print(idx);
                for l in idx:
                    if l[0]==l[1]:
                        RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files[l[0]],outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset);
                    else:   
                        RadarDS_BLPPI = multi_mmclx2cfrad(blppi_files[l[0]:l[1]],outdir,scan_name='BLPPI',gzip_flag=True,azimuth_offset=azimuth_offset);
        except:
            pass
       
        current_date = next_halfhour
    # Vertically pointing files for whole day
    try:
        vpt_files = find_mmclxfiles(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', indir,gzip_flag=True);
        if (len(vpt_files)>0):
            RadarDS_VPT = multi_mmclx2cfrad(vpt_files,outdir,scan_name='VPT',gzip_flag=True,azimuth_offset=azimuth_offset);
    except:
        pass
    # VAD files for whole day
    vad_dt_start = datetime.datetime.strptime(datestr,"%Y%m%d");
    vad_dt_end = vad_dt_start+datetime.timedelta(days=1);
    try:
        vad_files = find_mmclx_vad_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),80,90, indir,gzip_flag=True);
        if (len(vad_files)>0):
            RadarDS_VAD = multi_mmclx2cfrad(vad_files,outdir,scan_name='VAD',gzip_flag=True,azimuth_offset=azimuth_offset);
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
