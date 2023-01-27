#!/usr/bin/env python
# coding: utf-8

# ==========================================================================
# Module for processing mmclx radar files from Kepler (MIRA-35) radar
# Author: Chris Walden, UK Research & Innovation and
#                       National Centre for Atmospheric Science
# Last modified: 27-01-2023
# ==========================================================================

"""Module for processing mmclx radar data from Kepler (MIRA-35) radar."""


import datetime

import netCDF4 as nc4
import numpy as np

from pyart.config import FileMetadata, get_fillvalue
from pyart.core.radar import Radar
from pyart.io.common import _test_arguments, make_time_unit_str

from io import StringIO


def read_mira35_mmclx(filename, **kwargs):
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

    # time, range, fields, metadata, scan_type, latitue, longitude, altitude, altitude_agl,
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index, sweep_end_ray_index, rays_per_sweep,
    # target_scan_rate, rays_are_indexed, ray_angle_res,
    # azimuth, elevation, gate_x, gate_y, gate_z, gate_longitude, gate_latitude, projection, gate_altitude,
    # scan_rate, antenna_transition, 
    # None rotation, tilt, roll, drift, heading, pitch
    # ?? georefs_applied
    # instrument_parameters
    # radar_calibration
    # OK ngates
    # OK nrays
    # OK nsweeps
    
    # The following are not required for a fixed platform
    rotation = None;
    tilt = None;
    roll = None;
    drift = None;
    heading = None;
    pitch = None;
    
    georefs_applied = None;

    antenna_transition = None;

     
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
    ncobj = nc4.Dataset(filename)
    nrays = len(ncobj.dimensions["time"])
    ngates = len(ncobj.dimensions["range"])
    nsweeps = 1;
    
    ncvars = ncobj.variables

    # --------------------------------
    # latitude, longitude and altitude
    # --------------------------------
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    z = StringIO(ncobj.getncattr('Latitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lat','zs2'])
    if z1['zs2']==b'S' and z1['lat']>0: 
        latitude['data'] = -z1['lat']
    else:
        latitude['data'] = z1['lat']  
        
    z = StringIO(ncobj.getncattr('Longitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lon','zs2'])
    
    if z1['zs2']==b'W' and z1['lon']>0: 
        longitude['data'] = -z1['lon']
    else:
        longitude['data'] = z1['lon']  
    
    z = StringIO(ncobj.getncattr('Altitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['alt','zs2'])
    altitude['data'] = z1['alt']


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


    # metadata
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


    # sweep_mode, fixed_angle
    sweep_modes = {'ppi' : 'ppi', 'rhi' : 'rhi','vert' : 'vertical_pointing','man' : 'manual_rhi'}

    sweep_mode = filemetadata("sweep_mode")

    print(filename.lower());

    scan_type = None;
    sweep_mode["data"] = np.array(1 * [None]);

    for key, value in sweep_modes.items():
        print(key)
        if key in filename.lower(): 
            scan_type = value;
            sweep_mode["data"] = np.array(1 * [value]);
            sweep_mode["data"][0] = value;
            break;

    fixed_angles = {'ppi' : ncvars['elv'][0], 'rhi' : ncvars['azi'][0], 'vert' : ncvars['elv'][0], "man" : ncvars['azi'][0]}

    fixed_angle = filemetadata("fixed_angle")

    if scan_type is not None:
        fixed_angle["data"] = np.array(1 * [fixed_angles[scan_type]])
    else:
        fixed_angle["data"] = np.array(1 * [None])

    # time
    # interpolate between the first and last timestamps in the Time variable
    time = filemetadata('time')
    nctime = ncvars['time']
    
    dtime = nc4.num2date(ncvars['time'][:],'seconds since 1970-01-01 00:00:00')

    for idx, x in np.ndenumerate(dtime):
        dtime[idx]=x.replace(microsecond=ncvars['microsec'][idx])
            
    base_time = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    
    time['units'] = make_time_unit_str(base_time)  
    time['data']  = nc4.date2num(dtime,time['units']);

    # range
    _range = filemetadata('range')
    _range['data'] = ncvars['range'][:]
    _range['metres_to_centre_of_first_gate'] = _range['data'][0]
    # assuming the distance between all gates is constant, may not
    # always be true.
    _range['metres_between_gates'] = (_range['data'][1] - _range['data'][0])

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = ncvars['azi'][:]
    elevation['data'] = ncvars['elv'][:]

    # fields
    #field_name = ncobj.TypeName

    #field_data = np.ma.array(dset.variables[field_name][:])
    #if "MissingData" in dset.ncattrs():
    #    field_data[field_data == dset.MissingData] = np.ma.masked
    #if "RangeFolded" in dset.ncattrs():
    #    field_data[field_data == dset.RangeFolded] = np.ma.masked

    fields = {}
    #fields = {field_name: filemetadata(field_name)}
    #fields[field_name]["data"] = field_data
    #fields[field_name]["units"] = dset.variables[field_name].Units
    #fields[field_name]["_FillValue"] = get_fillvalue()


    try:
        ncvars['Zg']
        field_name = filemetadata.get_field_name('Zg')
        field_dic = filemetadata(field_name)
        #field_dic['_FillValue'] = ncvars['Zg']._FillValue
        field_dic['units'] = ncvars['Zg'].units
        field_dic['data'] = 10.0*np.log10(ncvars['Zg'][:]);
        #field_dic['applied_calibration_offset'] = ncvars['ZED_HC'].applied_calibration_offset
        fields[field_name] = field_dic
    except KeyError:
        print("Zg does not exist")





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

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
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
        radar_calibration = radar_calibration
    )





def read_mmclx(filename, **kwargs):
    """
    Read a netCDF mmclx file from MIRA-35 radar.

    Parameters
    ----------
    filename : str
        Name of netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.
    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('mmclx')

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # general parameters
    nrays = ncvars['azi'].shape[0]
    scan_type = 'ppi'

    # time
    # interpolate between the first and last timestamps in the Time variable
    time = filemetadata('time')
    nctime = ncvars['time']
    
    dtime = nc4.num2date(ncvars['time'][:],'seconds since 1970-01-01 00:00:00')

    for idx, x in np.ndenumerate(dtime):
        dtime[idx]=x.replace(microsecond=ncvars['microsec'][idx])
            
    base_time = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    
    
    time['units'] = make_time_unit_str(base_time)  
    time['data']  = date2num(dtime,time['units']);

    # range
    _range = filemetadata('range')
    _range['data'] = ncvars['range'][:]
    _range['metres_to_centre_of_first_gate'] = _range['data'][0]
    # assuming the distance between all gates is constant, may not
    # always be true.
    _range['metres_between_gates'] = (_range['data'][1] - _range['data'][0])

    # fields
    # files contain a single corrected reflectivity field
    fields = {}
    
    try:
        ncvars['Zg']
        field_name = filemetadata.get_field_name('Zg')
        field_dic = filemetadata(field_name)
        #field_dic['_FillValue'] = ncvars['Zg']._FillValue
        field_dic['units'] = 'dBZ';
        field_dic['data'] = 10.0*np.log10(ncvars['Zg'][:]);
        #field_dic['applied_calibration_offset'] = ncvars['ZED_HC'].applied_calibration_offset
        fields[field_name] = field_dic
    except KeyError:
        print("Zg does not exist")
    
    try:
        ncvars['Zg']
        field_name = filemetadata.get_field_name('Zg')
        field_dic = filemetadata(field_name)
        #field_dic['_FillValue'] = ncvars['Zg']._FillValue
        field_dic['units'] = ncvars['Zg'].units
        field_dic['data'] = 10.0*np.log10(ncvars['Zg'][:]);
        #field_dic['applied_calibration_offset'] = ncvars['ZED_HC'].applied_calibration_offset
        fields[field_name] = field_dic
    except KeyError:
        print("Zg does not exist")

    # -----------------------------
    # latitude, longitude, altitude
    # -----------------------------
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    
    # metadata
    metadata = filemetadata('metadata')
    for k in ['institution', 'title', 'used_algorithms']:
        if k in ncobj.ncattrs(): 
            metadata[k] = ncobj.getncattr(k)
    
    from io import StringIO

    z = StringIO(ncobj.getncattr('Latitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lat','zs2'])
    if z1['zs2']==b'S' and z1['lat']>0: 
        latitude['data'] = -z1['lat']
    else:
        latitude['data'] = z1['lat']  
        
    z = StringIO(ncobj.getncattr('Longitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['lon','zs2'])
    
    if z1['zs2']==b'W' and z1['lon']>0: 
        longitude['data'] = -z1['lon']
    else:
        longitude['data'] = z1['lon']  
    
    z = StringIO(ncobj.getncattr('Altitude'))
    
    z1 = np.genfromtxt(z, dtype=None, names=['alt','zs2'])
    altitude['data'] = z1['alt']
    
    
    # sweep parameters
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode   = filemetadata('sweep_mode')
    fixed_angle  = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index   = filemetadata('sweep_end_ray_index')

    # We only store single sweeps 
    sweep_number['data'] = np.arange(1, dtype='int32')
    sweep_mode['data'] =  np.array(1 * ['sector'])  # THIS NEEDS TO BE DETERMINED

    fixed_angle['data'] = np.array([np.round(ncvars['elv'][0],2)], dtype='float32')#np.array([0], dtype='float32')
    sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    sweep_end_ray_index['data'] = np.array([nrays-1], dtype='int32')

    # ------------------
    # azimuth, elevation
    # ------------------
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = ncvars['azi'][:]
    elevation['data'] = ncvars['elv'][:]#np.array([0.], dtype='float32')


    metadata['instrument_name']='ncas-radar-mobile-ka-band-1'

    # ---------------------
    # instrument parameters
    # ---------------------
    instrument_parameters = None
    


    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation, 
        instrument_parameters=instrument_parameters)

# ===================
# CONVERSION ROUTINES
# ===================

def convert_kepler_mmclx2l0b(infile,outfile,yaml_project_file,tracking_tag):

    """This routine converts mmclx data from the NCAS Mobile Ka-band Radar (Kepler) to Level 0b (cfradial) data.
    Metadata are added using information in the yaml_project_file.

    :param infile: Full path of NetCDF Level 0a mmclx data file, e.g. `<path-to-file>/20220907_071502.vert.mmclx`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-mobile-ka-band-1_cao_20220907-071502_fix_l0b_v1.0.nc`
    :type outfile: str
    """

    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

#    project = projects["tracking_tag"];

    radar = "ncas-mobile-radar-ka-band-1";

    for n in project["ncas_instruments"]:
        if radar in n:
            ncas_instrument = n[radar];

    print(ncas_instrument);
            
    DSin = nc4.Dataset(infile);

    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    # FIXME: need to add the microsecs which are stored separately in mmclx files

    # Should modify below to write to cfradial using pyart.
    # Then will modify the files to make them NCAS compliant.

    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')
    print("Creating {}".format(outfile));

    # ------------------------
    # Set up global attributes
    # ------------------------
    DSout.product_version = "v1.0" ;
    DSout.processing_level = "0b" ;

    DSout.licence = project["data_licence"];
    DSout.acknowledgement = project["acknowledgement"];
    DSout.platform = "Chilbolton Atmospheric Observatory" ;
    DSout.platform_type = "stationary_platform" ;
    DSout.title = ncas_instrument["title"];
    DSout.creator_name = "Chris Walden" ;
    DSout.creator_email = "chris.walden@ncas.ac.uk" ;
    DSout.creator_url = "https://orcid.org/0000-0002-5718-466X" ;
    DSout.institution = "National Centre for Atmospheric Science (NCAS)";
    DSout.instrument_name = "ncas-radar-camra-1";
    DSout.instrument_software = "radar-camra-rec" ;
    DSout.instrument_software_version = "1.4 Rev 58" ;

    DSout.references = "";
    DSout.source = "NCAS Mobile Ka-band Radar (Kepler)";
    DSout.comment = "";
    DSout.project = project["project_name"];
    DSout.project_principal_investigator = project["principal_investigator"]["name"];
    DSout.project_principal_investigator_email = project["principal_investigator"]["email"];
    DSout.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DSout.processing_software_url = "";
    DSout.processing_software_version = "1.0";

    DSout.scantype = "vertical_pointing";

    DSout.time_coverage_start = datetime.strftime(dt_start,'%Y-%m-%dT%H:%M:%SZ');
    DSout.time_coverage_end = datetime.strftime(dt_end,'%Y-%m-%dT%H:%M:%SZ');
    DSout.geospatial_bounds = "51.1450N -1.4384E";

    DSout.pulse_compression = "false";

    DSout.ADC_bits_per_sample = np.int(DSin.ADC_bits_per_sample);
    DSout.ADC_channels        = np.int(DSin.ADC_channels);

    # ----------------
    # Scalar variables
    # ----------------

    varin = DSin['latitude'];
    varout = DSout.createVariable('latitude',varin.datatype);
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degree_north';
    varout[:]=51.1450;

    varin = DSin['longitude'];
    varout = DSout.createVariable('longitude',varin.datatype);
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degree_east';
    varout[:]=-1.4384;

    varin = DSin['height'];
    varout = DSout.createVariable('altitude',varin.datatype);
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    varout.units = 'm';
    varout[:]=146.7;

    varout = DSout.createVariable('altitude_agl',varin.datatype);
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the elevation axis above ground';
    varout.units = 'm';
    varout[:]=16.0;

    varin = DSin['frequency'];
    varout = DSout.createVariable('frequency',varin.datatype);
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=varin[:];

    varin = DSin['prf'];
    varout = DSout.createVariable('prf',varin.datatype);
    varout.long_name = 'pulse repetition frequency';
    varout.units = 'Hz';
    varout[:]=varin[:];

    varin = DSin['beamwidthH'];
    varout = DSout.createVariable('beamwidthH',varin.datatype);
    varout.long_name = 'horizontal angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['beamwidthV'];
    varout = DSout.createVariable('beamwidthV',varin.datatype);
    varout.long_name = 'vertical angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['antenna_diameter'];
    varout = DSout.createVariable('antenna_diameter',varin.datatype);
    varout.long_name = 'antenna diameter';
    varout.units = 'm';
    varout[:]=varin[:];

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 9.0;

    varout = DSout.createVariable('antenna_focus_radial_location','f4');
    varout.long_name = 'distance along boresight from elevation axis to antenna focus';
    varout.units = 'm';
    varout[:] = 14.34;

    varin = DSin['pulse_period'];
    varout = DSout.createVariable('pulse_width',varin.datatype);
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:]=varin[:];

    varin = DSin['transmit_power'];
    varout = DSout.createVariable('transmit_power',varin.datatype);
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:]=varin[:];

    varin = DSin['clock'];
    varout = DSout.createVariable('clock',varin.datatype);
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:]=varin[:];

    varout = DSout.createVariable('clock_divfactor','i4');
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=DSin['clock'].clock_divfactor;

    varout = DSout.createVariable('delay_clocks','i4');
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:]=DSin.delay_clocks;

    varout = DSout.createVariable('samples_per_pulse','i4');
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:]=DSin.samples_per_pulse;

    varout = DSout.createVariable('pulses_per_daq_cycle','i4');
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';
    varout[:]=DSin.pulses_per_daq_cycle;

    varout = DSout.createVariable('pulses_per_ray','i4');
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:]=DSin.pulses_per_ray;

    varout = DSout.createVariable('radar_constant','f4');
    varout.long_name = 'radar constant';
    varout.units = 'dB';
    varout[:]=DSin.radar_constant;

    varout = DSout.createVariable('receiver_gain','f4');
    varout.long_name = 'receiver gain';
    varout.units = 'dB';
    varout[:]=DSin.receiver_gain;

    varout = DSout.createVariable('cable_losses','f4');
    varout.long_name = 'cable losses';
    varout.units = 'dB';
    varout[:]=DSin.cable_losses;

    varout = DSout.createVariable('extra_attenuation','f4');
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:]=DSin.extra_attenuation;


    # ---------------
    # Copy dimensions
    # ---------------
    the_dim = DSin.dimensions['time'];
    DSout.createDimension('time', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['pulses'];
    DSout.createDimension('pulse', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['samples'];
    DSout.createDimension('range', len(the_dim) if not the_dim.isunlimited() else None)

    # --------------------
    # Coordinate variables
    # --------------------
    varin = DSin['time'];
    varout = DSout.createVariable('time',varin.datatype,('time'));
    varout.standard_name = 'time';
    varout.long_name = 'time at the end of each recorded ray';
    varout.units = varin.units;
    varout[:]=varin[:];

    varin = DSin['range'];
    varout = DSout.createVariable('range',varin.datatype,('range'));
    varout.long_name = varin.long_name;
    varout.range_offset_applied = np.float32(varin.range_offset);
    varout.units = varin.units;
    varout[:]=varin[:];

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varin = DSin['elevation'];
    varout = DSout.createVariable('elevation',varin.datatype,'time');
    varout.long_name = 'elevation angle of the antenna boresight above the horizon';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);
    varout[:] = varin[:];

    varin = DSin['azimuth'];
    varout = DSout.createVariable('azimuth',varin.datatype,'time');
    varout.long_name = 'azimuth angle of the antenna boresight clockwise from the grid north';
    varout.comment = 'More generally this is the azimuth angle of the plane perpendicular to the elevation axis, which remains defined even when the elevation is 90 degree';
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);
    varout[:] = varin[:];

    # ---------------
    # Field variables
    # ---------------
    # These are:
    # SNR, VEL, RMS, LDR, NPK, SNRg, VELg, NPKg, RHO, DPS, RHOwav, DPSwav, 
    # HSDco, HSDcx, Ze, Zg, ISDRco, ISDRcx
    # The ones to use are:
    # SNR, VEL, RMS, LDR, NPK, RHO, DPS, HSDco, HSDcx, Ze, ISDRco, ISDRcx

    varin = DSin['Ze'];
    varout = DSout.createVariable('ZED',varin.datatype,('time','range'),zlib=True);
    varout.long_name = '';
    varout.units = 'dBZ';
    varout[:]=varin[:];
    comment_string  = ""
    varout.comment = comment_string;


    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: chilbolton_ts_utils.py convert_camra_ts_l0a2l0b"
    + " version:" + str(module_version));

    DSout.history = history + "\n" + DSin.history;

    DSout.last_revised_date = datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();

    return


def convert_galileo_ts_l0b2l1(infile,outfile,dBZ_offset,range_offset,data_version):

    """This routine converts raw (Level 0b) time series data from the Chilbolton 94GHz Cloud Radar (Galileo) to Level 1 data.
    Level 0b data have been generated from as-recorded Level 0a data by splitting the files along the time dimension and converting to NetCDF4.
    Processing by this routine involves removing redundant dimensions, and removing bias from the ADC samples of the received I and Q.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0b data file, e.g. `<path-to-file>/radar-galileo_<YYYYddmmHHMMSS>-<YYYYddmmHHMMSS>_fix-ts.nc4`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-w-band-1_cao_20201210-212823_fix-ts_l1_v1.0.nc`
    :type outfile: str

    :param dBZ_offset: dB calibration offset to apply to reflectivity.  This is converted to linear units and I and Q values are multiplied by the square root of this value.
    :type dBZoffset: float

    :param range_offset: Range offset to apply in metres.
    :type range_offset: float

    :param data_version: Version of data product in the format `N.m`, where `N` denotes the major verion and `m` a minor revision.
    :type data_version: str
    """

    # -----------------------------------------------------
    # Read NetCDF spectra file with embedded IQ time series
    # -----------------------------------------------------
    DSin = nc4.Dataset(infile);

    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    # -----------------------------
    # Set up NetCDF file for output
    # -----------------------------
    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')
    print(outfile)

    # ------------------------
    # Set up global attributes
    # ------------------------
    DSout.product_version = "v{}".format(data_version);
    DSout.processing_level = "1" ;

    DSout.licence = "This dataset is released for use under a Creative Commons Attribution 4.0 International (CC-BY 4.0) license (see https://creativecommons.org/licenses/by/4.0/ for terms and conditions)."
#    DSout.licence = "Data usage licence - UK Open Government Licence agreement: \n http://www.nationalarchives.gov.uk/doc/open-government-licence" ;
    DSout.acknowledgement = "This dataset was developed as part of the activity \"Doppler Wind Radar Science Performance Study (WIVERN-2)\", funded by the European Space Agency under Contract no. 4000130864/20/NL/CT.  Users should acknowledge UK Research and Innovation as the data provider (in partnership with the National Centre for Atmospheric Science)." ;
#    DSout.acknowledgement = "Acknowledgement is required of UK Research and Innovation as the data provider (in partnership with the National Centre for Atmospheric Science) whenever and wherever these data are used." ;
    DSout.platform = "Chilbolton Atmospheric Observatory" ;
    DSout.platform_type = "stationary_platform" ;
    DSout.title = "Time series from 94 GHz Galileo radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";
    DSout.creator_name = "Chris Walden" ;
    DSout.creator_email = "chris.walden@ncas.ac.uk" ;
    DSout.creator_url = "https://orcid.org/0000-0002-5718-466X" ;
    DSout.institution = "National Centre for Atmospheric Science (NCAS)";
    DSout.instrument_name = "ncas-radar-w-band-1";
    DSout.instrument_software = "radar-galileo-rec" ;
    DSout.instrument_software_version = "" ;

    DSout.references = "";
    DSout.source = "94GHz Galileo Radar";
    DSout.comment = "";
    DSout.project = "WIVERN-2 Doppler Wind Radar Science Performance Study";
    DSout.project_principal_investigator = "Anthony Illingworth";
    DSout.project_principal_investigator_email = "a.j.illingworth@reading.ac.uk";
    DSout.project_principal_investigator_url = "https://orcid.org/0000-0002-5774-8410";

    DSout.processing_software_url = "https://github.com/longlostjames/wivern_chilbolton_utils.git";
    DSout.processing_software_version = "1.0";

    DSout.scantype = "vertical_pointing";

    DSout.time_coverage_start = datetime.strftime(dt_start,'%Y-%m-%dT%H:%M:%SZ');
    DSout.time_coverage_end = datetime.strftime(dt_end,'%Y-%m-%dT%H:%M:%SZ');
    DSout.geospatial_bounds = "51.1447N -1.4384E";

    DSout.pulse_compression = "false";

    DSout.ADC_bits_per_sample = np.int(12);
    DSout.ADC_channels        = np.int(8);

    # ----------------
    # Scalar variables
    # ----------------

    varout = DSout.createVariable('latitude','f4');
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degree_north';
    varout[:]=51.1447;

    varout = DSout.createVariable('longitude','f4');
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degree_east';
    varout[:]=-1.4385;

    varout = DSout.createVariable('altitude','f4');
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the antenna above the geoid (WGS84)';
    varout.units = 'm';
    varout[:] = 131.6;

    varout = DSout.createVariable('altitude_agl','f4');
    varout.long_name = 'altitude of the antenna above ground';
    varout.units = 'm';
    varout[:] = 0.9;

    varout = DSout.createVariable('frequency','f4');
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=94.008;

    varin = DSin['prf'];
    varout = DSout.createVariable('prf',varin.datatype);
    varout.long_name = 'pulse repetition frequency';
    varout.units = 'Hz';
    varout[:]=varin[:];

    varin = DSin['beamwidthH'];
    varout = DSout.createVariable('beamwidthH',varin.datatype);
    varout.long_name = 'horizontal angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['beamwidthV'];
    varout = DSout.createVariable('beamwidthV',varin.datatype);
    varout.long_name = 'vertical angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['antenna_diameter'];
    varout = DSout.createVariable('antenna_diameter',varin.datatype);
    varout.long_name = 'antenna diameter';
    varout.units = 'm';
    varout.comment = "Refers to diameter of each antenna in bistatic pair.  Separation of antennae is 0.66m."
    varout[:] = varin[:];

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 0.18;

    varin = DSin['pulse_period'];
    varout = DSout.createVariable('pulse_width',varin.datatype);
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:] = varin[:];

    varin = DSin['transmit_power'];
    varout = DSout.createVariable('transmit_power',varin.datatype);
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:] = varin[:];

    varin = DSin['clock'];
    varout = DSout.createVariable('clock',varin.datatype);
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:] = varin[:];

    clock_divfactor = DSin['clock'].clock_divfactor;
    varout = DSout.createVariable('clock_divfactor','i4');
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=clock_divfactor;

    delay_clocks = DSin.getncattr('delay_clocks');
    varout = DSout.createVariable('delay_clocks','i4');
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:] = delay_clocks;

    samples_per_pulse = DSin.getncattr('samples_per_pulse');
    varout = DSout.createVariable('samples_per_pulse','i4');
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:] = samples_per_pulse;

    pulses_per_daq_cycle = DSin.getncattr('pulses_per_daq_cycle');
    varout = DSout.createVariable('pulses_per_daq_cycle','i4');
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';
    varout[:] = pulses_per_daq_cycle;

    pulses_per_ray = DSin.getncattr('pulses_per_ray');
    varout = DSout.createVariable('pulses_per_ray','i4');
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:] = pulses_per_ray;

    varout = DSout.createVariable('radar_constant','f4');
    varout.long_name = 'radar constant';
    varout.units = 'dB';

    varout = DSout.createVariable('receiver_gain','f4');
    varout.long_name = 'receiver gain';
    varout.units = 'dB';

    varout = DSout.createVariable('cable_losses','f4');
    varout.long_name = 'cable losses';
    varout.units = 'dB';

    varout = DSout.createVariable('extra_attenuation','f4');
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:] = 0.0;


    # ---------------
    # Copy dimensions
    # ---------------
    the_dim = DSin.dimensions['time'];
    DSout.createDimension('time', len(the_dim) if not the_dim.isunlimited() else None)

    fft_bin_dim = DSin.dimensions['fft_bin_dim'];
    spectra_number_dim = DSin.dimensions['spectra_number_dim'];
    DSout.createDimension('pulse', len(fft_bin_dim)*len(spectra_number_dim));

    the_dim = DSin.dimensions['range'];
    DSout.createDimension('range', len(the_dim) if not the_dim.isunlimited() else None)

    print(DSout.dimensions['pulse']);

    # --------------------
    # Coordinate variables
    # --------------------
    varin = DSin['time'];
    varout = DSout.createVariable('time',varin.datatype,('time'));
    varout.standard_name = 'time';
    varout.long_name = 'time at the start of each recorded ray';
    varout.units = varin.units;
    varout[:]=varin[:];

    varin = DSin['range'];
    varout = DSout.createVariable('range',varin.datatype,('range'));
    varout.long_name = 'distance from the antenna to the middle of each range gate';
    varout.range_offset_applied = np.float32(range_offset);
    varout.units = varin.units;
    varout[:]=varin[:]+range_offset-varin.range_offset;

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varout = DSout.createVariable('elevation','f4','time');
    varout.long_name = "elevation angle above the horizon of the antenna boresight";
    varout.comment   = "assumes transmit and receive antenna boresights are aligned";
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);

    varout = DSout.createVariable('azimuth','f4','time');
    varout.long_name   = "azimuth angle clockwise from grid north of the plane containing the antenna boresight and zenith vectors";
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);
    varout.comment = "assumes transmit and receive antenna boresights are aligned"

    # --------------------------------
    # Determine bias-corrected I and Q
    # --------------------------------
    # Input Level 0b file has I and Q dependent on the following dimensions (time,spectra,pulses,sample).
    # We rearrange this into a three-dimensional array dependent on (time,pulse,range);

    Ico_in = DSin['IPF_HH'][:,:,:,:];
    Qco_in = DSin['QPF_HH'][:,:,:,:];

    Icx_in = DSin['IPF_HV'][:,:,:,:];
    Qcx_in = DSin['QPF_HV'][:,:,:,:];

    nray     = Ico_in.shape[0];
    nspectra = Ico_in.shape[1];
    npulse   = Ico_in.shape[2];
    ngate    = Ico_in.shape[3];

    Ico = Ico_in.reshape(nray,npulse*nspectra,ngate);
    Qco = Qco_in.reshape(nray,npulse*nspectra,ngate);

    Ico_mean = np.mean(Ico,axis=1);
    Qco_mean = np.mean(Qco,axis=1);

    Ico[:,:,:] = Ico[:,:,:]-Ico_mean[:,None,:];
    Qco[:,:,:] = Qco[:,:,:]-Qco_mean[:,None,:];

    Zcal=dBZ_offset+10*np.log10(DSout['prf'][:]/512);

    Ico=Ico*np.sqrt(10**(Zcal/10.0));
    Qco=Qco*np.sqrt(10**(Zcal/10.0));

    Icx = Icx_in.reshape(nray,npulse*nspectra,ngate);
    Qcx = Qcx_in.reshape(nray,npulse*nspectra,ngate);

    Icx_mean = np.mean(Icx,axis=1);
    Qcx_mean = np.mean(Qcx,axis=1);

    Icx[:,:,:] = Icx[:,:,:]-Icx_mean[:,None,:];
    Qcx[:,:,:] = Qcx[:,:,:]-Qcx_mean[:,None,:];

    Icx=Icx*np.sqrt(10**(Zcal/10.0));
    Qcx=Qcx*np.sqrt(10**(Zcal/10.0));

    # ---------------
    # Field variables
    # ---------------

    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Ico[:]);
#    scale_factor = (np.max(Ico[:])-np.min(Ico[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((Ico[:,:,:] - add_offset)/scale_factor));
    #varout.scale_factor = np.float32(scale_factor);
    #varout.add_offset = np.float32(add_offset);
    #varout[:] = packed_data;
    varout[:] = Ico;

    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Qco);
#    scale_factor = (np.max(Qco)-np.min(Qco)) / (2**16-1)
#    packed_data = np.int16(np.rint((Qco - add_offset)/scale_factor));
    #varout.scale_factor = np.float32(scale_factor);
    #varout.add_offset = np.float32(add_offset);
    #varout[:] = packed_data;
    varout[:] = Qco;

    varout = DSout.createVariable('Icx','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Icx);
#    scale_factor = (np.max(Icx)-np.min(Icx)) / (2**16-1);
#    packed_data = np.int16(np.rint((Icx - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
#    varout[:] = packed_data;
    varout[:] = Icx;

    varout = DSout.createVariable('Qcx','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Qcx);
#    scale_factor = (np.max(Qcx)-np.min(Qcx)) / (2**16-1);
#    packed_data = np.int16(np.rint((Qcx - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
#    varout[:] = packed_data;
    varout[:] = Qcx;

    varout = DSout.createVariable('qc_flag','u1',('time','pulse','range'));
    varout.is_quality = 'true';
    varout.qualified_variables = 'I Q Icx Qcx';
    varout.long_name = 'Quality control flag';
    varout.flag_values = np.uint8(0),np.uint8(1), np.uint8(2), np.uint8(3), np.uint8(4), np.uint8(255);
    varout.flag_meanings = 'not_used good_data probably_good_data bad_data data_in_blind_range no_qc_performed';
    varout[:] = 2;

    blind_range = np.arange(14);
    DSout['qc_flag'][:,:,blind_range] = 4;

    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: kepler_utils.convert_mmclx2l0b"
    + " version:" + str(module_version));

    DSout.history = history + "\n" + DSin.history;

    DSout.last_revised_date = datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();

    return

def process_kepler(datestr,inpath,outpath,yaml_project_file,tracking_tag):

    pattern = '*{}*_fix-ts.nc'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    tsfiles = [];
    tsdirs = [];

    for root,dirs,files in os.walk(datepath):
        tsfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        tsdirs += dirs;

    data_version = "1.0";

    dBZ_offset = 9.0;
    range_offset = -865.56+864.0;

    l0bpath = os.path.join(outpath,'L0b',datestr);
    l1path = os.path.join(outpath,'L1',datestr);

    os.makedirs(l0bpath,exist_ok=True);
    os.makedirs(l1path,exist_ok=True);

    print(tsdirs);
    for dir in tsdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l0bpath,dir),exist_ok=True);
        os.makedirs(os.path.join(l1path,dir),exist_ok=True);

    for f in tsfiles:
        outfile_splits = os.path.split(f);

        outfile_string = outfile_splits[1].replace('.nc','_l0b.nc');
        splits = outfile_string.split('_');
        instrument_name =splits[0].replace('radar-camra','ncas-radar-camra-1');
        platform = 'cao';
        datestr = splits[1][0:8];
        timestr = splits[1][8:];
        level = splits[3].split('.')[0];
        l0bfile = '{}_{}_{}-{}_{}_{}_v{}.nc'.format(instrument_name,platform,datestr,timestr,splits[2],level,data_version)

        l0bfile.replace('radar-camra','ncas-radar-camra-1_cao');

        event = outfile_splits[0].split('/')[-1];

        if not event in tsdirs:
            l0bfile = os.path.join(outpath,'L0b',datestr,l0bfile);
        else:
            l0bfile = os.path.join(outpath,'L0b',datestr,event,l0bfile);

        convert_camra_ts_l0a2l0b(f,l0bfile,yaml_project_file,tracking_tag);

        l1file = l0bfile.replace('l0b','l1');
        l1file = l1file.replace('L0b','L1');

        print(l0bfile);
        print(l1file);

        convert_camra_ts_l0b2l1(l0bfile,l1file,dBZ_offset,ZDR_offset,range_offset,data_version,yaml_project_file,tracking_tag);

    return


# ------------------------------------------------------------------------
# Define function to produce cfradial files from mmclx files
# ------------------------------------------------------------------------
def mmclx_to_cfradial(mclxfile):

    user = getpass.getuser()

    print('Opening NetCDF file ' + mmclxfile)
    DS = nc4.Dataset(ncfile,'r+')

    scantype = DS.getncattr('scantype');

    dataset.close();

    # --------------------------------
    # Open NetCDF file using arm_pyart
    # --------------------------------
    if (scantype=='RHI'):
        radar = pyart.aux_io.read_camra_rhi(ncfile);
    elif (scantype=='PPI'):
        radar = pyart.aux_io.read_camra_ppi(ncfile);

    # -----------------------------------
    # Write cfradial file using arm_pyart
    # -----------------------------------
    cfradfile=ncfile.replace(".nc","-cfrad.nc");

    pyart.io.cfradial.write_cfradial(cfradfile, radar, format='NETCDF4',
        time_reference=False)
