Software description
====================

.. _raw-data-conversion:

Raw data conversion
-------------------

The 3 GHz, 35 GHz and 94 GHz radars involved in the WIVERN-2 campaign provide
raw data (designated as Level 0a) in different formats.  For the 3 GHz and
94 GHz radars the files are in NetCDF-3 format, but with differing content.
The 35 GHz radar Level 0a data are in binary format.

.. automodule:: chilbolton_ts_utils

For the 3 GHz radar the Level 0a data are first processed to produce Level 0b
files.  This is accomplished using the following routine:

.. autofunction:: chilbolton_ts_utils.convert_camra_ts_l0a2l0b

94 GHz radar time series are recorded as hourly files. Due to the large size of
these files the following function uses the Python `pynco`_ module to split them
into smaller, more manageable Level 0b files.

.. autofunction:: chilbolton_ts_utils.split_ncfile

Where required the following function is used to remove empty rays of data at
the beginning or end of a file.

.. autofunction:: chilbolton_ts_utils.trim_ncfile

Conversion from Level 0b to Level 1 is achieved using the following functions:

.. autofunction:: chilbolton_ts_utils.convert_camra_ts_l0b2l1

.. autofunction:: chilbolton_ts_utils.convert_galileo_ts_l0b2l1

Binary format 35 GHz radar files are converted directly to Level 1 using

.. autofunction:: chilbolton_ts_utils.convert_copernicus_ts_l0a2l1


.. _pynco: https://pynco.readthedocs.io/


.. _quicklook-generation:

Quicklook generation
-------------------

