data
======

.. contents::

input
------

dataset for_verification
""""""""""""""""""""""""""
‘Dataset for_verification’ contains time series of meteorological data. In this case, half-hourly data are provided. It is possible to work with any time interval, but due to the thermal inertia of the soil, the calculation of soil temperature may not be accurate when the time interval is longer than three hours.
It is recommended to name your own dataset ‘dataset sitename or projectname’.
The directory contains the following compulsory files (all in ASCII format):

-	A time vector (``t_.dat``):  a vertical array of time values, in decimal days of the year [1:366.99]. For example, 10 January 2009, 12:00 would be: 10.5.  All other files (see below) should correspond to this time vector (and thus have the same size).
-	A year vector (``year_.dat``):  the year corresponding to the time vector.  For example, 10 January 2009, 12:00 would be: 2009
-	TOC incoming shortwave radiation (``Rin_.dat``):  a broadband (0.3 to 2.5 μm) measurement of incoming shortwave radiation (W m-2), perpendicular to the surface.
-	TOC incoming long wave radiation (``Rli_.dat``):  a broadband (2.5 to 50 μm) measurement of incoming long wave radiation (W m-2), perpendicular to the surface.
-	Air pressure (``p_.dat``):  air pressure (hPa or mbar)
-	Air temperature measured above the canopy (``Ta_.dat``):  air temperature above the canopy in °C.
-	Vapour pressure measured above the canopy (``e_.dat``):  vapour pressure above the canopy (hPa or mbar).
-	Wind speed (``u_.dat``):  wind speed measured above the canopy (m s-1)

The following additional files (not compulsory) can be added:

-	Carbon dioxide concentration measured above the canopy (mg  m-3)

And the following tables (not compulsory):

-	Leaf area index (``LAI_.dat``)
-	Measurement height (``z_.dat``) (m)
-	Vegetation height (``h_.dat``) (m)
-	Maximum carboxylation capacity (``Table_Vcmax_.dat``)
-	Chlorophyll content file (``Table_Cab_.dat``)

If a table is not present, then the corresponding a priori value specified in the file input_data.xlsx file is used instead. It is only useful to create the tables LAI_dat etc. if the leaf area index, measurement height, vegetation height etc. change with time during the measurement period.

A table has a slightly different format than the other input files. A table has two columns: the first column contains the decimal DOY, the second column contains the value of the variable. The reason why tables have a different format is that the variables in the table are usually not measured at the same time interval as the meteorological input. For example, the LAI may be measured only once per month.

An example of a table can be found in ‘dataset for_verification’.
The measurement height is only relevant for wind speed, vapour pressure and the carbon dioxide concentration. It is currently not possible to specify separate measurement heights for each of these variables.

The carbon dioxide concentration must be provided in mg m-3. This is a commonly used unit in most data sets. SCOPE automatically converts this to ppm and to umol m-3 internally. If the carbon dioxide concentration file is not provided, SCOPE assumes a constant concentration corresponding to 380 ppm.

.. note::
 It is important that all files except for the tables have equal length, and that all measurements correspond to the time vector. A Julian calendar is used. The time zone should be provided (the difference between the local time in the file and UTC or GMT. Input files should be comma separated, space separated or tab separated ASCII files. They should not contain empty lines or comment lines.

In case SCOPE is run in individual mode, then the meteorological input files are not used. In that case, all meteorological data are taken from the spreadsheet.



directional
""""""""""""

The input in the directory ‘directional’ is only used for multi-angle simulations (if the option ‘directional’ is switched on in parameters. In this directory one can provide the observer’s zenith and azimuth angles. The files in this directory have two columns: the first column is the observer’s zenith angle, the second the observer’s azimuth (relative to that of the sun, counterclockwise), both in degrees. If the option ``directional`` is switched on, SCOPE will calculate the radiance spectrum in all directions given in the input file.

fluspect_parameters
""""""""""""""""""""

In this directory, absorption spectra of different leaf components are provided, according to PROSPECT 3.1, as well as Fluspect input: standard spectra for PSI and PSII.

leafangles
"""""""""""

In this folder, example leaf inclination distribution data are provided. It is possible to use these distributions instead of the leafinclination model of Verhoef (1998 :cite:`Verhoef1998`) with the two parameters LIDFa and LIDFb. In that case, provide a filename in the ``input_data.xlsx`` tab **filenames** or the file ``filenames.m``.

radiationdata
""""""""""""""

``RTMo.m`` calculates spectra based on MODTRAN5 outputs (T-18 system, :cite:`Verhoef2018`). One .atm (atmospheric) file is provided in the data, 12 more are provided separately in a different .zip folder (in order to minimize the size of the SCOPE package, these are not provided standard). Note that in the input data  (files as well as the spreadsheet), the broadband input radiation may be provided. SCOPE linearly scales the input spectra of the optical and the thermal domain in such a way, that the spectrally integrated input shortwave and long wave radiation matches with the measured values. A limitation of this approach is that the same shape of the input spectrum is used independent on the atmospheric conditions. If this scaling is not wanted, then leave ‘Rin’ and ‘Rli’ empty in the spreadsheet.

.. Note::
    In earlier versions of the model (1.34 and older), two input spectra of solar and sky radiation were provided (rad.txt and rad2.txt) in this directory. The data were calculated with MODTRAN4. The ASCII file in this directory consisted of three columns containing the following. The first column contained the wavelength in nm, the second column the solar radiation in W m-2 μm-1, and the third column the sky radiation in W m-2 μm-1. These data are now obsolete (since version 1.40).


soil_spectrum
""""""""""""""
In this directory, the soil spectrum is provided. The ASCII file in this directory consists of two columns containing the following: The first column contains the wavelength in μm, the following columns reflectance spectra. Note that it is also possible to simulate a soil reflectance spectrum with the BSM model. In that case the values for the BSM model parameters are taken from the input data, and the archived spectra in this folder are not used.

measured
---------

The validation data are stored in directory ‘measured’. It is up to the user to organize this directory.