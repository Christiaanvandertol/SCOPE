In each simulation
===================

.. contents::

aerodyn.dat
---------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **raa**
      - s m-1
      - total aerodynamic resistance above canopy
    * - **rawc**
      - s m-1
      - canopy total aerodynamic resistance below canopy
    * - **raws**
      - s m-1
      - soil total aerodynamic resistance below canopy
    * - **ustar**
      - m s-1
      - friction velocity


BOC_irradiance.dat
---------------------

BOC - bottom of canopy (61st layer)

rows - timestep

First 2162 columns: shaded fraction.

Last 2162 columns: average BOC irradiance.

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Emin_(61, :)**
      - W m-2 um-1
      - irradiance at the bottom of the canopy for the shaded fraction
    * - **\Emin_(61, :) + \Esun_(61, :) * gap.Ps(61, :)**
      - W m-2 um-1
      - average BOC irradiance (sunlit + shaded fraction)


fluxes.dat
------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **timestep**
      - \-
      - time step counter
    * - **counter**
      - \-
      - number of iterations in energy balance
    * - **year**
      - \-
      - year
    * - **T**
      - \-
      - decimal day of year (DOY)
    * - **Rntot**
      - W m-2
      - total net radiation
    * - **lEtot**
      - W m-2
      - total latent heat flux
    * - **Htot**
      - W m-2
      - total sensible heat
    * - **Rnctot**
      - W m-2
      - net radiation of canopy
    * - **lEctot**
      - W m-2
      - latent heat flux of canopy
    * - **Hctot**
      - W m-2
      - sensible heat of canopy
    * - **Actot**
      - umol m-2 s-1
      - net photosynthesis of canopy
    * - **Rnstot**
      - W m-2
      - net radiation of soil
    * - **lEstot**
      - W m-2
      - latent heat flux of soil
    * - **Hstot**
      - W m-2
      - sensible heat of soil
    * - **Gtot**
      - W m-2
      - soil heat flux
    * - **Resp**
      - umol m-2 s-1
      - soil respiration rate
    * - **aPAR_Cab**
      - umol m-2 s-1
      - absorbed PAR by chlorophylls a, b
    * - **aPAR**
      - umol m-2 s-1
      - absorbed PAR by leaves
    * - **fPAR**
      - \-
      - fraction of absorbed PAR by canopy, excluding soil
    * - **aPAR_energyunits**
      - W m-2
      - absorbed PAR
    * - **iPAR**
      - W m-2
      - incident PAR
    * - **fluortot**
      - W m-2
      - hemispherically and spectrally integrated chlorophyll fluorescence at the top
    * - **fluor_yield**
      - W W-1
      - Fluortot / aPAR_energyunits


irradiance_spectra.dat
------------------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Rin * (fEsuno + fEskyo)**
      - W m-2 um-1
      - spectrum of incoming radiation used in the simulation


pars_and_input.dat
----------------------

rows - timestep

columns - all input parameters from ``input_data.xlxs``


pars_and_input_short.dat
-----------------------------

rows - timestep

columns - ``Cab, Vcmo, LAI, hc,	zo,	d,	z,	Rin, Ta,	Rli,	p,	ea,	u,	Ca,	tts,	SMC``


radiation.dat
---------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **timestep**
      - \-
      - time step counter
    * - **year**
      - \-
      - year
    * - **T**
      - \-
      - decimal day of year (DOY)
    * - **ShortIn (Rin)**
      - W m-2
      - Incoming shortwave radiation (copy from input)
    * - **LongIn (Rli)**
      - W m-2
      - Incoming longwave radiation (copy from input)
    * - **HemisOutShort (Eouto)**
      - W m-2
      - hemispherical outgoing shortwave radiation
    * - **HemisOutLong (Eoutt + Eoutte)**
      - W m-2
      - hemispherical outgoing longwave radiation
    * - **HemisOutTot (Eouto + Eoutt + Eoutte)**
      - W m-2
      - total hemispherical outgoing radiation
    * - **Net (Rntot)**
      - W m-2
      - total net radiation


reflectance.dat
-----------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **\Lo_ * pi  / (\Esun_ + \Esky_)**
      - \-
      - fraction of radiation in observation direction \* pi / irradiance


spectrum_hemis_optical.dat
------------------------------

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Eout_**
      - W m-2 um-1
      - hemispherical outgoing radiation spectrum

spectrum_obsdir_optical.dat
-----------------------------

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Lo_**
      - W m-2 um-1 sr-1
      - radiance spectrum in observation direction

surftemp.dat
----------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **timestep**
      - \-
      - time step counter
    * - **year**
      - \-
      - year
    * - **T**
      - \-
      - decimal day of year (DOY)
    * - **Ta**
      - ºC
      - Air temperature above the canopy
    * - **Tss(1)**
      - ºC
      - Surface temperature of shaded soil
    * - **Tss(2)**
      - ºC
      - Surface temperature of sunlit soil
    * - **Tcave**
      - ºC
      - canopy weighted average temperature
    * - **Tsave**
      - ºC
      - soil weighted average temperature

wl.dat
---------

single row (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **wl**
      - nm
      - wavelengths of the spectral output files
