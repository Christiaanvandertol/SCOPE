In each simulation
===================

.. contents::


pars_and_input_short.csv
-----------------------------

rows - timestep

columns - values of parameters which values where changing during the simulation.


vegetation.csv
----------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **simulation_number**
      - \-
      - time step counter
    * - **year**
      - \-
      - year
    * - **DoY**
      - \-
      - decimal day of year (DOY)
    * - **aPAR**
      - umol m-2 s-1
      - absorbed PAR by leaves
    * - **aPARbyCab**
      - umol m-2 s-1
      - absorbed PAR by chlorophylls a, b
    * - **aPARbyCab(energyunits)**
      - W m-2
      - absorbed PAR
    * - **Photosynthesis**
      - umol m-2 s-1
      - net photosynthesis of canopy (Actot)
    * - **Electron_transport**
      - umol m-2 s-1
      - electron transport rate (Ja)
    * - **NPQ_energy**
      - W m-2
      - non-photochemical quenching (energy)
    * - **LST**
      - K
      - land surface temperature

fluxes.csv
------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **simulation_number**
      - \-
      - time step counter
    * - **nu_iterations**
      - \-
      - number of iterations in energy balance
    * - **year**
      - \-
      - year
    * - **DoY**
      - \-
      - decimal day of year (DOY)
    * - **Rnctot**
      - W m-2
      - net radiation of canopy
    * - **lEctot**
      - W m-2
      - latent heat flux of canopy (transpiration)
    * - **Hctot**
      - W m-2
      - sensible heat of canopy
    * - **Actot**
      - umol m-2 s-1
      - net photosynthesis of canopy
    * - **Tcave**
      - ºC
      - 'average' canopy temperature
    * - **Rnstot**
      - W m-2
      - net radiation of soil
    * - **lEstot**
      - W m-2
      - latent heat flux of soil (evaporation)
    * - **Hstot**
      - W m-2
      - sensible heat of soil
    * - **Gtot**
      - W m-2
      - soil heat flux
    * - **Tsave**
      - ºC
      - 'average' soil temperature
    * - **Rntot**
      - W m-2
      - total net radiation
    * - **lEtot**
      - W m-2
      - total latent heat flux
    * - **Htot**
      - W m-2
      - total sensible heat
    * - **rss**
      - s m-1
      - soil resistance to evaporation

radiation.csv
---------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **simulation_number**
      - \-
      - time step counter
    * - **year**
      - \-
      - year
    * - **DoY**
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

Reflectance and its quantiles
'''''''''''''''''''''''''''''''
For the meaning of reflectance factors, please, refer to :ref:`my_proposal/brdf:Definition`

reflectance.csv
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

rsd.csv
-----------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **rsd**
      - \-
      - directional-hemispherical reflectance factor

rdd.csv
-----------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **rdd**
      - \-
      - bi-hemispherical reflectance factor

rso.csv
-----------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **rso**
      - \-
      - bi-directional reflectance factor

rdo.csv
-----------------

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **rdo**
      - \-
      - hemispherical-directional reflectance factor

Radiation per wavelength
'''''''''''''''''''''''''''''''

Eout_spectrum.csv
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

Lo_spectrum.csv
-----------------------------

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Lotot_**
      - W m-2 um-1 sr-1
      - radiance spectrum in observation direction excluding fluorescence


Esun.csv
------------------------------

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Esun_**
      - W m-2 um-1 sr-1
      - direct top of canopy irradiance irradiance

Esky.csv
------------------------------

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Esky_**
      - W m-2 um-1 sr-1
      - diffuse top of canopy irradiance irradiance
