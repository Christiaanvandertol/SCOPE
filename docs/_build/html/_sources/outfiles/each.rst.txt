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
    * - **Photosynthesis**
      - umol m-2 s-1
      - net photosynthesis of canopy (Actot)
    * - **Electron_transport**
      - umol m-2 s-1
      - electron transport rate (Ja)
    * - **NPQ_energy**
      - W m-2
      - non-photochemical quenching [energy units]
    * - **NPQ_photon**
      - umol m-2 s-1
      - non-photochemical quenching
    * - **canopy_level_FQE**
      - umol photons (umol photons)-1
      - fluorescence quantum yield at canopy level
    * - **LST**
      - K
      - land surface temperature
    * - **emis**
      - \-
      - thermal energy emissivity

aPAR.csv
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
    * - **iPAR**
      - umol m-2 s-1
      - intercepted PAR
    * - **iPARE**
      - W m-2
      - intercepted PAR [energy units]
    * - **LAIsunlit**
      - m2 m-2
      - leaf area index of sunlit leaves
    * - **LAIshaded**
      - m2 m-2
      - leaf area index of shaded leaves
    * - **aPARtot**
      - umol m-2 s-1
      - absorbed PAR by all leaves
    * - **aPARsun**
      - umol m-2 s-1
      - absorbed PAR by sunlit leaves
    * - **aPARsha**
      - umol m-2 s-1
      - absorbed PAR by shaded leaves
    * - **aPARCabtot**
      - umol m-2 s-1
      - absorbed PAR by chlorophylls a,b of all leaves
    * - **aPARCabsun**
      - umol m-2 s-1
      - absorbed PAR by chlorophylls a,b of sunlit leaves
    * - **aPARCabsha**
      - umol m-2 s-1
      - absorbed PAR by chlorophylls a,b of shaded leaves
    * - **aPARCartot**
      - umol m-2 s-1
      - absorbed PAR by carotenoinds of all leaves
    * - **aPARCarsun**
      - umol m-2 s-1
      - absorbed PAR by carotenoinds of sunlit leaves
    * - **aPARCarsha**
      - umol m-2 s-1
      - absorbed PAR by carotenoinds of shaded leaves
    * - **aPARtotE**
      - W m-2
      - absorbed PAR by all leaves [energy units]
    * - **aPARsunE**
      - W m-2
      - absorbed PAR by sunlit leaves [energy units]
    * - **aPARshaE**
      - W m-2
      - absorbed PAR by shaded leaves [energy units]
    * - **aPARCabtotE**
      - W m-2
      - absorbed PAR by chlorophylls a,b of all leaves [energy units]
    * - **aPARCabsunE**
      - W m-2
      - absorbed PAR by chlorophylls a,b of sunlit leaves [energy units]
    * - **aPARCabshaE**
      - W m-2
      - absorbed PAR by chlorophylls a,b of shaded leaves [energy units]
    * - **aPARCartotE**
      - W m-2
      - absorbed PAR by carotenoinds of all leaves [energy units]
    * - **aPARCarsunE**
      - W m-2
      - absorbed PAR by carotenoinds of sunlit leaves [energy units]
    * - **aPARCarshaE**
      - W m-2
      - absorbed PAR by carotenoinds of shaded leaves [energy units]


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


resistances.csv
-----------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **raa**
      - s m-1
      - aerodynamic resistance above the canopy
    * - **raws**
      - s m-1
      - aerodynamic resistance within the soi
    * - **rss**
      - s m-1
      - soil resistance to evaporation
    * - **ustar**
      - m s-1
      - friction velocity u*


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
      - W m-2 um-1
      - direct top of canopy irradiance

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
      - W m-2 um-1
      - diffuse top of canopy irradiance
