``options.calc_fluo``
=======================

fluorescence_scalars.csv
---------------------------

rows - time (simulation number)

columns - variables

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **F_1stpeak**
      - W m-2 um-1 sr-1
      - rad.F685
    * - **wl_1stpeak**
      - nm
      - rad.wl685
    * - **F_2ndpeak**
      - W m-2 um-1 sr-1
      - rad.F740
    * - **wl_2ndpeak**
      - nm
      - rad.wl740
    * - **F687**
      - W m-2 um-1 sr-1
      - rad.F687
    * - **F760**
      - W m-2 um-1 sr-1
      - rad.F760
    * - **LFtot**
      - W m-2 um-1 sr-1
      - rad.LoutF
    * - **EFtot**
      - W m-2
      - rad.EoutF
    * - **EFtot_RC**
      - W m-2
      - rad.EoutFrc

fluorescence.csv
-------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl 640:1:850 nm

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF_**
      - W m-2 um-1 sr-1
      - fluorescence per wavelength in observation direction


sigmaF.csv
-------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl 640:1:850 nm

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **sigmaF**
      - \-
      - escape probability

fluorescence_hemis.csv
------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl 640:1:850 nm

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **EoutF_**
      - W m-2 um-1
      - TOC hemispherically integrated fluorescence


fluorescence_ReabsCorr.csv
-----------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl 640:1:850 nm

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **EoutFrc_**
      - W m-2 um-1
      - fluorescence_spectrum 640:1:850 nm reabsorption corrected
	  
	  
fluorescence_AllLeaves.csv
-----------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl 640:1:850 nm

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Femleaves_**
      - W m-2 um-1
      - fluorescence_spectrum 640:1:850 nm emission by all leaves


Lo_spectrum_inclF.csv
-----------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Lototf_**
      - W m-2 um-1 sr-1
      - upwelling radiance in observation direction including fluorescence



apparent_reflectance.csv
---------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - wl

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **(\Lo_ + \LoF_) * pi  / (\Esun_ + \Esky_)**
      - \-
      - fraction of radiation in observation direction + emitted fluorescence \* pi / irradiance