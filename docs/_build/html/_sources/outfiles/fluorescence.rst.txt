``options.calc_fluo``
=======================

fluorescence.dat
-------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - fluorescence from both photosystems in observation direction

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF_**
      - W m-2 um-1 sr-1
      - fluorescence per wavelength in observation direction

fluorescence_emitted_by_all_leaves.dat
-----------------------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - total emitted fluorescence by all leaves.
Within canopy scattering / re-absorption is omitted.
Within leaf scattering / re-absorption is taken into account.

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Fem_**
      - W m-2 um-1
      - hemispherical emitted fluorescence by all leaves


fluorescence_emitted_by_all_photosystems.dat
-----------------------------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - total emitted fluorescence by all photosystems for wavelengths
Within canopy scattering / re-absorption is omitted.
Within leaf scattering / re-absorption is omitted.


.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Femtot**
      - W m-2 um-1
      - hemispherical emitted fluorescence by all photosystems per wavelengths (excluding leaf and canopy re-absorption and scattering)

fluorescence_hemis.dat
------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - top of canopy (TOC) hemispherical fluorescence

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Fhem_**
      - W m-2 um-1
      - TOC hemispherical fluorescence

fluorescence_scattered.dat
-----------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - top of canopy (TOC) fluorescence contribution from leaves and soil after scattering

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **sum(LoF_scattered) + sum(LoF_soil)**
      - W m-2 um-1 sr-1
      - TOC directional fluorescence from leaves and soil after scattering


fluorescence_shaded.dat
--------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - top of canopy (TOC) fluorescence contribution from shaded leaves in observer direction per wavelengths

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF_shaded**
      - W m-2 um-1 sr-1
      - TOC fluorescence from shaded leaves in observer direction


fluorescence_sunlit.dat
-------------------------

.. Note:: ``options.calc_fluor``

rows - time (simulation number)

columns - top of canopy (TOC) fluorescence contribution from sunlit leaves in observer direction per wavelengths

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF_sunlit**
      - W m-2 um-1 sr-1
      - TOC fluorescence from sunlit leaves in observer direction


fluorescencePSI.dat
----------------------

.. Note:: ``options.calc_fluor`` && ``options.calc_PSI``

rows - time (simulation number)

columns - fluorescence of photosystem I (PSI) per wavelength in observation direction

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF1_**
      - W m-2 um-1 sr-1
      - fluorescence of PSI per wavelength in observation direction


fluorescencePSII.dat
----------------------

.. Note:: ``options.calc_fluor`` && ``options.calc_PSI``

rows - time (simulation number)

columns - fluorescence of photosystem II (PSII) per wavelength in observation direction

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LoF2_**
      - W m-2 um-1 sr-1
      - fluorescence of PSII per wavelength in observation direction
