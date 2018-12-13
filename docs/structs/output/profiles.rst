profiles
=========

Vertical profiles of temperatures and fluxes

.. Note::
    This is an optional struct that requires ``options.calc_vert_profiles`` to be calculated and output

Calculated
""""""""""""

:func:`.RTMo`

Variations
""""""""""""

:func:`.RTMf` if ``options.calc_fluor``

:func:`.ebal` if ``options.calc_ebal``

Output file
""""""""""""

if ``options.calc_vert_profiles``

- :ref:`outfiles/profiles:layer_aPAR.dat`
- :ref:`outfiles/profiles:layer_aPAR_Cab.dat`

if ``options.calc_fluor`` & ``options.calc_vert_profiles``

- :ref:`outfiles/profiles:layer_fluorescence.dat`

if ``options.calc_ebal`` & ``options.calc_vert_profiles``

- :ref:`outfiles/profiles:leaftemp.dat`
- :ref:`outfiles/profiles:layer_h.dat`
- :ref:`outfiles/profiles:layer_le.dat`
- :ref:`outfiles/profiles:layer_a.dat`
- :ref:`outfiles/profiles:layer_NPQ.dat`
- :ref:`outfiles/profiles:layer_rn.dat`

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``etah, etau``
      - :func:`.RTMf`
    * - ``Knh, Khu``
      - :func:`.RTMz`

Fields
"""""""

Fields calculated in :func:`.RTMo` if ``options.calc_vert_profiles``

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Pn1d**
      - umol m-2 s-1
      - [60 x 1] double
      - absorbed photosynthetically active radiation (aPAR) per leaf layer
    * - **Pn1d_Cab**
      - umol m-2 s-1
      - [60 x 1] double
      - aPAR per leaf layer

Fields added in :func:`.RTMf` if ``options.calc_fluor``

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **fluorescence**
      - W m-2
      - [60 x 1] double
      - upward fluorescence per layer

Fields added in :func:`.ebal` if ``options.clc_ebal``

.. Warning::
    Averaging temperatures is not physically accurate

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **etah**
      - \-
      - [60 x 1] double
      - Fs / Fo ratio for shaded leaves
    * - **etau**
      - \-
      - [13 x 36 x 60] double
      - Fs / Fo ratio for sunlit leaves
    * - **Tchave**
      - ºC
      - double
      - mean temp shaded leaves
    * - **Tch**
      - ºC
      - [60 x 1] double
      - leaf temperature of shaded leaves, per layer
    * - **Tcu1d**
      - ºC
      - [60 x 1] double
      - leaf temperature of sunlit leaves, per layer
    * - **Tc1d**
      - ºC
      - [60 x 1] double
      - weighted average leaf temperature, per layer
    * - **Hc1d**
      - W m-2
      - [60 x 1] double
      - mean sensible heat leaves, per layer
    * - **lEc1d**
      - W m-2
      - [60 x 1] double
      - mean latent heat leaves, per layer
    * - **A1d**
      - umol m-2 s-1
      - [60 x 1] double
      - mean photosynthesis leaves, per layer
    * - **Rn1d**
      - W m-2
      - [60 x 1] double
      - net radiation per leaf layer
    * - **F_Pn1d**
      -
      - [60 x 1] double
      - mean fluorescence leaves, per layer
    * - **qE**
      -
      - [60 x 1] double
      - average NPQ = 1-(fm-fo)/(fm0-fo0), per layer
    * - **Knu**
      -
      - [13 x 36 x 60] double
      - NPQ of sunlit leaves
    * - **Knh**
      -
      - [60 x 1] double
      - NPQ of shaded leaves
