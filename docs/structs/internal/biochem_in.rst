biochem_in
===========

.. Note::
    This is an internal struct of :func:`.ebal`. It is only available within :func:`.ebal` in debug mode.

Input of the biochemical routine :func:`.biochemical` or :func:`.biochemical_MD12` for **leaf** photosynthesis and fluorescence.

Initialized
""""""""""""

:func:`.ebal`

Variations
""""""""""""

For sunlit leaves size is [13 x 36 x 60] for shaded [60 x 1]



Fields
"""""""

Fields initialized in :func:`.ebal`


.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Fluorescence_model**
      - \-
      - bool
      - ``options.Fluorescence_model``
    * - **Type**
      - \-
      - char
      - photosynthesis type C3 or C4
    * - **p**
      - hPa
      - double
      - atmospheric pressure
    * - **m**
      -
      - double
      - Ball-Berry stomatal conductance parameter
    * - **O**
      - per mille
      - double
      - atmospheric O2 concentration
    * - **Rdparam**
      - \-
      - double
      - fraction of respiration
    * - **T**
      - ºC
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - leaf temperature per canopy layer
    * - **eb**
      - hPa
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - leaf water (ea) per canopy layer
    * - **Cs**
      - ppm
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - leaf CO2 concentration per canopy layer
    * - **Vcmo**
      - umol m-2 s-1
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - maximum carboxylation rate per canopy layer
    * - **Q**
      - W m-2
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - absorbed photosynthetically active radiation (PAR) by chlorophylls per canopy layer
    * - **A**
      - umol m-2 s-1
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - photosynthesis (CO2 assimilation rate) per canopy layer

These parameters are specific for :func:`.biochemical`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **tempcor**
      - \-
      - double
      - ``options.apply_T_corr``
    * - **Tparams**
      - K
      - [5 x 1] double
      - the temperature response of fluorescence
    * - **stressfactor**
      - \-
      - double
      - stress factor to reduce Vcmax

These parameters are specific for :func:`.biochemical_MD12`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Tyear**
      - ºC
      - double
      - mean annual temperature
    * - **beta**
      - \-
      - double
      - fraction of photons partitioned to PSII
    * - **kNPQs**
      - s-1
      - double
      - rate constant of sustained thermal dissipation
    * - **qLs**
      - \-
      - double
      - fraction of functional reaction centres