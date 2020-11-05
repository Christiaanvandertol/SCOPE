directional
================

.. Note::
    This is an optional struct that requires ``options.calc_directional`` & ``options.calc_ebal``
    to be output and calculated

Calculated
""""""""""""

:func:`.calc_brdf`

Variations
""""""""""""

``psi, tto`` initialized in ``SCOPE.m`` are extended in :func:`.calc_brdf`

``options.calc_planck`` enables ``Lot_`` calculation (disables of ``Lot`` and ``BrightnessT``)

``options.calc_fluor`` enables ``LoF_`` calculation

Output file
""""""""""""

- :ref:`outfiles/directional:Directional/Angles (SunAngle x.xx degrees).dat`

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``noa``
      - :func:`.calc_brdf`

Fields
"""""""

Fields initialized in ``SCOPE.m`` (read from :ref:`directories/input:directional`)

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **noa**
      - \-
      - double
      - number of observation angles
    * - **psi**
      - deg
      - [333 x 1] double
      - observation relative azimuth angles
    * - **tto**
      - deg
      - [333 x 1] double
      - observation zenith angles


Fields calculated in :func:`.calc_brdf`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **psi**
      - deg
      - [364 x 1] double
      - observation relative azimuth angles
    * - **tto**
      - deg
      - [364 x 1] double
      - observation zenith angles
    * - **brdf_**
      - \-
      - [2162 x 364] double
      - bidirectional reflectance distribution function
    * - **Eoutte**
      -
      - [1 x 364] double
      -
    * - **Lot**
      -
      - [161 x 364] double
      -
    * - **BrightnessT**
      - K
      - [1 x 364] double
      - brightness temperature
    * - ``options.calc_planck``
      -
      -
      -
    * - **Lot_**
      - W m-2 um-1 sr-1
      - [161 x 364] double
      - outgoing thermal radiation in observation direction
    * - ``options.calc_fluor``
      -
      -
      -
    * - **LoF_**
      - W m-2 um-1 sr-1
      - [211 x 364] double
      - outgoing fluorescence radiation in observation direction
