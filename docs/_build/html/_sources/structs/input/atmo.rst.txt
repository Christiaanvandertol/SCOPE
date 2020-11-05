atmo
======

Atmospheric transfer functions from standard FLEX atmospheres

Initialized
""""""""""""

``SCOPE.m`` loaded from ./input/:ref:`directories/input:radiationdata` and aggregated by :func:`.aggreg`.

Filename is specified on ``filenames`` sheet, ``atmos_file`` cell of ``input_data.xlsx``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``M, Ta``
      - :func:`.RTMo`


Fields
"""""""

Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **M**
      - \-
      - [2162 x 6] double
      - from ``FLEX-S3_std.atm``
      - atmospheric transmittance functions T1, T3, T4, T5, T12, T16
    * - **Ta**
      - ºC
      - double
      - 20.0 (== meteo.Ta)
      - air temperature
