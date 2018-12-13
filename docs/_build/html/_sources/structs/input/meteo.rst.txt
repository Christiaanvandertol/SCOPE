meteo
======
Meteorological variables

Initialized
""""""""""""

:func:`.select_input`

Variations
""""""""""""

:func:`.ebal` uses ``max(meteo.u, 0.2)``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``z``
      - :func:`.load_timeseries`
    * - | ``Ta, ea, Ca, p``
        | ``u, z`` -> :ref:`structs/internal/resist_in:Resist_in`
        | ``Oa`` -> :ref:`structs/internal/biochem_in:Biochem_in`
      - :func:`.ebal`
    * - ``Rin, Rli``
      - :func:`.RTMo`
    * - ``Rin, Rli``
      - :func:`.output_data`


Fields
"""""""

Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **z**
      - m
      - double
      - 10.0
      - measurement height of meteorological data
    * - **Rin**
      - W m-2
      - double
      - 600.0
      - broadband incoming shortwave radiation (0.4-2.5 um)
    * - **Ta**
      - ºC
      - double
      - 20.0
      - air temperature
    * - **Rli**
      - W m-2
      - double
      - 300.0
      - broadband incoming longwave radiation (2.5-50 um)
    * - **p**
      - hPa
      - double
      - 970.0
      - air pressure
    * - **ea**
      - hPa
      - double
      - 15.0
      - atmospheric vapour pressure
    * - **u**
      - m s-1
      - double
      - 2.0
      - wind speed at height z
    * - **Ca**
      - ppm
      - double
      - 380.0
      - atmospheric CO2 concentration
    * - **Oa**
      - per mille
      - double
      - 209.0
      - atmospheric O2 concentration
