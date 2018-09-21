xyt
====
Geographical location and time of the project

Initialized
""""""""""""

:func:`.select_input`
:func:`.load_timeseries`

Variations
""""""""""""

``SCOPE`` overwrites xyt.t_, xyt.year_ if ``options.simulation != 1``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``t, timezn, LON, LAT``
      - :func:`.load_timeseries`
    * - ``t``
      - :func:`.ebal`
    * - ``t, startDOY, endDOY``
      - ``SCOPE``
    * - ``year, t``
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
    * - **startDOY**
      - \-
      - double
      - 169.0
      - Julian day (decimal) of start of simulations
    * - **endDOY**
      - \-
      - double
      - 170.0
      - Julian day (decimal) of end of simulations
    * - **LAT**
      - decimal deg
      - double
      - 52.25
      - Latitude
    * - **LON**
      - decimal deg
      - double
      - 5.69
      - Longitude
    * - **timezn**
      - hours
      - double
      - 1.0
      - east of Greenwich

Fields initialized in :func:`.load_timeseries`

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - .. _xyt.t:

        **t**
      - \-
      - [n x 1] double
      - [214 x 1]
      - Julian day (decimal)
    * - .. _xyt.year:

        **year**
      - \-
      - [m x 1] int
      - [216 x 1]
      - Calendar year