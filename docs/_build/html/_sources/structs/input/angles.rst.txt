angles
=======

Solar and observation zenith and azimuth angles

Initialized
""""""""""""

:func:`.select_input`

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - | ``tts``
        | ``tto, psi`` -> directional_angles
      - :func:`.calc_brdf`
    * - ``tto, psi``
      - | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
    * - ``tts, tto, psi``
      - | :func:`.RTMf`
        | :func:`.RTMo`
        | :func:`.RTMz`
    * - ``tts``
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
    * - **tts**
      - deg
      - double
      - 30.0
      - solar zenith angle
    * - **tto**
      - deg
      - double
      - 0.0
      - observer zenith angle
    * - **psi**
      - deg
      - double
      - 90.0
      - | azimuthal difference between solar and observation angle
        | relative azimuth angle
