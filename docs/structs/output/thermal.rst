thermal
========
Leaf, soil and air temperatures


Fields are initialized by :func:`.initialize_output_structures`

Calculated
""""""""""""

:func:`.ebal`

Output file
""""""""""""

- :ref:`outfiles/each:fluxes.csv`
- :ref:`outfiles/each:aerodyn.dat`

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``Tcu, Tch, Ts(1)``
      - | :func:`.calc_brdf` (if ``options.calc_planck``) -> :func:`.RTMt_planck`
        | :func:`.calc_brdf` (else) -> :func:`.RTMt_sb`
    * - ``Tcu, Tch, Ts(1), Ts(2)``
      - ``SCOPE.m`` (if ``options.calc_planck``) -> :func:`.RTMt_planck`


Fields
"""""""

Fields added in :func:`.ebal`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Tcave**
      - ºC
      - double
      - canopy weighted average temperature
    * - **Tsave**
      - ºC
      - double
      - soil weighted average temperature
    * - **raa**
      - s m-1
      - double
      - total aerodynamic resistance above canopy
    * - **rawc**
      - s m-1
      - double
      - canopy total aerodynamic resistance below canopy
    * - **raws**
      - s m-1
      - double
      - soil total aerodynamic resistance below canopy
    * - **ustar**
      - m s-1
      - double
      - friction velocity
    * - **Ts**
      - ºC
      - [2 x 1] double
      - sunlit and shaded soil temperature
    * - **Ta**
      - ºC
      - double
      - air temperature as in input
    * - **Tcu**
      - ºC
      - [13 x 36 x 60] double
      - sunlit leaves canopy temperature
    * - **Tch**
      - ºC
      - [60 x 1] double
      - shaded leaves canopy temperature

