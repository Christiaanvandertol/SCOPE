gap
======

Probabilities of levels being observed, illuminated per layer

.. Note::
    This is an optional struct that requires ``options.calc_vert_profiles`` to be output, however it will be calculated anyway

Calculated
""""""""""""

:func:`.RTMo`

Output file
""""""""""""""

if ``options.calc_vert_profiles``

:ref:`outfiles/profiles:gap.dat`

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user

    * - ``Ps``
      - | :func:`.ebal`
        | ``SCOPE.m``
    * - ``Pso``
      - :func:`.RTMo`
    * - ``Ps, Po, Pso``
      - | :func:`.RTMf`
        | :func:`.RTMz`
        | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
    * - ``K``
      - | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`


Fields
"""""""

Fields initialized in :func:`.RTMo`

.. list-table::
    :widths: 10 20 20 50

    * - variable
      - units
      - type
      - description
    * - **Ps**
      - \-
      - [61 x 1] double
      - fraction of sunlit leaves per layer
    * - **Po**
      - \-
      - [61 x 1] double
      - fraction of observed leaves per layer
    * - **Pso**
      - \-
      - [61 x 1] double
      - fraction of sunlit and (at the same time) observed per layer
    * - **K**
      - \-
      - double
      - extinction coefficient in direction of observer integrated over leaf angles
    * - **k**
      - \-
      - double
      - extinction coefficient in direction of sun integrated over leaf angles
