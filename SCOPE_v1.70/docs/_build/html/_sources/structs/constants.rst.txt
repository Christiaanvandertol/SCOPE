constants
==========
Physical constants

Initialized
""""""""""""
:func:`.define_constants`

Variations
""""""""""""


Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user

    * - ``sigmaSB``
      - :func:`.calc_brdf`
    * - ``MH2O, Mair, rhoa, cp, g, kappa, sigmaSB``
      - :func:`.ebal`
    * - ``kappa``
      - | :func:`.zo_and_d`
        | :func:`.resistances`
    * - ``rhoa, cp, MH2O, R``
      - :func:`.heatfluxes`
    * - ``deg2rad, Mair, MCO2, rhoa``
      - :func:`.load_timeseries`
    * - ``sigmaSB, C2K``
      - | :func:`.Brightness_T`
        | :func:`.RTMt_sb`
    * - ``deg2rad``
      - | :func:`.RTMf`
        | :func:`.RTMz`
        | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
    * - ``A, h, c``
      - :func:`.RTMo`


Fields
"""""""

Fields initialized in :func:`.define_constants`

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description

    * - **A**
      - mol-1
      - double
      - 6.02214E23
      - Constant of Avogadro
    * - **h**
      - J s
      - double
      - 6.6262E-34
      - Planck's constant
    * - **c**
      - m s-1
      - double
      - 299792458
      - Speed of light
    * - **cp**
      - J kg-1 K-1
      - double
      - 1004
      - Specific heat of dry air
    * - **R**
      - J mol-1K-1
      - double
      - 8.314
      - Molar gas constant
    * - **rhoa**
      - kg m-3
      - double
      - 1.2047
      - Specific mass of air
    * - **g**
      - m s-2
      - double
      - 9.81
      - Gravity acceleration
    * - **kappa**
      - \-
      - double
      - 0.4
      - Von Karman constant
    * - **MH2O**
      - g mol-1
      - double
      - 18
      - Molecular mass of water
    * - **Mair**
      - g mol-1
      - double
      - 28.96
      - Molecular mass of dry air
    * - **MCO2**
      - g mol-1
      - double
      - 44
      - Molecular mass of carbon dioxide
    * - **sigmaSB**
      - W m-2 K-4
      - double
      - 5.67E-8
      - Stefan Boltzman constant
    * - **deg2rad**
      - rad
      - double
      - pi/180
      - Conversion from deg to rad
    * - **C2K**
      - K
      - double
      - 273.15
      - Melting point of water
