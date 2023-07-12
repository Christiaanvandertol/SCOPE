soil
=====
Soil properties (such as soil moisture, emissivity and the reflectance spectrum)

Initialized
""""""""""""
:func:`.select_input`

``SCOPE.m``

Variations
""""""""""""
soil.Tsold_ may be changed by :func:`.ebal` if ``options.soil_heat_method < 2`` (default case)

soil.rss_, soil.rbs_ may be calculated by :func:`.calc_rssrbs` if ``options.calc_rss_rbs`` is selected

soil.GAM_ produced by :func:`.Soil_Inertia0` or :func:`.Soil_Inertia1` if ``options.soil_heat_method``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user

    * - ``spectrum, rs_thermal``
      - ``SCOPE.m``
    * - ``cs, rhos, lambdas``
      - :func:`.Soil_Inertia0`
    * - ``SMC``
      - | :func:`.Soil_Inertia1`
        | :func:`.BSM`
        | :func:`.calc_rssrbs`
    * - ``CSSOIL``
      - :func:`.zo_and_d`
    * - ``refl``
      - | :func:`.RTMf`
        | :func:`.RTMo`
        | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
        | :func:`.RTMz`
    * - | ``Ts, Tsold, GAM, rss``
        | ``rbs, rss`` -> :ref:`structs/internal/resist_in:Resist_in`
      - :func:`.ebal`


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

    * - **spectrum**
      - \-
      - int
      - 1
      - spectrum number (column in the database soil_file)
    * - .. _soil.rss:

        **rss**
      - s m-1
      - double
      - 500.0
      - soil resistance for evaporation from the pore space
    * - **rs_thermal**
      - \-
      - double
      - 0.06
      - broadband soil reflectance in the thermal range ``1 - emissivity``
    * - **cs**
      - J kg-1 K-1
      - double
      - 1180.0
      - specific heat capacity of the soil
    * - **rhos**
      - kg m-3
      - double
      - 1800.0
      - specific mass of the soil
    * - **CSSOIL**
      - \-
      - double
      - 0.01
      - Drag coefficient for soil Verhoef et al. (1997) *(from Aerodynamic)*
    * - **lambdas**
      - J m-1 K-1
      - double
      - 1.55
      - heat conductivity of the soil
    * - .. _soil.rbs:

        **rbs**
      - s m-1
      - double
      - 10.0
      - soil boundary layer resistance *(from Aerodynamic)*
    * - **SMC**
      - \-
      - double
      - 0.25
      - volumetric soil moisture content in the root zone
    * - **BSMBrightness**
      - \-
      - double
      - 0.5
      - BSM model parameter for soil brightness
    * - **BSMlat**
      - \-
      - double
      - 25.0
      - BSM model parameter 'lat'
    * - **BSMlon**
      - \-
      - double
      - 45.0
      - BSM model parameter  'long'

Derived variables


.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - .. _soil.GAM:

        **GAM**
      - J m-2 s-0.5 K-1
      - double
      - ~1814.4 :func:`.Soil_Inertia0`,:func:`.Soil_Inertia1`
      - soil thermal inertia

Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description

    * - .. _soil.Tsold:

        **Tsold**
      - ºC?
      - [12 x 2] double
      - 20.0
      - only if ``options.soil_heat_method < 2``
    * - **refl**
      - \-
      - [2162 x 1] double
      - [2162 x 1] double
      - | soil reflectance
        | in fact length == nwl
    * - **Ts**
      - ºC?
      - [2 x 1] double
      - [~15; ~15]
      - initial soil surface temperature