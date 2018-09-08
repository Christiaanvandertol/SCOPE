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

soil.rss_, soil.rbs_ may be changed by :func:`.calc_rssrbs` if ``options.calc_rss_rbs`` is selected

Used
"""""

Most of the values are used by :func:`.ebal`, :func:`.Soil_Inertia0`, :func:`.Soil_Inertia1`

soil.rs_ value is used by :func:`.RTMf`, :func:`.RTMo`, :func:`.RTMt_planck`, :func:`.RTMt_sb`, :func:`.RTMz`

soil.CSSOIL_ is used by :func:`.zo_and_d`

Fields
"""""""

Fields initialized in ``SCOPE.m``

.. _soil.Tsold:

:Tsold: only if ``options.soil_heat_method < 2``

    :units: ?
    :type: [12 x 2] double
    :default: 20.0

.. _soil.rs:

:refl: soil reflectance

    :units: ?
    :type: [nwl x 1] double
    :default: [2162 x 1] double

:Ts: initial soil surface temperature

    :units: ?
    :type: [2 x 1] double
    :default: ~15

|

Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

:spectrum: Spectrum number (column in the database soil_file)

    :units: \-
    :type: double
    :default: 1.0

.. _soil.rss:

:rss: soil resistance for evaporation from the pore space

    :units: s m-1
    :type: double
    :default: 500.0

:rs_thermal: broadband soil reflectance in the thermal range (1-emissivity)

    :units: \-
    :type: double
    :default: 0.06

:cs: specific heat capacity of the soil

    :units: J kg-1 K-1
    :type: double
    :default: 1180.0

:rhos: specific mass of the soil

    :units: kg m-3
    :type: double
    :default: 1800.0

.. _soil.CSSOIL:

:CSSOIL: Verhoef et al. (1997) Drag coefficient for soil *(from Aerodynamic)*

    :units: ?
    :type: double
    :default: 0.01


:lambdas: heat conductivity of the soil

    :units: J m-1 K-1
    :type: double
    :default: 1.55

.. _soil.rbs:

:rbs: soil boundary layer resistance *(from Aerodynamic)*

    :units: s m-1
    :type: double
    :default: 10.0

:SMC: volumetric soil moisture content in the root zone

    :units: ?
    :type: double
    :default: 0.25

:BSMBrightness: BSM model parameter for soil brightness

    :units: ?
    :type: double
    :default: 0.5

:BSMlat: BSM model parameter 'lat'

    :units: ?
    :type: double
    :default: 25.0

:BSMlon: BSM model parameter  'long'

    :units: ?
    :type: double
    :default: 45.0

|

Derived variables

:GAM: produced by :func:`.Soil_Inertia0` or :func:`.Soil_Inertia1` if ``options.soil_heat_method``

    :units: ?
    :type: double
    :default: ~1814. :func:`.Soil_Inertia0`


