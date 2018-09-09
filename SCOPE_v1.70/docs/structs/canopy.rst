Canopy
=======
Canopy parameters, such as leaf area index and leaf inclination distribution function

Initialized
""""""""""""
``SCOPE.m``

:func:`.select_input`

Variations
""""""""""""

canopy.lidf_ can be read from ``LIDF_file`` if its name is provided in the ``filenames`` sheet of ``input_data.xlsx``

.. note::
    ``LIDF_file`` must be located in ``/data/input/leafangles`` (:ref:`leafangles`) and have **3 header lines**.

canopy.zo_, canopy.d_ may be calculated by :func:`.zo_and_d` if ``options.calc_zo`` is selected

canopy.hc_ may be set in :func:`.load_timeseries`

Used
"""""

.. list-table::
    :header-rows: 0

    * - ``nlayers, kV, xl, LAI`` ``LAI, rwc, zo, d, hc, leafwidth, Cd``
      - :func:`.ebal`

    * - ``nlayers, LAI, litab, lazitab, lidf``
      - | :func:`.RTMf` :func:`.RTMt_planck` :func:`.RTMt_sb` :func:`.RTMz`
        | ``nlayers, LAI, litab, lazitab, lidf`` values is used by :func:`.RTMf`, :func:`.RTMo` (also ``x, hot, nlincl, nlazi``), :func:`.RTMt_planck`,
        | :func:`.RTMt_sb`, :func:`.RTMz`

``nlayers, nlincl, nlazi and lidf`` are used by :func:`.meanleaf`

Most of the values are used by :func:`.ebal` for calculations (``nlayers, kV, xl, LAI``) and to initialize the :ref:`Resist_in` structure (``LAI, rwc, zo, d, hc, leafwidth, Cd``)

``nlayers, LAI, litab, lazitab, lidf`` values is used by :func:`.RTMf`, :func:`.RTMo` (also ``x, hot, nlincl, nlazi``), :func:`.RTMt_planck`, :func:`.RTMt_sb`, :func:`.RTMz`

``nlayers, nlincl, nlazi and lidf`` are used by :func:`.meanleaf`

``LIDFa, LIDFb`` are used by :func:`.leafangles`


``CR, CD1, Psicor, LAI, hc`` are used by :func:`.zo_and_d`

``LAI, hc, zo, do`` are used by :func:`.load_timeseries`


Fields
"""""""

Fields initialized in ``SCOPE.m``

:nlayers: the number of layers in a canopy

    :units: \-
    :type: int
    :default: 60

:x: levels in canopy except for the top: *bottom = -1, top = -1/canopy.nlayers*

    :units: \-
    :type: [canopy.nlayers x 1] double
    :default: [60 x 1] equally spaced vector (0, -1]

:xl: levels in canopy (canopy.x) and the top = 0

    :units: \-
    :type: [(canopy.nlayers + 1) x 1] double
    :default: [61 x 1] equally spaced vector [0, -1]

:nlincl: number of leaf inclinations

    :units: \-
    :type: int
    :default: 13

:nlazi: number of leaf azimuth angles

    :units: \-
    :type: int
    :default: 36

:litab: SAIL leaf inclination angles

    :units: deg
    :type: [13 x 1] double
    :default: non-equally spaced vector [5, 89]

 .. warning::

    never change the angles in canopy.litab unless :func:`.leafangles` ('ladgen') is also adapted

:lazitab: leaf azimuth angles relative to the sun

    :units: deg
    :type: [1 x 36] double
    :default: equally spaced vector [5, 355]

.. _canopy.lidf:

:lidf: leaf inclination distribution function

    :units: ?
    :type: [13 x 1] double
    :default: calculated by :func:`.leafangles`

Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

:LAI: Leaf area index

    :units: m2 m-2
    :type: double
    :default: 3.0

.. _canopy.hc:

:hc: vegetation height

    :units: m
    :type: double
    :default: 2.0

:LIDFa: leaf inclination

    :units: \-
    :type: double
    :default: -0.35

:LIDFb: variation in leaf inclination

    :units: \-
    :type: double
    :default: -0.15

:leafwidth: leaf width

    :units: m
    :type: double
    :default: 0.1

:rb: leaf boundary resistance

    :units: s m-1
    :type: double
    :default: 10.0

:Cd: leaf drag coefficient

    :units: nan
    :type: double
    :default: 0.3

:CR: Verhoef et al. (1997)  Drag coefficient for isolated tree: Andrieu1997?

    :units: ?
    :type: double
    :default: 0.35

:CD1: Verhoef et al. (1997)  fitting parameter

    :units: ?
    :type: double
    :default: 20.6

:Psicor: Verhoef et al. (1997)  Roughness layer correction

    :units: ?
    :type: double
    :default: 0.2

:rwc: within canopy layer resistance

    :units: s m-1
    :type: double
    :default: 0.0

:kV: extinction coefficient for ``Vcmax`` in the vertical (maximum at the top). 0 for uniform ``Vcmax``

    :units: nan
    :type: double
    :default: 0.6396

.. _canopy.zo:

:zo: roughness length for momentum of the canopy

    :units: m
    :type: double
    :default: 0.246

.. _canopy.d:

:d: displacement height

    :units: m
    :type: double
    :default: 1.34

:hot: hotspot parameter ``canopy.leafwidth / canopy.hc``

    :units: \-
    :type: double
    :default: 0.05