canopy
=======
Canopy parameters, such as leaf area index and leaf inclination distribution function

Initialized
""""""""""""
``SCOPE.m``

:func:`.select_input`

Variations
""""""""""""

canopy.lidf_ can be read from ``LIDF_file`` if its name is provided in the ``filenames.csv``.

.. note::
    ``LIDF_file`` must be located in ``./input/leafangles`` (:ref:`directories/input:leafangles`) and have **3 header lines**.

canopy.zo_, canopy.d_ may be calculated by :func:`.zo_and_d` if ``options.calc_zo`` is selected

canopy.hc_ may be set in :func:`.load_timeseries`

.. warning::
    never change the angles in canopy.litab_ unless :func:`.leafangles` ('ladgen') is also adapted

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user


    * - ``nlayers, nlincl, nlazi, lidf``
      - :func:`.meanleaf`
    * - ``CR, CD1, Psicor, LAI, hc``
      - :func:`.zo_and_d`
    * - ``LAI, hc, zo, d``
      - :func:`.load_timeseries`
    * - ``LIDFa, LIDFb``
      - :func:`.leafangles`
    * - | ``nlayers, kV, xl, LAI``
        | ``LAI, rwc, zo, d, hc, leafwidth, Cd`` -> :ref:`structs/internal/resist_in:Resist_in`
      - :func:`.ebal`
    * - ``nlayers, lidf, litab, lazitab, LAI``
      - | :func:`.RTMf`
        | :func:`.RTMo`
        | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
        | :func:`.RTMz`
    * - ``nlincl, nlazi, x, hot``
      - :func:`.RTMo`
    * - ``x, nlayers, LAI``
      - ``SCOPE.m``


The meaning of LIDFa and LIDFb in relation to the leaf angle distribution (LAD) is presented in

.. figure:: ../../images/LAD.png
    :align: center


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
    * - **nlayers**
      - \-
      - int
      - 60
      - the number of layers in a canopy
    * - **x**
      - \-
      - [60 x 1] double
      - | (0 : -1]
        | equally spaced vector
      - | levels in canopy except for the top:
        | ``bottom = -1``,
        | ``top = -1/canopy.nlayers``
        | in fact length == canopy.nlayers + 1
    * - **xl**
      - \-
      - [61 x 1] double
      - | [0 : -1]
        | equally spaced vector
      - | levels in canopy and the top
        | [0, canopy.x]
        | in fact length == canopy.nlayers + 1
    * - **nlincl**
      - \-
      - int
      - 13
      - number of leaf inclinations
    * - **nlazi**
      - \-
      - int
      - 36
      - number of leaf azimuth angles
    * - .. _canopy.litab:

        **litab**
      - deg
      - [13 x 1] double
      - | [5 : 89]
        | *non-equally* spaced vector
      - SAIL leaf inclination angles
    * - **lazitab**
      - \-
      - [1 x 36] double
      - | [5 : 355]
        | equally spaced vector
      - leaf azimuth angles relative to the sun
    * - .. _canopy.lidf:

        **lidf**
      - ?
      - [13 x 1] double
      - :func:`.leafangles`
      - leaf inclination distribution function


Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **LAI**
      - m2 m-2
      - double
      - 3.0
      - Leaf area index
    * - .. _canopy.hc:

        **hc**
      - m
      - double
      - 2.0
      - vegetation height
    * - **LIDFa**
      - \-
      - double
      - -0.35
      - leaf inclination
    * - **LIDFb**
      - \-
      - double
      - -0.15
      - variation in leaf inclination
    * - **leafwidth**
      - m
      - double
      - 0.1
      - leaf width
    * - **rb**
      - s m-1
      - double
      - 10.0
      - leaf boundary resistance
    * - **Cd**
      - ?
      - double
      - 0.3
      - leaf drag coefficient
    * - **CR**
      - ?
      - double
      - 0.35
      - Verhoef et al. (1997)  Drag coefficient for isolated tree
    * - **CD1**
      - ?
      - double
      - 20.6
      - Verhoef et al. (1997)  fitting parameter
    * - **Psicor**
      - ?
      - double
      - 0.2
      - Verhoef et al. (1997)  Roughness layer correction
    * - **rwc**
      - s m-1
      - double
      - 0.0
      - within canopy layer resistance
    * - **kV**
      - ?
      - double
      - 0.6396
      - extinction coefficient for ``Vcmax`` in the vertical (maximum at the top). 0 for uniform ``Vcmax``
    * - .. _canopy.zo:

        **zo**
      - m
      - double
      - 0.246
      - roughness length for momentum of the canopy
    * - .. _canopy.d:

        **d**
      - m
      - double
      - 1.34
      - displacement height
    * - **hot**
      - ?
      - double
      - 0.05
      - hotspot parameter ``canopy.leafwidth / canopy.hc``
