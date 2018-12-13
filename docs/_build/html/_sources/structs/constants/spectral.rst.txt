spectral
=========
Wavelength ranges of MODTRAN, SCOPE, PAR, fluorescence


Initialized
""""""""""""
:func:`.define_bands`

``SCOPE.m``

Variations
""""""""""""


Used
"""""

.. list-table::
    :widths: 75 25

    * - variable
      - user

    * - | ``wlS, wlF, wlT`` as index
        | ``IwlT``
      - :func:`.calc_brdf`
    * - ``wlS``
      - :func:`.ebal`
    * - ``wlE, wlF, wlP``
      - | :func:`.fluspect_B_CX`
        | :func:`.fluspect_B_CX_PSI_PSII_combined`
    * - ``wlS``
      - :func:`.create_output_files`
    * - ``wlS, wlF``
      - :func:`.initialize_output_structures`
    * - ``wlS, wlT, wlF``
      - :func:`.output_data`
    * - ``wlS, wlF, wlE, IwlF``
      - :func:`.RTMf`
    * - ``wlS, wlP, wlT, wlPAR, IwlP, IwlT``
      - :func:`.RTMo`
    * - ``wlS, wlT, IwlT``
      - :func:`.RTMt_planck`
    * - ``wlS``
      - :func:`.RTMt_sb`
    * - ``wlS, wlZ``
      - :func:`.RTMz`
    * - ``SCOPEspec``
      - :func:`.aggreg`
    * - ``wlS, wlP, wlT, wlF, IwlP, IwlT, IwlF``
      - ``SCOPE.m``

Fields
"""""""

Fields initialized in :func:`.define_bands`

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description

    * - **wlS**
      - nm
      - [1 x 2162] int
      - | 400 : 1 : 2400
        | 2500 :  100 : 15000
        | 16000 : 1000 : 50000
      - SCOPE ranges
    * - **wlP**
      - nm
      - [1 x 2001] int
      - 400 : 1 : 2400
      - PROSPECT data range
    * - **wlE**
      - nm
      - [1 x 351] int
      - 400 : 1 : 750
      - excitation in E-F matrix
    * - **wlF**
      - nm
      - [1 x 211] int
      - 640 : 1 : 850
      - chlorophyll fluorescence in E-F matrix
    * - **wlO**
      - nm
      - [1 x 2001] int
      - 400 : 1 : 2400
      - optical part (== wlP)
    * - **wlT**
      - nm
      - [1 x 161] int
      - | 2500 :  100 : 15000
        | 16000 : 1000 : 50000
      - thermal part
    * - **wlZ**
      - nm
      - [1 x 101] int
      - 500 : 1 : 600
      - xanthophyll region
    * - **wlPAR**
      - nm
      - [1 x 301] int
      - 400 : 1 : 700
      - PAR range

Fields used by :func:`.aggreg` to read MODTRAN data ``spectral.SCOPEspec``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description

    * - **SCOPEspec.nreg**
      - \-
      - int
      - 3
      - number of regions
    * - **SCOPEspec.start**
      - nm
      - [1 x SCOPEspec.nreg] int
      - [400, 2500, 16000]
      - number of regions
    * - **SCOPEspec.end**
      - nm
      - [1 x SCOPEspec.nreg] int
      - [2400, 15000, 50000]
      - number of regions
    * - **SCOPEspec.res**
      - nm
      - [1 x SCOPEspec.nreg] int
      - [1, 100, 1000]
      - number of regions

Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description

    * - **IwlP**
      - \-
      - [1 x 2001] int
      - 1 : 2001
      - index of wlP in wlS
    * - **IwlT**
      - \-
      - [1 x 161] int
      - 2002 : 2162
      - index of wlT in wlS
    * - **IwlF**
      - \-
      - [1 x 211] int
      - 241 : 451
      - index of wlF in wlS
