Directional
============

.. Note::
    This is an optional output that requires ``options.calc_directional`` & ``options.calc_ebal``
    However, the folder will always be created

Angles (SunAngle x.xx degrees).dat
""""""""""""""""""""""""""""""""""""""

Contains the directions.

   * The 1st row gives the observation zenith angles
   * The 2nd row gives the observation azimuth angles
   * The 3rd row gives the solar zenith angles (constant from ``input_data.xlsx``)

columns - combination number (a set of direction used for simulation)

Columns in the output files correspond to the columns in Angles

BRDF (SunAngle x.xx degrees).dat
""""""""""""""""""""""""""""""""""""""

   * The 1st column gives the wl values corresponding to the BRDF values
   * Other columns give the BRDF values corresponding to the directions given by observation zenith angles (first column in the Angles file)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **wlS**
      - um
      - full wl range SCOPE
    * - **brdf_**
      - \-
      - bidirectional reflectance distribution function

Temperatures (SunAngle x.xx degrees).dat
""""""""""""""""""""""""""""""""""""""""""

   * The 1st column gives the wl values corresponding to the brightness temperature values (except for broadband)
   * Other columns give the brightness temperature (BT) values corresponding to the directions given by a column in the Angles file

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **BrightnessT**
      - K
      - brightness temperature
    * - ``options.calc_planck``
      -
      -
    * - **wlT**
      - um
      - thermal wl range SCOPE
    * - **Lot_**
      -
      -


Fluorescence (SunAngle x.xx degrees).dat
""""""""""""""""""""""""""""""""""""""""""

if ``options.calc_fluor``

   * The 1st column gives the wl values corresponding to the brightness temperature values (except for broadband)
   * Other columns give the fluorescence corresponding to the directions given by a column in the Angles file

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **wlF**
      - um
      - fluorescence wl range SCOPE
    * - **LoF_**
      -
      -

read me.txt
""""""""""""""""""""""""""""""""""""""""""

