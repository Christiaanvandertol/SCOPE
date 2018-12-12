BOC_irradiance.dat
====================

BOC - bottom of canopy (61st layer)

rows - timestep

First 2162 columns: shaded fraction.

Last 2162 columns: average BOC irradiance.

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Emin_(61, :)**
      - W m-2 um-1
      - irradiance at the bottom of the canopy for the shaded fraction
    * - **\Emin_(61, :) + \Esun_(61, :) * gap.Ps(61, :)**
      - W m-2 um-1
      - average BOC irradiance (sunlit + shaded fraction)
