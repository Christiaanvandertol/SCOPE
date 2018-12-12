leaftemp.dat
========================

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - leaf temperatures per layer (60 * 3) leaf temperature of sunlit leaves, shaded leaves, and weighted average leaf temperature per layer

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Tcu1d**
      - ºC
      - leaf temperature of sunlit leaves, per layer
    * - **Tch**
      - ºC
      - leaf temperature of shaded leaves, per layer
    * - **Tc1d**
      - ºC
      - weighted average leaf temperature, per layer
