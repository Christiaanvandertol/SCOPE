``options.calc_vert_profiles``
=================================

gap.dat
---------

.. Note:: ``options.calc_vert_profiles``

rows - time (simulation number)

columns - [Ps Po Pso] => 61 * 3 columns

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Ps**
      - \-
      - fraction of sunlit leaves per layer
    * - **Po**
      - \-
      - fraction of observed leaves per layer
    * - **Pso**
      - \-
      - fraction of sunlit and (at the same time) observed per layer


layer_a.dat
-------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - photosynthesis per layer, total soil respiration (60 + 1)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **A1d**
      - umol m-2 s-1
      - mean photosynthesis leaves, per layer
    * - **Resp**
      - umol m-2 s-1
      - soil respiration rate

layer_aPAR.dat
------------------

.. Note:: ``options.calc_vert_profiles``

rows - time (simulation number)

columns - absorbed photosynthetically active radiation (aPAR)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Pn1d**
      - umol m-2 s-1
      - aPAR per leaf layer

layer_aPAR_Cab.dat
----------------------

.. Note:: ``options.calc_vert_profiles``

rows - time (simulation number)

columns - absorbed photosynthetically active radiation (aPAR) by chlorophylls (Cab) per leaf layer

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Pn1d_Cab**
      - umol m-2 s-1
      - aPAR by Cab per leaf layer

layer_fluorescence.dat
---------------------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_fluor``

rows - time (simulation number)

columns - upward fluorescence per layer

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **fluorescence**
      - W m-2
      - upward fluorescence per layer

layer_h.dat
---------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - sensible heat flux per layer, total sensible heat of soil (60 + 1)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Hc1d**
      - W m-2
      - mean sensible heat leaves, per layer
    * - **Hstot**
      - W m-2
      - sensible heat of soil


layer_le.dat
--------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - latent heat flux per layer, total latent heat of soil (60 + 1)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **lEc1d**
      - W m-2
      - mean latent heat leaves, per layer
    * - **lEstot**
      - W m-2
      - latent heat of soil


layer_NPQ.dat
----------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - average NPQ = 1-(fm-fo)/(fm0-fo0), per layer (60)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **qE**
      -
      - average NPQ = 1-(fm-fo)/(fm0-fo0), per layer


layer_rn.dat
--------------

.. Note:: ``options.calc_vert_profiles`` & ``options.calc_ebal``

rows - time (simulation number)

columns - net radiation per leaf layer, total net radiation of soil (60 + 1)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Rn1d**
      - W m-2
      - net radiation per leaf layer
    * - **Rnstot**
      - W m-2
      - net radiation of soil


leaftemp.dat
----------------

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
