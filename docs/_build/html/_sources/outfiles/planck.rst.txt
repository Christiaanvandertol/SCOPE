``options.calc_ebal`` & ``options.calc_planck``
================================================

.. contents::

spectrum_hemis_thermal.dat
----------------------------

.. Note:: ``options.calc_ebal`` & ``options.calc_planck``

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Eoutte_**
      - W m-2 um-1
      - hemispherical outgoing thermal radiation

spectrum_obsdir_BlackBody.dat
--------------------------------

.. Note:: ``options.calc_ebal``

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **LotBB_**
      - W m-2 um-1 sr-1
      - thermal BlackBody emission spectrum in observation direction

spectrum_obsdir_thermal.dat
-----------------------------

.. Note:: ``options.calc_ebal`` & ``options.calc_planck``

rows - time (simulation number)

columns - wl number (2162)

.. list-table::
    :widths: 20 20 60

    * - variable
      - units
      - description
    * - **Lot_**
      - W m-2 um-1 sr-1
      - outgoing thermal radiation in observation direction
