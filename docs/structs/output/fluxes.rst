fluxes
========
Fluxes calculated by the model (turbulent heat exchange, radiation, CO2)

Fields are initialized by :func:`.initialize_output_structures`

Calculated
""""""""""""

:func:`.ebal`

Output file
""""""""""""

- :ref:`outfiles/each:fluxes.csv`

if ``options.calc_vert_profiles`` & ``options.calc_ebal`` soil fluxes will also be recorded, because soil is the 61st layer.

- :ref:`outfiles/profiles:leaftemp.dat`
- :ref:`outfiles/profiles:layer_h.dat`
- :ref:`outfiles/profiles:layer_le.dat`
- :ref:`outfiles/profiles:layer_a.dat`
- :ref:`outfiles/profiles:layer_rn.dat`

Variations
""""""""""""

Various absorbed photosynthetically active radiations (aPAR) can be calculated by :func:`.ebal`
(if ``options.clc_ebal``) or by ``SCOPE.m``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``aPAR_Cab_eta`` -> ``rad.Femtot``
      - ``SCOPE`` (if ``options.calc_fluor``)
    * - ``aPAR_Wm2`` -> fPAR
      - :func:`.output_data`

Fields
"""""""

Fields calculated in :func:`.ebal`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Rntot**
      - W m-2
      - double
      - total net radiation
    * - **lEtot**
      - W m-2
      - double
      - total latent heat flux
    * - **Htot**
      - W m-2
      - double
      - total sensible heat
    * - **Atot**
      - umol m-2 s-1
      - double
      - total net CO2 uptake (canopy + soil)
    * - **Rnctot**
      - W m-2
      - double
      - net radiation of canopy
    * - **lEctot**
      - W m-2
      - double
      - latent heat flux of canopy
    * - **Hctot**
      - W m-2
      - double
      - sensible heat of canopy
    * - **Actot**
      - umol m-2 s-1
      - double
      - net photosynthesis of canopy
    * - **Rnstot**
      - W m-2
      - double
      - net radiation of soil
    * - **lEstot**
      - W m-2
      - double
      - latent heat flux of soil
    * - **Hstot**
      - W m-2
      - double
      - sensible heat of soil
    * - **Gtot**
      - W m-2
      - double
      - soil heat flux
    * - **Resp**
      - umol m-2 s-1
      - double
      - soil respiration rate
    * - **Au**
      - umol m-2 s-1
      - [13 x 36 x 60] double
      - sunlit leaves net CO2 assimilation
    * - **Ah**
      - umol m-2 s-1
      - [60 x 1] double
      - shaded leaves net CO2 assimilation

Fields added by :func:`.ebal`  (if ``options.calc_ebal == 1``) or by ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **aPAR**
      - umol m-2 s-1
      - double
      - absorbed PAR by leaves
    * - **aPAR_Cab**
      - umol m-2 s-1
      - double
      - absorbed PAR by chlorophylls a, b
    * - **aPAR_Wm2**
      - W m-2
      - double
      - absorbed PAR
    * - **aPAR_Cab_eta**
      - umol m-2 s-1
      - double
      - green ePAR * relative fluorescence emission efficiency
