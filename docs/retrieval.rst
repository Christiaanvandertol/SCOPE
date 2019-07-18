Retrieval (Model inversion)
=============================

.. warning::
    Github of the project: https://github.com/Prikaziuk/retrieval_rtmo

.. contents::

.. toctree::
    :maxdepth: 1

    retrieval/hyperspectral
    retrieval/multispectral
    retrieval/plots
    retrieval/sensors
    retrieval/tricks


Highlights
-----------

1. Numerical optimization of RTMo (optical) module of SCOPE model
    * Top of canopy (TOC) reflectance
2. Sensors
    * ground (ASD, FLOX)
    * airborne (HyPlant)
    * satellite (Sentinel-2 MSI, Sentinel-3 OLCI)
    * custom - just add yourself into ``input/sensors.xlsx``
3. Available options
    * Hyperspectral instruments:
        * data format - text
        * parallel computing (parfor)
        * time-series (angles, radiation)
        * calculation of observation geometry from coordinates and time
        * plotting of validation data (1 : 1 plots)
        * output - .xlsx or .csv
    * Multispectral instruments (Satellites):
        * data format - NetCDF4
        * parallel computing (parfor)
        * retrieval on an image subset (K x K pixels)
        * output - .xlsx + NetCDF4

Definition
-----------

Typically you provide vegetation parameters and SCOPE gives you expected reflectance spectra.
For retrieval you provide the reflectance spectra and SCOPE gives you expected vegetation parameters, resulting in similar spectra.

That is why it is called model inversion. We somehow want to get the input (parameters) from the output (spectra).

There are two basic inversion approaches:
    * look-up table
        - spectra are produced beforehand in the desired ranges of values
        - for example you change chlorophyll content (Cab) from 0 to 100 ug cm-2 in steps of 5.
        - **advantage** fast
        - **disadvantage** coarse: if you tune 6 parameters giving only 10 steps for each you
          end up in a look-up table with 1.000.000 simulations
    * numerical optimization
        - spectra are calculated on the spot (for each measurement during the run)
        - algorithm is minimizing the cost-function (difference between modelled and measured data)
        - **advantage** sharp values
        - **disadvantage** slow, local minimum trap


Directory structure
---------------------

.. code-block:: none
    :emphasize-lines: 12, 13, 14, 22

    demonstrator-master.zip
    ├── src
    │   ├── +helpers
    │   ├── +io
    │   ├── +models
    │   ├── +plot
    │   ├── +sat
    │   ├── +to_sensor
    │   ├── +ts
    │   ├── COST_4SAIL_common.m
    │   ├── fit_spectra.m
    │   ├── Input_data.xlsx
    │   ├── main.m
    │   ├── main_sat.m
    │   └── propagate_uncertainty.m
    │
    ├── input
    │   ├── fluspect_data
    │   ├── radiationdata
    │   ├── soil_spectrum
    │   ├── PC_flu.xlsx
    │   └── sensors.xlsx
    │
    ├── measured
    │   ├── airborne
    │   ├── canopy
    │   ├── leaf
    │   ├── synthetic
    │   └── Sentinel-2_sample.nc
    │
    └── output
        ├── sample_S2
        └── sample_synthetic


Acknowledgement
-----------------

The project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 721995.

.. figure:: images/trustee_logo.png
    :align: center
    :target: http://www.trusteenetwork.eu/