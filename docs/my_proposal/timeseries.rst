Time series
================

``options.simulation== 1``

Definition
''''''''''''
A set of variables recorded together with time.

SCOPE
''''''''

Set option :ref:`options:``simulation``` to 1.

Provide ``startDOY`` and ``endDOY`` as time stamps in berkeley format: (FLUXNET standard since 2015) - YYYYMMDD[HHMMSS].
The data read from files in *"./input/dataset X"* will be subset to your startDOY, endDOY.

SCOPE will make runs for every time stamp in your data.

It is possible to interpolate variable that are measured less frequently than meteorological data (:ref:`directories/input:vegetation_retrieved_csv`).
This file can be created from 'Output' sheet of retrieval :ref:`retrieval/hyperspectral:Output` by adding a time stamp column.

.. figure:: ../images/ts.png
    :scale: 50%
    :align: center

    Results of time series run from "./input/:ref:`directories/input:dataset for_verification`".
    Each point - single SCOPE run.
    Left column - inputs, right - outputs from :ref:`outfiles/each:fluxes.csv`.
    Blue - standard 'verificationdata' dataset, red - dataset where Vcmax and Cab values are interpolated.
    As you see Cab=80 for standard and Cab changes with time (from 37.7 to 39.2) with interpolation.

.. warning::
    Sometimes time series produce Inf, -Inf, NaN in the output and error or warning of energy balance closure.
    We recommend checking energy balance closure for each time step.
