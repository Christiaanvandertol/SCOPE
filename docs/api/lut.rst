+lut
==========

These files were created for SCOPE per pixel implementation within the
`TRUSTEE project <http://www.trusteenetwork.eu/>`_ (H2020, TRuStEE, MSCA-ITN-2016 No 721995) tutorials.
Briefly, it creates a look up table (LUT) and compares two types of search within LUT -
based on root mean squared error (RMSE) and gaussian process regression (GPR).


.. module:: src.+lut

LUT input sampling within parameter borders adjusted in ``src/+lut/input_borders.csv``

.. autofunction:: generate_lut_input

SCOPE run in time series mode with lut_in.csv data ...

Training and testing of gaussian process regression (GPR)

.. autofunction:: train_gpr

Preparation of data: flattening of tiff files to csv.

.. autofunction:: image2csv

Preparation of validation data (random pixels from flattened image)

.. autofunction:: pick100

SCOPE run in time series mode with validation.csv data ...

Validation of methods

.. autofunction:: validate_gpr

.. autofunction:: validate_rmse

Usage of methods

.. autofunction:: use_gpr

.. autofunction:: use_rmse


Supporting functions


.. autofunction:: change_detection

.. autofunction:: compress_geotiff

.. autofunction:: csv2image_plane

.. autofunction:: lut_search

.. autofunction:: plot_1to1

.. autofunction:: plot_image

.. autofunction:: write_tiff

