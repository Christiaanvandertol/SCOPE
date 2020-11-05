Multispectral retrieval (satellite)
====================================

Minimal input
---------------

Input is provided in ``./src/Input_data.xlsx`` that has several sheets.

To run the model for multispectral instruments use ``main_sat.m``.

.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Input``
      - Parameters and ranges
      - | Set 1 in the column **tune** in front of the parameter you want to retrieve
        | Change **value** of fixed parameters (where tune == 0)

    * - ``Satellite``
      - Multispectral retrieval
      - | Provide **image_path** (NetCDF4)
        | Select **instrument_name**
        | Provide band names used in your NetCDF for angles and coordinates
        | Go to ``Bands`` sheet

    * - ``Bands``
      - Band names and wl
      - | Provide **your_names** of bands as they are in NetCDF
        | Provide **your_wl** of bands (for plotting)
        | Keep the order! => leave a gap in bands if you excluded any


Retrieval on image subset (K x K)
-----------------------------------

If your image is big but you are interested only in a small area around a certain point you can specify
the coordinates of that point and the number of pixels around it.

.. list-table::

    * - sheet
      - purpose
      - action
    * - ``Satellite``
      - Multispectral retrieval
      - | Provide **pix_lat, pix_lon**
        | Provide **K** (number of pixels around that point)


Parallel computing (parfor)
-----------------------------

Fully supported as for hyperspectral retrieval :ref:`retrieval/hyperspectral:Parallel computing (parfor)`.

The changes should be made in ``main_sat.m`` (currently lines 105-111)


Output
-------

Windows
'''''''''

If you run ``main_sat.m`` you provide input file (**image_path**) as NetCDF4, the output will have two files:

1. ``%Y-%m-%d_%H%M%S.xlsx`` (see hyperspectral output :ref:`retrieval/hyperspectral:Windows`)
2. NetCDF4: retrieved parameter (Cab, LAI etc.) are be written as separate bands

So you can have a map of your retrieval.

Linux
''''''''

Nobody forbids you to do the same trick that is used for hyperspectral retrieval, just copy these lines from ``main.m`` to ``main_sat.m``

.. code-block:: matlab
    :caption: main.m
    :lineno-start: 122

    %% start saving
    q = parallel.pool.DataQueue;
    if isunix
        path = io.create_output_folder(input_path, path, tab.variable);

        tmp_zeros_res.rmse = rmse_all;
        tmp_zeros_res.parameters = parameters;
        tmp_zeros_unc.std_params = parameters_std;

        tmp_zeros_res.refl_mod = refl_mod;
        tmp_zeros_res.soil_mod = refl_soil;
        tmp_zeros_res.sif = sif_rad;
        tmp_zeros_res.sif_norm = sif_norm;

        tmp_zeros_meas.refl = refl_mod;
        tmp_zeros_meas.wl = measured.wl;
        tmp_zeros_meas.i_sif = measured.i_sif;

        io.save_output_csv(0, tmp_zeros_res, tmp_zeros_unc, tmp_zeros_meas, path)
        afterEach(q, @(x) io.save_output_csv(x{1}, x{2}, x{3}, x{4}, path));
    else
        path = io.create_output_file(input_path, path, measured, tab.variable);
        afterEach(q, @(x) io.save_output_j(x{1}, x{2}, x{3}, x{4}, path));
    end


Multispectral sensor from .txt file
--------------------------------------

If you do not want to retrieve from NetCDF4 file you can provide reflectance as .txt file and a proper sensor name.

Input is provided in ``./src/Input_data.xlsx``.

To run the model for multispectral instruments from .txt use ``main.m``.

.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Input``
      - Parameters and ranges
      - | Set 1 in the column **tune** in front of the parameter you want to retrieve
        | Change **value** of fixed parameters (where tune == 0)

    * - ``Filenames``
      - Multispectral retrieval
      - | Provide path to measured **reflectance** (text)
        | Provide path to wavelength of measurements **reflectance_wl**
        | Select **instrument_name**
        | Chose the number of the spectrum you fit **c** or all (-999)
        | Provide geometry **tts, tto, psi** or **lat, lon, datetime, tz**
        | Go to ``Bands`` sheet

    * - ``Bands``
      - Band names and wl
      - | Provide **your_names** of bands
        | Provide **your_wl** of bands (same as  **reflectance_wl**)
        | Keep the order! => leave a gap in bands if you excluded any
