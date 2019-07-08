Retrieval (Model inversion)
=============================

.. warning::
    Github of the project: https://github.com/Prikaziuk/demostrator

.. contents::

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
    * Multispectral instruments (Satellites):
        * data format - NetCDF4
        * parallel computing (parfor)
        * retrieval on an image subset (K x K pixels)

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
        - spectra is calculated on the spot (for each measurement during the run)
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


Input (Input_data.xlsx)
------------------------

Input is provided in ``./src/Input_data.xlsx`` that has several sheets.


Minimal hyperspectral
'''''''''''''''''''''''

To run the model for hyperspectral instruments use ``main.m``.

.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Input``
      - Parameters and ranges
      - | Set 1 in the column **tune** in front of the parameter you want to retrieve
        | Change **value** of fixed parameters (where tune == 0)

    * - ``Filenames``
      - Hyperspectral retrieval
      - | Provide path to measured **reflectance** (text)
        | Provide path to wavelength of measurements **reflectance_wl**
        | Select **instrument_name** or provide **FWHM**
        | Chose the number of the spectrum you fit **c** or all (-999)
        | Provide geometry **tts, tto, psi** or **lat, lon, datetime, tz, summertime**


Minimal multispectral (Satellite)
''''''''''''''''''''''''''''''''''''

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

Parallel computing
--------------------

parfor
'''''''

Each spectra is optimized on a single core (CPU). It is possible to use more cores (3 on modern computers) to speed up the processing.

1. Find and uncomment the following lines in ``main.m`` (currently 151-157).
2. change **for** to **parfor** in ``main.m``(currently 169)

.. code-block:: matlab

    %% uncomment these lines, select N_proc you want, change for-loop to parfor-loop
    N_proc = 3;
    if isempty(gcp('nocreate'))
    %     prof = parallel.importProfile('local_Copy.settings');
    %     parallel.defaultClusterProfile(prof);
        parpool(N_proc, 'IdleTimeout', Inf);
    end

    ...

    %% fitting
    %% change to parfor if you like
    parfor j = c
        ...
    end

.. note::
    Although parfor loops are, of course, faster, writing data to file from each iteration is slower (at least in the current implementation).
    We suggest first running one spectra without parallel computing to make sure you would not fail.
    Then write results to file **after** the parfor loop with :func:`io.save_output`

Speed up optimization (not recommended)
''''''''''''''''''''''''''''''''''''''''

If you prefer quantity over quality you may provide additional :func:`lsqnonlin` parameters in the :func:`fit_spectra`.

With **stoptol** == 1E-6 one spectra takes ~12 seconds, with **stoptol** == 1E-6 ~3 seconds (currently line 115).

.. code-block:: matlab

        stoptol = 1E-3;  % we recommend e-6

.. Warning::
    With stoptol 1E-3 we were not able to reproduce with high quality even SCOPE own spectra without any additional noise.


Time series
-------------

.. Note::
    The number of columns in any file from ``TimeSeres`` sheet has to be equal to those in the **reflectance** file.


Different angles for different spectra
''''''''''''''''''''''''''''''''''''''''

Usually you have more than one spectra to fit and those spectra were probably recorded at different time and with different angles.
Sensitivity analysis shows that solar and observation angles are crucial for accurate reflectance simulation.


.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Filenames``
      - Hyperspectral retrieval
      - | Put 1 in **timeseries** cell (B23)

    * - ``TimeSeres``
      - Paths (hyperspectral only)
      - | Enable **timeseries** on ``Filenames`` sheet
        | Provide paths to files with angles **tts_path**, **tto_path**, **psi_path**

You can provide only one path (for instance tts_path), then values for tto and psi are taken from the ``Filenames`` sheet.

Solar zenith angle (tts) from date, time and location
'''''''''''''''''''''''''''''''''''''''''''''''''''''''

With handheld spectrometers such as ASD observation zenith angle (tto) is 0 (nadir), which makes relative azimuth angle (phi) not important,
however solar zenith angle (tts) is constantly changing.

Good news! We can calculate solar zenith angle from date, time and coordinates of the measurement.

.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Filenames``
      - Hyperspectral retrieval
      - | Put 1 in **timeseries** cell (B23)
        | ``Delete`` value for **tts** (leave empty cell)
        | Provide **lat, lon** of the measurements
        | Provide **tz** timezone of measurements (UTC+tz)

    * - ``TimeSeres``
      - Paths (hyperspectral only)
      - | Enable **timeseries** on ``Filenames`` sheet
        | Provide **datetime_path** to the file with date and time of the measurements
        | Datetime format is ``%Y-%m-%d %H:%M:%S`` ("2019-07-01 12:30:20")

Timezone is provided in relation to UTC, so Netherlands are tz == 1 in Winter and tz == 2 in Summer.
If your time is already UTC tz == 0.

**summertime** option simply increments **tz**. It is the same providing for the Netherlands tz = 1, summertime = 1 or tz = 2, summertime = 0.

Validation plots
------------------

To see how well the retrieval worked for your spectra you can provide **validation** on ``Filenames`` tab.

Validation file requirements:
    - first columns - parameter names as they are on ``Input`` tab
    - other columns - measured values of this parameter for each spectra.

Overall the number of columns in **validation** file equals to the number of columns in the **reflectance** file **+1** for parameter names.

Example:

.. list-table::

    * - names
      - spec_1
      - spec_2
      - spec_3

    * - LAI
      - 1.0
      - NA
      - 1.5

    * - Cab
      - 40
      - 30
      - 50

.. figure:: ../images/quality_of_synthetic.png

Prior value selection
'''''''''''''''''''''''

.. Warning::
    Advanced users only

If you know your research area well, you might want to "recommend" the retrieval algorithm to stay close to, say, mean values of the parameters.

For this purpose on the ``Input`` tab:

1. Set **value** column to the mean value of your parameter.
2. Set **uncertainty** as the standard deviation of your parameter.
3. Uncomment the line in :func:`COST4SAIL` (currently 88)

.. code-block:: matlab

    er2 = 0;
    er2 = (p - prior.Apm) ./ prior.Aps;

    %% total error
    er = [er1 ; 3E-2* er2];  % change value of 3E-2 to higher / lower

Output
----------

We were experimenting with various output formats to satisfy Linux user and comply with the requirements of parfor loop.

Windows
'''''''''

``Input_data.xlsx`` is copied into **output_path** directory renamed as "%Y-%m-%d_%H%M%S.xlsx" ("2019-06-09-181952.xlsx") and the following sheets are written:


.. list-table::

    * - sheet
      - output
      - workspace matrix [#1]_

    * - ``Output``
      - | RMSE of spectral fit
        | retrieved parameter values
        | propagated standard deviation from **reflectance_std**
      - | `rmse_all`
        | `parameters` [#2]_
        | `parameters_std`

    * - ``Rmeas``
      - | wavelength of measurements
        | measured reflectance from **reflectance** file
      - | `measured.wl`
        | `measured.refl`

    * - ``Rmod``
      - | wavelength of measurements
        | simulated (best-fit) reflectance
      - `refl_mod`

    * - ``Rsoilmod``
      - | wavelength of measurements
        | simulated soil reflectance (BSM model)
      - `refl_soil`

    * - ``Fluorescence`` [#3]_
      - | wavelength of fluorescence 640:850
        | **only if SIF_PC were tuned**
        | fluorescence in radiance units (W m-2 sr-1)
      - `sif_rad`

    * - ``Fluorescence_norm`` [#3]_
      - | wavelength of fluorescence 640:850
        | **only if SIF_PC were tuned**
        | fluorescence in reflectance units
      - `sif_norm`

Those sheets are already present in ``Input_data.xlsx`` but are written later.

.. [#1] Each matrix (besides measured.refl) is preallocated with zeros and each column corresponds to the column in measured.refl.

    In this way if you tune only, say, the spectrum number 5 (c == 5) and you have 10 spectra in your **reflectance** file
    all these matrices will have 10 columns, 9 filled with zeros and the column number 5 with your retrieved values.


.. [#2] The order of row of *parameters* corresponds to the row of ``tab`` table read from ``Input`` sheet of ``Input_data.xlsx``.

    In other words row names of *parameters* == *tab.variable*

.. [#3] Currently sun-induced fluorescence (SIF) is reconstructed as a liner combination of the four principal components (SIF_PC1-4) to speed-up the retrieval.

    Although it can improve the fit in red-NIR region do not trust the values too much.

.. Note::
    You can load all the results back to matlab from the output "%Y-%m-%d_%H%M%S.xlsx" file with :func:`io.read_output()`


Linux
'''''''

Matlab can read .xlsx files but can't write into this format on Linux.
We hope you can configure ``Input_data.xlsx`` at your Linux machine or have it configured elsewhere.

On Linux inside **output_path** directory one more directory is created named as "%Y-%m-%d_%H%M%S" ("2019-06-09-181952").

``Input_data.xlsx`` is copied into that subfolder. All sheets listed in Windows section with the same information are written as separate .csv files.

Satellite
''''''''''

If you run ``main_sat.m`` you provide input file (**image_path**) as NetCDF4, the output will have two files:

1. "%Y-%m-%d_%H%M%S.xlsx" (see Windows section)
2. NetCDF4: retrieved parameter (Cab, LAI etc.) are be written as separate bands

So you can have a map of your retrieval.

Custom sensors
----------------

Sensor collects radiance with a certain spectral sampling interval (SSI) and defined sensor response function (SRF).
In the simples case SRF is just a gaussian curve with the centre at the band wavelength and width defined with full width half maximum (FWHM).

We convolve irradiance from **atmfile** (resolution 0.01 nm) (or **Esun, Esky** (resolution 1 nm)) and simulated reflected radiance to the sensor parameters.

For satellites SRFs are known and can be added in ``input/sensors.xlsx``.

For multispectral sensors we recommend reconstructing the SRF as a gaussian curve from FWHM with :func:`to_sensor.gaussian_fwhm()`

Hyperspectral instruments typically have individual sensors for visible and SWIR regions with different FWHM.
However, taken that SSI is < 3 nm and FWHM is relatively small, we do not think it is necessary to reconstruct the gaussian curve for each spectral region.

.. figure:: ../images/FWHM.png

    From this figure you can see that FWHM is important for accurate radiance simulation, but not so important for reflectance.

Plots
-------

We provide several plots for the validation of retrieval quality.

Reflectance
''''''''''''

Modelled and measured reflectance curves show if the modelled curve is close enough to the measured curve.

.. figure:: ../images/modelled_vs_measured.png


As there might be many spectra we do not show all the figures immediately after the run finishes
but we save them all in the matlab graphical array named ``figures`` with :func:`plot.reflectance_hidden()`
To actually see them set visibility property of the j-th figure to 'on':

.. code-block:: matlab

    set(figures(1), 'Visible', 'on')    % show the plot for the 1st spectrum
    set(figures(5:7), 'Visible', 'on')  % show the plots for the 5, 6, 7 spectra
    set(figures, 'Visible', 'on')       % show the plots for all spectra

If you close the figure (by pressing the upper right red cross) it will also disappear from ``figures`` array;
in fact there will be "handle to the deleted figure".

.. Note::
    You can draw all these plots again from the "%Y-%m-%d_%H%M%S.xlsx" file with :func:`plot.replot_all()`


Jacobian
''''''''''

The algorithm of local minimum search :func:`lsqnonlin` uses the gradient descent method
where the next step is driven by the cost-function partial derivative matrix - Jacobian matrix.

For each spectra we calculate the Jacobian matrix at the last step of optimization to propagate standard deviation of reflectance to that of the retrieved parameters.

Jacobian can be plot with :func:`plot.jacobian` and its SVD (singular value decomposition) can be plot with :func:`plot.jacobian_svd` for j-th spectra.

.. Note::
    We do not write the Jacobian matrix to file.
    If you want to analyze it plot it immediately after the model run or save as the .mat file yourself.

.. figure:: ../images/Jac.png

.. .. figure:: ../images/svd.png

Acknowledgement
-----------------

The project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 721995.

.. figure:: ../images/trustee_logo.png
    :align: center
