Hyperspectral retrieval
=========================

Minimal input
---------------

Input is provided in ``./src/Input_data.xlsx`` that has several sheets.

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
        | Provide geometry **tts, tto, psi**


I do not know the observation geometry
----------------------------------------

With handheld spectrometers such as ASD observation zenith angle (tto) is 0 (nadir), which makes relative azimuth angle (phi) not important,
however solar zenith angle (tts) is constantly changing.

Good news! We can calculate solar zenith angle from date, time and coordinates of the measurement.

.. list-table::

    * - sheet
      - purpose
      - action

    * - ``Filenames``
      - Hyperspectral retrieval
      - | ``Delete`` value for **tts** (leave empty cell)
        | Provide **lat, lon** of the measurements
        | Provide **datetime**
        | Datetime format is ``%Y-%m-%d %H:%M:%S`` ("2019-07-01 12:30:20")
        | Provide **tz** timezone of measurements (UTC+tz)

Timezone is provided in relation to UTC, so Netherlands are tz == 1 in Winter and tz == 2 in Summer.
If your time is already UTC tz == 0.

**summertime** option simply increments **tz**. It is the same providing for the Netherlands tz = 1, summertime = 1 or tz = 2, summertime = 0.


Time series (different angles for different spectra)
------------------------------------------------------

.. Note::
    The number of columns in any file from ``TimeSeres`` sheet has to be equal to those in the **reflectance** file.


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
        | Or **datetime_path** to the file with date and time of the measurements
        | Datetime format is ``%Y-%m-%d %H:%M:%S`` ("2019-07-01 12:30:20")

You can provide only one path (for instance tts_path), then values for tto and psi are taken from the ``Filenames`` sheet.


Parallel computing (parfor)
----------------------------

Each spectra is optimized on a single core (CPU). It is possible to use more cores (3 on modern computers) to speed up the processing.

Obviously, it is useful only with Time series mode.

1. Find and uncomment the following lines in ``main.m`` (currently 149-154)

.. code-block:: matlab
    :lineno-start: 149
    :caption: main.m

    %% uncomment these lines, select N_proc you want, change for-loop to parfor-loop
    N_proc = 3;
    if isempty(gcp('nocreate'))
    %     prof = parallel.importProfile('local_Copy.settings');
    %     parallel.defaultClusterProfile(prof);
        parpool(N_proc, 'IdleTimeout', Inf);
    end

2. Change **for** to **parfor** in ``main.m`` (currently 166)

.. code-block:: matlab
    :lineno-start: 164
    :caption: main.m

    %% fitting
    %% change to parfor if you like
    parfor j = c
        ...
    end

.. note::
    Although parfor loops are, of course, faster, writing data to file from each iteration is slower (at least in the current implementation).
    We suggest first running one spectra without parallel computing to make sure you would not fail.
    Then write results to file **after** the parfor loop with :func:`io.save_output`


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
