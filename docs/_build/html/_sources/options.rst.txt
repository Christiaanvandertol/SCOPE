Options
========

Effectively this are (almost [#]_ [#]_) all capabilities of the SCOPE model.

.. contents::
    :local:


This is an input structure that controls the workflow.

The values have binary (or tertiary) logic thus equal to 0 or 1 (or 2).

Influence on the output files is highlighted in the corresponding section :ref:`outfiles:Output files`

.. Note:: Not all combinations can bring to the desired result

Initialized
""""""""""""

``SCOPE.m``: read from ``input_data.xlsx`` or ``setoptions.m``


Rules of input reading
""""""""""""""""""""""""

``simulation``
-----------------------

Defines rules of input reading

Switch in ``SCOPE.m`` (multiple)

**0**
    **individual run(s)**: specify one value for fixed input parameters, and an equal number (> 1) of values for all parameters that vary between the runs.

**1**
    | **time series** (uses text files with meteo input as time series from *"../data/input/dataset X"* with files similar to ../data/input/:ref:`directories/data:dataset for_verification` specified on the ``filenames`` sheet of ``input_data.xslx``
    | :func:`.load_timeseries`

**2**
    **Lookup-Table**: specify a number of values in the row of input parameters. All possible combinations of inputs will be used.

Let us illustrate what the difference is in details.

It is possible to specify several values in a row on ``inputdata`` sheet of ``input_data.xslx``. Suppose we have an the following combination of input parameters. Notice, we provide two values for Cab and Cca parameters.

.. figure:: ./images/simulation.bmp

If **individual run(s)** (``options.simulation == 0``) was chosen the given combination will end up in **two** simulations:

* Cab=80, Cca=20
* Cab=40, Cca=10


If  **Lookup-Table** (``options.simulation == 2``) was chosen the given combination will end up in **four** simulations:

* Cab=80, Cca=20
* Cab=80, Cca=10
* Cab=40, Cca=20
* Cab=40, Cca=10

-----------------------

Variations in input
"""""""""""""""""""""

``rt_thermal``
-----------------------

Leaf and soil emissivity in thermal range

Switch in ``SCOPE.m``

**0**

    provide emissivity values as input :ref:`structs/input/leafbio:leafbio` (rho_thermal, tau_thermal), :ref:`structs/input/soil:soil`.rs_thermal

**1**
    use values from fluspect and soil at 2400 nm for the TIR range


--------------------------------


``calc_zo``
-----------------------

roughness length for momentum of the canopy (zo) and displacement height (d)

Switch in :func:`.select_input` :func:`.load_timeseries`

**0**

     zo and d values provided in the inputdata :ref:`structs/input/canopy:canopy`

**1**
    calculate zo and d from the LAI, canopy height, CD1, CR, CSSOIL (recommended if LAI changes in time series) :func:`zo_and_d`


--------------------------------


``soilspectrum``
-----------------------

Calculate soil reflectance or use from a file in ../data/input/:ref:`directories/data:soil_spectrum`

Switch in ``SCOPE.m``

**0**

    | use soil spectrum from the file with :ref:`structs/input/soil:soil`.spectrum
    | default file is ``soilnew.txt``, can be changed on the ``filenames`` sheet ``soil_file`` cell
    | variable name is ``rsfile``

**1**
    simulate soil spectrum with the BSM model (:func:`BSM`) parameters are fixed in code



--------------------------------


``soil_heat_method``
-----------------------

Method of ground heat flux (G) calculation

Switch in ``SCOPE.m``, :func:`.select_input`, :func:`.ebal`

**0**

    | standard calculation of thermal inertia from soil characteristic
    | :func:`.Soil_Inertia0` in :func:`.select_input`

**1**
    | empirically calibrated formula from soil moisture content :func:`.Soil_Inertia1` in :func:`.select_input`

**2**
    | as constant fraction (0.35) of soil net radiation
    | :func:`.Soil_Inertia0` in :func:`.select_input`


--------------------------------


``calc_rss_rbs``
-----------------------

soil resistance for evaporation from the pore space (rss) and soil boundary layer resistance (rbs)


Switch in :func:`.select_input`

**0**

    use resistance rss and rbs as provided in inputdata :ref:`structs/input/soil:soil`

**1**
    calculate rss from soil moisture content and correct rbs for LAI :func:`.calc_rssrbs`


--------------------------------


Variations in output
"""""""""""""""""""""

:func:`.RTMo` (SAIL) is executed in any valid run. Other functions may be included with these options.

--------------------------------

``calc_ebal``
--------------

Switch in ``SCOPE.m``

**0**


    Only :func:`.RTMo` is run (with :func:`.RTMf` if ``options.calc_fluor``)

**1**

    Calculate the complete energy balance.

    .. Warning:: required for ``calc_planck``, ``calc_directional``, ``calc_xanthophyllabs``

--------------------------------


``calc_planck``
-----------------------

Calculate spectrum of thermal radiation with spectral emissivity instead of broadband

.. Warning:: only effective with ``calc_ebal == 1``

Switch in ``SCOPE.m``, :func:`.calc_brdf`

**0**

    :func:`.RTMt_sb` - broadband brightness temperature is calculated in accordance to Stefan-Boltzman’s equation.

**1**
    | :func:`.RTMt_planck` is launched in ``SCOPE.m`` and :func:`.calc_brdf` (if ``calc_directional``).
    | Calculation is done per each wavelength thus takes more time than Stefan-Boltzman.


--------------------------------


``calc_directional``
-----------------------

Calculate BRDF and directional temperature for many angles specified in the file: :ref:`directories/data:directional`.

.. Warning::
    - only effective with ``calc_ebal == 1``
    - Be patient, this takes some time

Switch in ``SCOPE.m``, :func:`.calc_brdf`

**0**

    -

**1**
    | struct :ref:`structs/output/directional:directional` is loaded from the file :ref:`directories/data:directional`
    | :func:`.calc_brdf` is launched in ``SCOPE.m``



--------------------------------


``calc_xanthophyllabs``
-----------------------

Calculate dynamic xanthopyll absorption (zeaxanthin) for simulating PRI (photochemical reflectance index)

.. Warning::
    - only effective with ``calc_ebal == 1``

Switch in ``SCOPE.m``

**0**

    -

**1**
    :func:`.RTMz` is launched in ``SCOPE.m`` and :func:`.calc_brdf` (if ``calc_directional``)


--------------------------------

``calc_vert_profiles``
-----------------------

Calculation of vertical profiles (per 60 canopy layers).

Corresponding structure :ref:`structs/output/profiles:profiles`

Switch in ``SCOPE.m``, :func:`.RTMo` and :func:`.ebal`

**0**

    Profiles are not calculated

**1**
    | Photosynthetically active radiation (PAR) per layer is calculated in :func:`.RTMo`
    | Energy, temperature and photosynthesis fluxes per layer are calculated in :func:`.ebal`
    | Fluorescence fluxes are calculated in :func:`.RTMf` if (``calc_fluor``)


--------------------------------

``calc_fluor``
-----------------------

Calculation of fluorescence

Switch in ``SCOPE.m``, :func:`.calc_brdf`

**0**

    No fluorescence

**1**
    | :func:`.RTMf` is launched in ``SCOPE.m`` and :func:`.calc_brdf` (if ``calc_directional``)
    | total emitted fluorescence is calculated by ``SCOPE.m``


--------------------------------

``calc_PSI``
-----------------------

Separate fluorescence of photosystems I and II (PSI, PSII) or not

Switch in ``SCOPE.m``, :func:`.select_input`

**0**

    | **recommended**
    | treat the whole fluorescence spectrum as one spectrum (new calibrated optipar)
    | fluspect version :func:`.fluspect_B_CX_PSI_PSII_combined`

**1**
    | differentiate PSI and PSII with Franck et al. spectra (of SCOPE 1.62 and older)
    | fluspect version :func:`.fluspect_B_CX`
    | fluorescence quantum efficiency of PSI is set to 0.2 of PSII in :func:`.select_input`


--------------------------------


``Fluorescence_model``
-----------------------

Fluorescence model

Switch in :func:`.ebal`

**0**

    empirical, with sustained NPQ (fit to Flexas' data)

**1**
    empirical, with sigmoid for Kn: :func:`.biochemical` (Berry-Van der Tol)

**2**
    :func:`.biochemical_MD12` (von Caemmerer-Magnani)


--------------------------------



``apply_T_corr``
-----------------------

correct Vcmax and rate constants for temperature

.. Warning::
    only effective with ``Fluorescence_model != 2`` i.e. for :func:`.biochemical`

Switch in :func:`.ebal`

**0**
    -

**1**
    correction in accordance to Q10 rule


--------------------------------

For users' comfort
"""""""""""""""""""""""""

``verify``
-----------------------

verify the results (compare to saved 'standard' output) to test the code for the first time

Switch in ``SCOPE.m``

**0**
    -

**1**
    runs :func:`.output_verification`


--------------------------------


``save_headers``
-----------------------

write header lines in output files

Switch in :func:`.create_output_files`

**0**
    -

**1**
    runs additional section in :func:`.create_output_files` which writes two lines (names, units) in output files


--------------------------------


``makeplots``
-----------------------

plot the results

Switch in ``SCOPE.m``

**0**
    -

**1**
    launches :func:`.plots` for the results of the last run

.. [#] extra output variables that are not saved to files (see :ref:`structs:Structs`) are available in the workspace after the model run.
.. [#] model can be varied by user, please, consult :ref:`api:API` to learn signatures of functions