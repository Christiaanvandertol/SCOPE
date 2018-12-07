Options
========

Simulation options, such as time series or look-up tables, fluorescence.

The values have binary (or tertiary) logic thus equal to 0 or 1 (or 2).

Influence on the output files is highlighted in the corresponding section :ref:`output_files:Output files`

.. Note:: Not all combinations are possible

Initialized
""""""""""""

``SCOPE.m``: read from ``input_data.xlsx`` or ``setoptions.m``


Effects
""""""""""

:func:`.RTMo` (SAIL) is executed in any valid run. Other functions may be executed with this options.


``calc_ebal``
--------------

Switch in ``SCOPE.m``

**0**


    Only :func:`.RTMo` is run (with :func:`.RTMf` if ``options.calc_fluor``)

**1**

    Calculate the complete energy balance.

    .. Warning:: required for ``calc_planck``, ``calc_directional``, ``calc_xanthophyllabs``

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


``calc_fluor``
-----------------------

Calculation of fluorescence

Switch in ``SCOPE.m``, :func:`.calc_brdf`

**0**

    No fluorescence

**1**
    | :func:`.RTMf` is launched in ``SCOPE.m`` and :func:`.calc_brdf` (if ``calc_directional``)
    | total emitted fluorescence is calculated by ``SCOPE.m``


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


``rt_thermal``
-----------------------

Leaf and soil emissivity in thermal range

Switch in ``SCOPE.m``

**0**

    provide emissivity values as input :ref:`structs/input/leafbio:leafbio` (rho_thermal, tau_thermal), :ref:`structs/input/soil:soil`.rs_thermal

**1**
    use values from fluspect and soil at 2400 nm for the TIR range


``calc_zo``
-----------------------

roughness length for momentum of the canopy (zo) and displacement height (d)

Switch in :func:`.select_input` :func:`.load_timeseries`

**0**

     zo and d values provided in the inputdata :ref:`structs/input/canopy:canopy`

**1**
    calculate zo and d from the LAI, canopy height, CD1, CR, CSSOIL (recommended if LAI changes in time series) :func:`zo_and_d`


``soilspectrum``
-----------------------

Calculate soil reflectance or use from a file in ../data/input/:ref:`directories/data:soil_spectrum`

Switch in ``SCOPE.m``

**0**

    | use soil spectrum from the file with :ref:`structs/input/soil:soil`.spectrum
    | default file is ``soilnew.txt``, can be changed on the ``filenames`` sheet ``soil_file`` cell

**1**
    simulate soil spectrum with the BSM model (:func:`BSM`) parameters are fixed in code



``soil_heat_method``
-----------------------

Method of ground heat flux (G) calculation

Switch in ``SCOPE.m``, :func:`.select_input`, :func:`.ebal`

.. Error:: In ebal it is either 1 or 2, 0 is redirected to 2. Depends of ``options.simulation``

**0**

    | standard calculation of thermal inertia from soil characteristic
    | :func:`.Soil_Inertia0` in :func:`.select_input`

**1**
    | empirically calibrated formula from soil moisture content :func:`.Soil_Inertia1` in :func:`.select_input`

**2**
    | as constant fraction of soil net radiation
    | :func:`.Soil_Inertia0` in :func:`.select_input`



``Fluorescence_model``
-----------------------

Fluorescence model

.. Error:: 0 == 2. Maybe, not fluorescence but biochemical model?

Switch in :func:`.ebal`

**0**

    empirical, with sustained NPQ (fit to Flexas' data)

**1**
    empirical, with sigmoid for Kn: :func:`.biochemical` (Berry-Van der Tol)

**2**
    :func:`.biochemical_MD12` (von Caemmerer-Magnani)

