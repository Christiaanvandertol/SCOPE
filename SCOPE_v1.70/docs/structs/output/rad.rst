rad
====

Radiation fluxes: both input (MODTRAN) and output

A large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.

Fields are initialized by :func:`.initialize_output_structures`

Calculated
""""""""""""
:func:`.RTMo`

:func:`.RTMf`

:func:`.RTMt_planck`

:func:`.RTMt_sb`

:func:`.ebal`

``SCOPE.m``

Output file
""""""""""""

- :ref:`outfiles/radiation:radiation.dat`
- :ref:`outfiles/spectrum_obsdir_optical:spectrum_obsdir_optical.dat`
- :ref:`outfiles/spectrum_hemis_optical:spectrum_hemis_optical.dat`
- :ref:`outfiles/irradiance_spectra:irradiance_spectra.dat`
- :ref:`outfiles/reflectance:reflectance.dat`
- :ref:`outfiles/BOC_irradiance:BOC_irradiance.dat`

if ``options.calc_ebal``

- :ref:`outfiles/spectrum_obsdir_BlackBody:spectrum_obsdir_BlackBody.dat`

if ``options.calc_planck``

- :ref:`outfiles/spectrum_hemis_thermal:spectrum_hemis_thermal.dat`
- :ref:`outfiles/spectrum_obsdir_thermal:spectrum_obsdir_thermal.dat`

if ``options.calc_fluor``

- :ref:`outfiles/fluorescence:fluorescence.dat`
- :ref:`outfiles/fluorescencePSI:fluorescencePSI.dat`
- :ref:`outfiles/fluorescencePSII:fluorescencePSII.dat`
- :ref:`outfiles/fluorescence_hemis:fluorescence_hemis.dat`
- :ref:`outfiles/fluorescence_emitted_by_all_leaves:fluorescence_emitted_by_all_leaves.dat`
- :ref:`outfiles/fluorescence_emitted_by_all_photosystems:fluorescence_emitted_by_all_photosystems.dat`
- :ref:`outfiles/fluorescence_sunlit:fluorescence_sunlit.dat`
- :ref:`outfiles/fluorescence_shaded:fluorescence_shaded.dat`
- :ref:`outfiles/fluorescence_scattered:fluorescence_scattered.dat`

Variations
""""""""""""

if ``options.calc_PSI`` fluorescence (``LoF_``) is partitioned between photosystems ``LoF1_, LoF2_``


Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``Lot LoF_``
      - :func:`.calc_brdf`
    * - | ``Rnuc, Rnhct, Rnuct, Rnhst, Rnust, Rnhc, Rnuc, Rnhs, Rnus``
        | ``Pnh_Cab, Pnu_Cab`` -> :ref:`structs/biochem_in:Biochem_in`
        | ``Pnh, Pnu, Pnh_PAR, Pnu_PAR``
        | ``Eoutte``
      - :func:`.ebal`
    * - ``vb, vf, Esun_, Emin_, Eplu``
      - | :func:`.RTMf`
        | :func:`.RTMz`
    * - ``Pnh, Pnu, Pnh_Cab, Pnu_Cab, Rnh_PAR, Rnu_PAR``
      - ``SCOPE.m``


Fields
"""""""

Fields initialized in :func:`.RTMo`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **rsd**
      - \-
      - [2162 x 1] double
      - conical-hemispherical reflectance factor (specular in -> diffuse out)
    * - **rdd**
      - \-
      - [2162 x 1] double
      - bihemispherical reflectance factor (diffuse in -> diffuse out)
    * - **rdo**
      - \-
      - [2162 x 1] double
      - hemispherical-conical reflectance factor (diffuse in -> specular out)
    * - **rso**
      - \-
      - [2162 x 1] double
      - biconical reflectance factor (specular in -> specular out)
    * - **vb**
      - \-
      - [2162 x 1] double
      -
    * - **vf**
      - \-
      - [2162 x 1] double
      -
    * - **Esun_**
      - mW m-2 um-1
      - [2162 x 1] double
      - incident solar spectrum
    * - **Esky_**
      - mW m-2 um-1
      - [2162 x 1] double
      - incident sky spectrum
    * - **PAR**
      - mol m-2 s-1
      - double
      - incident spectrally integrated PAR
    * - **fEsuno**
      - \-
      - [2162 x 1] double
      - fraction of direct light (optical)
    * - **fEskyo**
      - \-
      - [2162 x 1] double
      - fraction of diffuse light (optical)
    * - **fEsunt**
      - \-
      - [2162 x 1] double
      - fraction of direct light (thermal)
    * - **fEskyt**
      - \-
      - [2162 x 1] double
      - fraction of diffuse light (thermal)
    * - **Eplu_**
      - mW m-2 um-1
      - [61 x 2162] double
      - upward diffuse radiation in the canopy
    * - **Emin_**
      - mW m-2 um-1
      - [61 x 2162] double
      - downward diffuse radiation in the canopy
    * - **Lo_**
      - mW m-2 um-1 sr-1
      - [2162 x 1] double
      - top of canopy (TOC) radiance in observation direction
    * - **Eout_**
      - mW m-2 um-1
      - [2162 x 1] double
      - top of canopy (TOC) upward radiation
    * - **Eouto**
      - W m-2
      - double
      - spectrally integrated upward optical radiation
    * - **Eoutt**
      - W m-2
      - double
      - spectrally integrated upward thermal radiation
    * - **Rnhs**
      - W m-2
      - double
      - net radiation of shaded soil
    * - **Rnus**
      - W m-2
      - double
      - net radiation of sunlit soil
    * - **Rnhc**
      - W m-2
      - [60 x 1] double
      - net radiation of shaded leaves
    * - **Rnuc**
      - W m-2
      - [13 x 36x 60] double
      - net radiation of sunlit leaves
    * - **Pnh**
      - mol n-2 s-1
      - [60 x 1] double
      - net PAR of shaded leaves
    * - **Pnu**
      - mol n-2 s-1
      - [13 x 36x 60] double
      - net PAR of sunlit leaves
    * - **Pnh_Cab**
      - mol n-2 s-1
      - [60 x 1] double
      - net PAR absorbed by Cab of shaded leaves
    * - **Pnu_Cab**
      - mol n-2 s-1
      - [13 x 36x 60] double
      - net PAR absorbed by Cab of sunlit leaves
    * - **Pnh_PAR**
      - W m-2
      - [60 x 1] double
      - net PAR of shaded leaves (W m-2)
    * - **Pnu_PAR**
      - W m-2
      - [13 x 36x 60] double
      - net PAR of sunlit leaves (W m-2)
    * - **Etoto**
      -
      - double
      -

Fields initialized in :func:`.RTMf`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Fem_**
      - W m-2 um-1
      - [211 x 1] double
      -
    * - **Fhem_**
      - W m-2 um-1
      - [211 x 1] double
      -
    * - **LoF_**
      - W m-2 um-1 sr-1
      - [211 x 1] double
      - fluorescence per wavelength
    * - **LoF1_**
      - W m-2 um-1 sr-1
      - [211 x 1] double
      - fluorescence from photosystem I (PSI) per wavelength
    * - **LoF2_**
      - W m-2 um-1 sr-1
      - [211 x 1] double
      - fluorescence from photosystem II (PSII) per wavelength
    * - **Fhem_**
      - W m-2 um-1
      - [211 x 1] double
      -
    * - **Fmin_**
      - W m-2 sr-1
      - [211 x 61] double
      -
    * - **Fplu_**
      - W m-2 sr-1
      - [211 x 61] double
      -
    * - **LoF_sunlit**
      - W m-2 um-1 sr-1
      - [211 x 2] double
      - TOC fluorescence contribution from sunlit leaves in observer direction per wavelengths
    * - **LoF_shaded**
      - W m-2 um-1 sr-1
      - [211 x 2] double
      - TOC fluorescence contribution from shaded leaves in observer direction per wavelengths
    * - **LoF_scattered**
      - W m-2 um-1 sr-1
      - [211 x 2] double
      - TOC fluorescence contribution from leaves after scattering
    * - **LoF_soil**
      - W m-2 um-1 sr-1
      - [211 x 2] double
      - TOC fluorescence contribution from soil after scattering
    * - **Eoutf**
      - W m-2 sr-1
      - double
      -
    * - **Eminf_**
      - W m-2 sr-1
      - [61 x 21] double
      -
    * - **Epluf_**
      - W m-2 sr-1
      - [61 x 21] double
      -

Fields initialized in :func:`.RTMt_planck`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Lot_**
      -
      - double
      -
    * - **Eoutte_**
      -
      - double
      -
    * - **Eplut_**
      -
      - [61 x 1] double
      -
    * - **Emint_**
      -
      - [61 x 1] double
      -

Fields initialized in :func:`.RTMt_sb`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Lot**
      -
      - double
      -
    * - **Eoutte**
      -
      - double
      -
    * - **Eplut**
      -
      - [61 x 1] double
      -
    * - **Emint**
      -
      - [61 x 1] double
      -
    * - **Rnuct**
      -
      - [13 x 36 x 60] double
      -
    * - **Rnhct**
      -
      - [60 x 1] double
      -
    * - **Rnust**
      -
      - double
      -
    * - **Rnhst**
      -
      - double
      -

Fields added in :func:`.ebal`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **LotBB_**
      - W m-2 sr-1
      - [2162 x 1] double
      - blackbody radiance

Fields added in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Femtot**
      -
      - [211 x 1] double
      - total emitted fluorescence by all photosystems
