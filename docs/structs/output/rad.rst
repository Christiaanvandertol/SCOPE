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

- :ref:`outfiles/each:radiation.csv`
- :ref:`outfiles/each:Eout_spectrum.csv`
- :ref:`outfiles/each:Lo_spectrum.csv`
- :ref:`outfiles/each:Esun.csv`
- :ref:`outfiles/each:Esky.csv`
- :ref:`outfiles/each:reflectance.csv`
- :ref:`outfiles/each:rsd.csv`
- :ref:`outfiles/each:rdd.csv`
- :ref:`outfiles/each:rso.csv`
- :ref:`outfiles/each:rdo.csv`

if ``options.calc_ebal``

- :ref:`outfiles/planck:spectrum_obsdir_BlackBody.dat`

if ``options.calc_planck``

Changes will be seen in

- :ref:`outfiles/each:Eout_spectrum.csv`
- :ref:`outfiles/each:Lo_spectrum.csv`

if ``options.calc_fluor``

- :ref:`outfiles/fluorescence:fluorescence_scalars.csv`
- :ref:`outfiles/fluorescence:fluorescence.csv`
- :ref:`outfiles/fluorescence:sigmaF.csv`
- :ref:`outfiles/fluorescence:fluorescence_hemis.csv`
- :ref:`outfiles/fluorescence:Lo_spectrum_inclF.csv`

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
        | ``Pnh_Cab, Pnu_Cab`` -> :ref:`structs/internal/biochem_in:Biochem_in`
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
      - directional back scattering coefficient for diffuse incidence
    * - **vf**
      - \-
      - [2162 x 1] double
      - directional forward scattering coefficient for diffuse incidence
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

.. Note:: Model simulated fluorescence at 3 levels:

    - level of photosystems individually (PSI, PSII) or together
    - level of leaves
    - level of canopy
        - in observation direction (reaching sensor) (typically starts with **Lo**)
        - hemispherically integrated

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **Fem_**
      - W m-2 um-1
      - [211 x 1] double
      - total emitted fluorescence by all leaves, excluding within canopy scattering / re-absorption
    * - **Fhem_**
      - W m-2 um-1
      - [211 x 1] double
      - TOC hemispherically integrated fluorescence
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
      - W m-2 um-1
      - [211 x 61] double
      - downward fluorescence flux profile
    * - **Fplu_**
      - W m-2 um-1
      - [211 x 61] double
      - upward fluorescence flux profile
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
      - TOC fluorescence contribution after scattering from leaves
    * - **LoF_soil**
      - W m-2 um-1 sr-1
      - [211 x 2] double
      - TOC fluorescence contribution after scattering from soil
    * - **Eoutf**
      - W m-2
      - double
      - hemispherically and spectrally integrated TOC fluorescence
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
      - W m-2 um-1
      - [211 x 1] double
      - total emitted fluorescence by all photosystems per wavelengths (excluding leaf and canopy re-absorption and scattering)
