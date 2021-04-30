leafbio
========
Leaf biochemical parameters

Initialized
""""""""""""
:func:`.select_input`

``SCOPE.m``

Variations
""""""""""""
leafbio.Cca_ may be calculated as 25% of Cab: ``options.Cca_function_of_Cab``

leafbio.fqe_ may be double (PSII only) or [2 x 1] double (PSI = 0.2 * PSII, PSII) if ``options.calc_PSI``

leafbio.V2Z_ can be set to 0 with ``options.calc_PSI``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``Cab, Cca, V2Z, Cw, Cdm, Cs, Cant, Cp, Cbc, N, fqe``
      - | :func:`.fluspect_B_CX` if ``options.calc_PSI``
        | :func:`.fluspect_B_CX_PSI_PSII_combined`
    * - ``Type, BallBerrySlope, Rdparam, Tyear, Kn0, Knalpha, Knbeta, beta, qLs, kNPQs, stressfactor, Tparam, Vcmax25`` -> :ref:`structs/internal/biochem_in:Biochem_in`
      - :func:`.ebal`
    * - ``Vcmax25, Cab``
      - :func:`.load_timeseries`
    * - ``rho_thermal, tau_thermal, fqe``
      - ``SCOPE.m``

Fields
"""""""

Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **Cab**
      - ug cm-2
      - double
      - 80.0
      - Chlorophyll AB content
    * - .. _leafbio.Cca:

        **Cca**
      - ug cm-2
      - double
      - 20.0
      - Carotenoid content. Usually 25% of Cab if ``options.Cca_function_of_Cab``
    * - **Cdm**
      - g cm-2
      - double
      - 0.012
      - Dry matter content
    * - **Cw**
      - cm
      - double
      - 0.009
      - leaf water equivalent layer
    * - **Cs**
      - \-
      - double
      - 0.0
      - senescent material fraction
    * - **Cant**
      - ug cm-2
      - double
      - 0.0
      - Anthocyanins
    * - **Cp**
      - ug cm-2
      - double
      - 0.0
      - protein
    * - **Cbc**
      - ug cm-2
      - double
      - 0.0
      - ?
    * - **N**
      - \-
      - double
      - 1.4
      - leaf thickness parameters
    * - **Vcmax25**
      - umol m-2 s-1
      - double
      - 60.0
      - maximum carboxylation capacity (at optimum temperature of 25C, former Vcmo)
    * - **BallBerrySlope**
      - -
      - double
      - 8.0
      - slope of Ball-Berry stomatal conductance model (former m)
    * - **BallBerry0**
      - -
      - double
      - 0.01
      - intercept of Ball-Berry stomatal conductance model
    * - **Kn0**
      - -
      - double
      - 2.48
      - Kn0: parameter for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha);
    * - **Knalpha**
      - -
      - double
      - 2.83
      - alpha parameter for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha)
    * - **Knbeta**
      - -
      - double
      - 0.114
      - beta parameter for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha)
    * - **Type**
      - \-
      - int => char
      - 0 ('C3')
      - Photochemical pathway: 0 => 'C3', 1 => 'C4'
    * - **Tparam**
      - ºK
      - [5 x 1] double
      - [0.2, 0.3, 281, 308, 328]
      - See ``PFT.xls``. These are five parameters specifying the temperature response.
    * - .. _leafbio.fqe:

        **fqe**
      - \-
      - [1 x 1] | [2 x 1] double
      - 0.01
      - fluorescence quantum yield efficiency at photosystem level
    * - **Rdparam**
      - \-
      - double
      - 0.015
      - ``Respiration = Rdparam * Vcmcax``
    * - **rho_thermal**
      - \-
      - double
      - 0.01
      - broadband thermal reflectance
    * - **tau_thermal**
      - \-
      - double
      - 0.01
      - broadband thermal transmittance
    * - **Tyear**
      - ºC
      - double
      - 15.0
      - mean annual temperature
    * - **beta**
      - \-
      - double
      - 0.507
      - fraction of photons partitioned to PSII (**0.507** for C3, **0.4** for C4; Yin et al. 2006 :cite:`YIN2006`; Yin and Struik 2012 :cite:`Yin2012`)
    * - **kNPQs**
      - s-1
      - double
      - 0.0
      - rate constant of sustained thermal dissipation (Porcar-Castell 2011 :cite:`Porcar-Castell2011`)
    * - **qLs**
      - \-
      - double
      - 1.0
      - fraction of functional reaction centres (Porcar-Castell 2011 :cite:`Porcar-Castell2011`)
    * - **stressfactor**
      - \-
      - double
      - 1.0
      - optional input: stress factor to reduce ``Vcmax`` (for example soil moisture, leaf age)


Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - .. _leafbio.V2Z:

        **V2Z**
      - \-
      - double
      - 1.0
      - violaxantine to zeaxantine ratio. 0 if ``options.calc_PSI``