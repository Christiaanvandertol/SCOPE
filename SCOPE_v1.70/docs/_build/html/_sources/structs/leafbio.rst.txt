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

Used
"""""

Most of the values are used by :func:`.fluspect_B_CX` (if ``options.calc_PSI``) OR :func:`.fluspect_B_CX_PSI_PSII_combined`,

By :func:`.ebal` to initialize :ref:`Biochem_in` structure

:func:`.load_timeseries`

==================================================================

Fields
"""""""

Fields initialized in ``SCOPE.m``

:V2Z: violaxantine to zeaxantine ratio. 0 if ``options.calc_PSI``

    :units: \-
    :type: double
    :default: 1

Fields initialized in :func:`.select_input` (read from ``input_data.xlsx``)

:Cab: Chlorophyll AB content

    :units: ug cm-2
    :type: double
    :default: 80.0

.. _leafbio.Cca:

:Cca: Carotenoid content. Usually 25% of Cab if ``options.Cca_function_of_Cab``

    :units: ug cm-2
    :type: double
    :default: 20.0

:Cdm: Dry matter content

    :units: g cm-2
    :type: double
    :default: 0.012

:Cw: leaf water equivalent layer

    :units: cm
    :type: double
    :default: 0.009

:Cs: scenecent material fraction

    :units: \-
    :type: double
    :default: 0.0

:Cant: Anthocyanins

    :units: ug cm-2
    :type: double
    :default: 0.0

:N: leaf thickness parameters

    :units: \-
    :type: double
    :default: 1.4

:Vcmo: maximum carboxylation capacity (at optimum temperature)

    :units: umol m-2 s-1
    :type: double
    :default: 60.0

:m: Ball-Berry stomatal conductance parameter

    :units: ?
    :type: double
    :default: 8.0

:Type: Photochemical pathway: 0 => 'C3', 1 => 'C4'

    :units: \-
    :type: int => char
    :default: 0 ('C3')

:Tparam: See ``PFT.xls``. These are five parameters specifying the temperature response.

    :units: ºK
    :type: [5 x 1] double
    :default: [0.2, 0.3, 281, 308, 328]

.. _leafbio.fqe:

:fqe: fluorescence quantum yield efficiency at photosystem level

    :units: \-
    :type: double | [2 x 1] double if ``options.calc_PSI``
    :default: 0.01

:Rdparam: ``Respiration = Rdparam * Vcmcax``

    :units: \-
    :type: double
    :default: 0.015

:rho_thermal: broadband thermal reflectance

    :units: \-
    :type: double
    :default: 0.01

:tau_thermal: broadband thermal transmittance

    :units: \-
    :type: double
    :default: 0.01

:Tyear: mean annual temperature

    :units: ºC
    :type: double
    :default: 15.0

:beta: fraction of photons partitioned to PSII (**0.507** for C3, **0.4** for C4; Yin et al. 2006 :cite:`YIN2006`; Yin and Struik 2012 :cite:`Yin2012`)

    :units: \-
    :type: double
    :default: 0.507

:kNPQs: rate constant of sustained thermal dissipation (Porcar-Castell 2011 :cite:`Porcar-Castell2011`)

    :units: s-1
    :type: double
    :default: 0.0

:qLs: fraction of functional reaction centres (Porcar-Castell 2011 :cite:`Porcar-Castell2011`)

    :units: \-
    :type: double
    :default: 1.0

:stressfactor: optional input: stress factor to reduce ``Vcmax`` (for example soil moisture, leaf age).

    :units: \-
    :type: double
    :default: 1.0
