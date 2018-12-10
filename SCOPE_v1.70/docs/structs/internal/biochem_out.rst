Biochem_out
============

.. Note::
    This is an internal struct of :func:`.ebal`. It is only available within :func:`.ebal` in debug mode.

Output of the biochemical routine :func:`.biochemical` or :func:`.biochemical_MD12` for **leaf** photosynthesis and fluorescence.

Initialized
""""""""""""

:func:`.ebal`

Variations
""""""""""""

For sunlit leaves size is [13 x 36 x 60] for shaded [60 x 1]


Fields
"""""""

Output specific for :func:`.biochemical`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **A**
      - umol m-2 s-1
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - photosynthesis (CO2 assimilation rate) per canopy layer
    * - **Ci**
      - ppm
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - within leaf CO2 concentration per canopy layer
    * - **eta**
      - \-
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - Fs / Fo
    * - **rcw**
      - \-
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - stomatal resistance
    * - **qE**
      - \-
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - non-photochemical quenching 1 - (Fm - Fo) / (Fm0 -Fo0)
    * - **Kn**
      - \-
      - | [13 x 36 x 60] double
        | [60 x 1] double
      - NPQ = (Fm - Fm')/Fm' = Kn/(Kf + Kd)

Output specific specific for :func:`.biochemical_MD12`

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
