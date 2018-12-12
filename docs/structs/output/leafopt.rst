leafopt
========

Leaf optical properties

.. Note::
    - This output is used, available in the workspace after a SCOPE run but is **not** written to any output file
    - **leafoptZ** is exactly the same struct calculated by selected version of Fluspect but when all violaxanthin is transformed to zeaxanthin (``V2Z == 1``)


Calculated
""""""""""""

:func:`.fluspect_B_CX_PSI_PSII_combined`

Variations
""""""""""""

:func:`.fluspect_B_CX` if ``options.calc_PSI`` =>
    | ``Mb`` is partitioned to ``MbI, MbII``
    | ``Mf`` is partitioned to ``MfI, MfII``


Output file
""""""""""""

- not saved

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``refl, tran``
      - | :func:`.RTMf`
        | :func:`.RTMo`
        | :func:`.RTMt_planck`
        | :func:`.RTMt_sb`
        | :func:`.RTMz`
    * - ``Mb, Mf``
      - :func:`.RTMf`
    * - ``kChlrel``
      - :func:`.RTMo`
    * - ``reflZ, tranZ``
      - :func:`.RTMz`


Fields
"""""""

Fields calculated in :func:`.fluspect_B_CX_PSI_PSII_combined` or
:func:`.fluspect_B_CX` if ``options.calc_PSI``

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **refl**
      - \-
      - [2162 x 1] double
      - leaf reflectance
    * - **tran**
      - \-
      - [2162 x 1] double
      - leaf transmittance
    * - **kChlrel**
      - \-
      - [2001 x 1] double
      -
    * - **Mb**
      - \-
      - [211 x 351] double
      - backward scattering fluorescence matrix
    * - **Mf**
      - \-
      - [211 x 351] double
      - forward scattering fluorescence matrix
    * - **Mbl_rc**
      - \-
      - [211 x 351] double
      - backward scattering fluorescence matrix without re-absorption
    * - **Mfl_rc**
      - \-
      - [211 x 351] double
      - forward scattering fluorescence matrix without re-absorption


These values are added from **leafoptZ**, calculated with :ref:`structs/input/leafbio:leafbio`.V2Z == 1.

.. list-table::
    :widths: 10 10 20 60

    * - variable
      - units
      - type
      - description
    * - **reflZ**
      - \-
      - [2162 x 1] double
      - leaf reflectance with only zeaxanthin
    * - **tranZ**
      - \-
      - [2162 x 1] double
      - leaf transmittance with only zeaxanthin