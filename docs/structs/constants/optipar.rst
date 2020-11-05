optipar
========

Leaf optical parameters: specific absorption coefficients (SAC) of leaf chemical components

Concentration * SAC

Initialized
""""""""""""

``SCOPE.m``: read from :ref:`directories/input:fluspect_parameters`

Variations
""""""""""""

Different files (corresponding to PROSPECT versions) from :ref:`directories/input:fluspect_parameters` can be selected in ``filenames.csv`` (``leaf_files`` cell) or uncomment lines directly in ``SCOPE.m`` (~ 170)

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - | ``nr, Kdm, Kab, Kw, Ks, Kant``
        | ``Kca (or KcaV, KcaZ)``
      - | :func:`.fluspect_B_CX`
        | :func:`.fluspect_B_CX_PSI_PSII_combined`
    * - | ``phiI, phiII``  (if ``options.calc_PSI``)
        | ``phi``
      - | ``SCOPE.m`` :func:`.fluspect_B_CX` if ``options.calc_PSI``
        | ``SCOPE.m`` :func:`.fluspect_B_CX_PSI_PSII_combined`

Fields
"""""""

Fields loaded in ``SCOPE.m``

.. list-table::
    :widths: 10 20 20 50

    * - variable
      - units
      - type
      - description
    * - **wl**
      - nm
      - [2001 x 1] double
      - SAC wavelength range 400-2400
    * - **nr**
      - nm-1
      - [2001 x 1] double
      - refractive index
    * - **Kab**
      - nm-1
      - [2001 x 1] double
      - SAC of chlorophylls a + b
    * - **Kca**
      - nm-1
      - [2001 x 1] double
      - SAC of carotenoids (violaxanthin + zeaxanthin)
    * - **Ks**
      - nm-1
      - [2001 x 1] double
      - SAC of senescent material
    * - **Kw**
      - nm-1
      - [2001 x 1] double
      - SAC of water
    * - **Kdm**
      - nm-1
      - [2001 x 1] double
      - SAC of dry matter
    * - **phiI**
      - nm-1
      - [2001 x 1] double
      - SAC of PSI
    * - **phiII**
      - nm-1
      - [2001 x 1] double
      - SAC of PSI
    * - **phi**
      - nm-1
      - [2001 x 1] double
      - SAC of PSI + PSII
    * - **KcaV**
      - nm-1
      - [2001 x 1] double
      - SAC of violaxanthin
    * - **KcaZ**
      - nm-1
      - [2001 x 1] double
      - SAC of zeaxanthin
    * - **Kant**
      - nm-1
      - [2001 x 1] double
      - SAC of anthocyanins
    * - **KcaV2**
      - nm-1
      - [2001 x 1] double
      - former SAC of violaxanthin
    * - **GSV**
      - nm-1
      - [2001 x 3] double
      - Global Soil Vectors spectra
    * - **nw**
      - nm-1
      - [2001 x 1] double
      - water refraction index spectrum
