soilemp
========

.. Note::
    This is an optional struct created for :func:`.BSM` with ``options.soilspectrum == 1``


Initialized
""""""""""""

``SCOPE.m``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``SMC, film``
      - :func:`.BSM`

Fields
"""""""

Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 10 20 50

    * - variable
      - units
      - type
      - default
      - description
    * - **SMC**
      -
      - double
      - 25
      - soil moisture capacity parameter
    * - **film**
      -
      - double
      - 0.015
      - single water film optical thickness
