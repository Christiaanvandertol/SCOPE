iter
=====
Numerical parameters, such as the number of iterations needed to reach energy balance closure

Initialized
""""""""""""

``SCOPE.m``

Variations
""""""""""""

``counter`` is incremented in :func:`.ebal`

``Wc`` is set to 0.2 if ``counter`` > 50 in :func:`.ebal`


Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - ``maxit, maxEBer, Wc``
      - :func:`.ebal`
    * - ``counter``
      - | :func:`.initialize_output_structures`
        | :func:`.output_data`

Fields
"""""""

Fields initialized in ``SCOPE.m``

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **maxit**
      - \-
      - int
      - 100
      - maximum number of iterations
    * - **maxEBer**
      - W m-2
      - double
      - 1.0
      - maximum accepted error in energy balance
    * - **Wc**
      - \-
      - double
      - 1.0
      - weight coefficient for iterative calculation of Tc
    * - **counter**
      - \-
      - int
      - 0
      - counter, changed in :func:`.ebal`
