V
=====

Variable names and values.

.. Note::
    - This is an array of 63 structs
    - This is a transitional data structure between ``input_data.xlsx`` file and :ref:`structs/input/input_index:input structs`

Initialized
""""""""""""

:func:`.assignvarnames` -> ``Name``

``SCOPE.m`` or :func:`.load_timeseries` -> ``Val``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - all
      - :func:`.select_input`

Fields
"""""""

Fields initialized in ``SCOPE.m``. Each of 63 structs in this array has these fields.

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **Name**
      - \-
      - char
      - defined in ``SCOPE.m``
      - SCOPE file identifiers
    * - **Val**
      - corresponding
      - char
      - \-
      - values from ``input_data.xlsx`` or files in timeseries


