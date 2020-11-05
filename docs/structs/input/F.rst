F
====

Filenames from ``filenames.csv``.

The files are located in ./:ref:`directories/input:input`

.. Note:: This is an array of 22 structs

Initialized
""""""""""""

``SCOPE.m``

Used
"""""
.. list-table::
    :widths: 75 25

    * - variable
      - user
    * - | ``soil_file, leaf_file, atmos_file``
        | ``LIDF_file`` if provided
      - ``SCOPE.m``
    * - ``Simulation_Name``
      - :func:`.create_output_files`
    * - other structs from this
      - :func:`.load_timeseries`

Fields
"""""""

Fields initialized in ``SCOPE.m``. Each of 22 structs in this array has these fields.

.. list-table::
    :widths: 10 10 20 10 50

    * - variable
      - units
      - type
      - default
      - description
    * - **FileID**
      - \-
      - char
      - defined in ``SCOPE.m``
      - SCOPE file identifiers
    * - **FileNames**
      - \-
      - [1 x 1] cell
      - \-
      - filenames from ``filenames``

