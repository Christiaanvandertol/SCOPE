SCOPE_v1.70
============

.. contents::

Files
'''''''''''''

``input_data.xlsx``
---------------------

Main input file is ``input_data.xlsx`` with 4 sheets. In the documentation we refer to this file, although text alternatives are also possible.

.. Note::
    If Excel is not available, it is possible to use input from text files (.m and .txt). See **alternative**.

    To specify which input to use (text or excel) comment / uncomment lines in ``set_parameter_filenames.m` with ``%`` sign.

.. list-table::
    :widths: 15 70 15
    :header-rows: 1
    :stub-columns: 1

    * - sheet (tab)
      - content
      - alternative
    * - readme
      - | sheets description of ``input_data.xlsx``
        | explanation of leaf inclination distribution function (LIDF) parameters
        | recommended values for plan functional types (PFTs)
        | some parameter ranges
      - \-
    * - options
      - :ref:`options:Options`
      - ``setoptions.m``
    * - filenames
      - filenames for current simulation and for time-series
      - ``filenames.m``
    * - inputdata
      - values for :ref:`structs/input/input_index:input structs`
      - ``inputdata.txt``


Directories
'''''''''''''

output
-------

The function :func:`.output_data` saves the output of SCOPE in an output directory.

In SCOPE, output_data is called after each calculation.

The data are stored in the following directory:
``SiteName_yyyy-mm-dd-hh-mm``

:In which:
    ``yyyy`` refers to the Julian year,

    ``mm`` to the month,

    ``dd`` the day,

    ``hh`` the hour and

    ``mm`` the minutes

    of the time when the simulation was started.

for files see :ref:`outfiles:Output files`

src
----

.m files with the code.

* :ref:`api/equations:+equations`
* :ref:`api/helpers:+helpers`
* :ref:`api/io:+io (input output)`
* :ref:`api/plot:+plot`
* :ref:`api/not_used:not_used`

