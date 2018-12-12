Getting started
====================

.. contents::

0. Software requirements
--------------------------

The model SCOPE_v1.70 is written in Matlab R2015b running on a Windows operating system. We took care not to use functions that are available in all recent Matlab versions, but we cannot give any warranty that it works under other operating systems and other Matlab versions.

.. warning::
    If you do **not** have Matlab on your computer you can use ``SCOPE.exe`` with `Matlab Runtime`_ **only R2015b (version 9.0)**


.. _Matlab Runtime: https://nl.mathworks.com/products/compiler/matlab-runtime.html

SCOPE consists of several scripts and functions (modules), which can be used separately or as parts of the integrated SCOPE model (``SCOPE.m``).

When the modules are used separately, then it is important to provide input in the structures specified in :ref:`structs:structs`.

When the integrated model is called, then the input is automatically loaded from the spreadsheet :ref:`directories/scope:``input_data.xlsx``` and from the files specified in ./:ref:`directories/data:data`/input.

Basic knowledge of the use of Matlab is required to operate the model.

The application of the model involves the following steps:


1.	Unpack the zip file
-------------------------------
Unpack the model, and **leave the directory structure intact**.


2.	Run the model once
------------------------------
Run the model once, before modifying the parameters and input. It will check whether the software works under your system. The model runs with an example data set (``options.verify``), and the output is automatically compared to output that it should produce. If there is any difference in the results, messages will show up.

* Navigate to the directory where the matlab code is
    ./SCOPE_v1.70/:ref:`directories/scope:src`
* Open ``SCOPE.m`` in Matlab
* in Matlab command window type:
    .. code-block:: matlab

        SCOPE


Running the model may take a while because almost all options are switched on. If the output of the model is not as expected, then messages will appear. There will also be graphs appearing showing the freshly produced output together with the expected output. If all is ok then no graphs or warnings are produced.


3.	Set the input in ``input_data.xlsx``
---------------------------------------------

Main input file i``input_data.xlsx`` with 4 sheets is located in ./SCOPE_v1.70. In the documentation we refer to this file, although text alternatives are also possible.

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

To find out ranges and units of input parameters take a look into :ref:`structs/input/input_index:input structs`.

Pay extra attention to the :ref:`options:``simulation```


4. Analyse the output
-------------------------

All output files and their content (variables, units) are available at :ref:`outfiles:Output files`.

Some output files are available for each run, the others can be written with various :ref:`options:Options`.

To plot the output either select ``options.makeplots`` or use function from :func:`.plots`

.. Note::
    Radiation, spectral and fluorescence output usually has two quantiles:

    * outgoing diffuse light (**hemispherical**) W m-2 um-1
    * outgoing light in observation directions (**directional**, the one that actually reaches the sensor) W m-2 um-1 sr-1

    To get further information see: :ref:`my_proposal/brdf:Definition`


5.	Going further
---------------------------------------------

``SCOPE.m`` is a script, thus after a run all matlab structures that were generated during the run (input, output, constants) are available in the workspace. You can get some extra variables that are not written to output files. You can find out available variables at :ref:`structs:Structs`.

All functions are documented within the code and also at :ref:`api:API`.

For any questions, please, use SCOPE_model SCOPE_model_ group.

.. _SCOPE_model: https://groups.google.com/forum/?fromgroups#!forum/scope_model
