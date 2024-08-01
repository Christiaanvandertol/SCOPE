Getting started
====================

Video tutorials are available on |yt|

.. |yt| image:: ./images/yt_logo_rgb_light.png
   :height: 3ex
   :class: no-scaled-link
   :target: https://youtube.com/playlist?list=PLMKVJ8XOixyTEcUkYDw1YTgKHi-FgimF8&si=E-f3O-j3PrMREL97
   

.. contents::

0. Software requirements
--------------------------

The model SCOPE is written in Matlab R2015b running on a Windows operating system. We took care not to use functions that are available in all recent Matlab versions, but we cannot give any warranty that it works under other operating systems and other Matlab versions.


SCOPE consists of several scripts and functions (modules), which can be used separately or as parts of the integrated SCOPE model (``SCOPE.m``).

When the modules are used separately, then it is important to provide input in the structures specified in :ref:`structs:structs`.

When the integrated model is called, then the input is automatically loaded from the csv and other files located in ./:ref:`directories/input:input`/

Basic knowledge of the use of Matlab is required to operate the model.

The application of the model involves the following steps:


1.	Unpack the zip file
-------------------------------
Unpack the model, and **leave the directory structure intact**.


2.	Run the model once
------------------------------
Run the model once, before modifying the parameters and input. It will check whether the software works under your system. The model runs with an example data set (``options.verify``), and the output is automatically compared to output that it should produce. If there is any difference in the results, messages will show up.

* Navigate to the directory where the ``SCOPE.m`` is (root directory)
* Open ``SCOPE.m`` in Matlab
* in Matlab command window type:
    .. code-block:: matlab

        SCOPE


Running the model may take a while because almost all options are switched on. If the output of the model is not as expected, then messages will appear. There will also be graphs appearing showing the freshly produced output together with the expected output. If all is ok then no graphs or warnings are produced.


3.	Set the input in csv files
---------------------------------------------

Main input files - ``filenames.csv, input_data.scv, setoptions.csv`` (former excel sheets) are located in ``./input``.

.. list-table::
    :widths: 30 70
    :header-rows: 1
    :stub-columns: 1

    * - file
      - content
    * - options
      - :ref:`options:Options`
    * - filenames
      - filenames for current simulation and for time-series
    * - inputdata
      - values for :ref:`structs/input/input_index:input structs`
    * - mSCOPE
      - leaf traits per canopy layer: *optional* only used when ``options.mSCOPE == 1``

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

For any questions, please, use github issues.
