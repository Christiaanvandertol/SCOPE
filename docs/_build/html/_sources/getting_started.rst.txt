Getting started
====================

.. contents::

Software requirements
''''''''''''''''''''''''

The model SCOPE_v1.70 is written in Matlab R2015b running on a Windows operating system. We took care not to use functions that are available in all recent Matlab versions, but we cannot give any warranty that it works under other operating systems and other Matlab versions.

.. warning::
    If you do not have Matlab on your computer you can use ``SCOPE.exe`` with `Matlab Runtime`_ **only R2015b (version 9.0)**


.. _Matlab Runtime: https://nl.mathworks.com/products/compiler/matlab-runtime.html

SCOPE consists of several scripts and functions (modules), which can be used separately or as parts of the integrated SCOPE model (``SCOPE.m``).

When the modules are used separately, then it is important to provide input in the structures specified in :ref:`structs:structs`.

When the integrated model is called, then the input is automatically loaded from the spreadsheet (``input_data.xlsx``) and from the files specified therein (by default located in ../data/:ref:`directories/data:input`.

Basic knowledge of the use of Matlab is required to operate the model.

The application of the model involves the following steps:

1.	Unpack the zip file
-------------------------------
Unpack the model, and **leave the directory structure intact**.

2.	Run the model once
------------------------------
Run the model once, before modifying the parameters and input. It will check whether the software works under your system. The model runs with an example data set (``options.verify``), and the output is automatically compared to output that it should produce. If there is any difference in the results, messages will show up.

The model is executed by opening Matlab, navigating to the directory where the matlab code is (‘cd ./SCOPE_v1.70/:ref:’), and running ‘SCOPE’. Running the model may take a while because almost all options are switched on. If the output of the model is not as expected, then messages will appear. There will also be graphs appearing showing the freshly produced output together with the expected output. If all is ok then no graphs or warnings are produced.

3.	Evaluate and complete the spreadsheet
-----------------------------------------------
‘input_data.xlsx’ or the input files ‘filenames’, ‘inputdata’ and ‘setoptions’
The required input is specified in the spreadsheet file ‘input_data.xlsx’. Open this file. It has three sheets:

-	Readme:  here information about the simulation can be entered
-	Options: specify the simulation options here
-	Filenames:  specify the name of the simulation, the soil and leaf optical property files, and the file names of meteorological input time series.
-	InputData:  specify all the parameters and input variables. Meteorological data specified here will be overwritten by values in the input files if these files have been specified (filenames) and the series option is switched on (options).

If you prefer not to use the spreadsheet, you can provide the input as a text file, and the options and file names as ‘.m’ file (which can be edited with any text editor like notepad). Examples are given. You can set the option whether to use the spreadsheet or text input in the file ‘set_parameter_filenames.m’, by commenting out (with ‘%’) the option that is not wanted. The filenames can be specified here.

4.	Simulation option ‘Individual runs’
---------------------------------------------
The last simulation option is important: to run the model for a few cases only, choose the option: simulation = 0. In that case the model runs for the input specified in the InputData sheet. It is possible to specify more than one value for one input variable, by filling in values in the next column. The model will run as many simulations as there are columns in the input data spreadsheet, say n runs. For run i it will select the data from column i for all variables that have n values. For all other variables, it will select the first value only. For example:
Cab  	10 	20 	30
Cdm 	0.012
N 	1.5 	2
It will do three runs, the first time with Cab = 10, Cdm = 0.012, and N = 1.5;  the second time with Cab = 20, Cdm = 0.012, and N = 1.5;  and a third time with Cab = 30, Cdm = 0.012 and N = 1.5.  The value of N = 2 is ignored and the run cycle ends.
The output is the same as for the time series (see below), except that two additional files are produced: ‘pars_and_input.dat’ and ‘pars_and_input_short.dat’. Both files always have a header. The first file lists the values of all parameters and input variables (that are part of the structure ‘v’) that were used in the simulations, one row for each simulation. The second file lists only the parameters that were varied. Suppose that, for example, if 3 parameters were given 10 different values, while the other parameters were given only 1 single value for each simulation. In that case the pars_and_input_short.dat output file contains three columns with the parameter values corresponding to teach simulations.

5.	Simulation option Time series
------------------------------------------

For the time series run, set simulation = 1. SCOPE now uses the meteorological input as saved in the ascii files specified in the sheet: ‘filenames’. SCOPE runs as many times as there are values in the ascii files. For all input that is not in files, it uses the first value specified in the ‘InputData’ sheet. A value for an input variable in the spreadsheet is overwritten by the value in the time series file of that variable, if this file is provided.
Note: In version 1.70, it is possible to leave a few meteorological input data files blank. In that case, this variable will be the (constant) value in the input data spreadsheet or the general input data text file.

6.	Simulation option Lookup Table
----------------------------------------

For the LUT option, specify ‘simulation = 2’. This option is similar to the ‘individual runs’, except that the model runs over all possible combinations of parameters. For example:
Cab  	10 	20 	30
Cdm 	0.012
N 	1.5 	2
It will do six runs, the first time with Cab = 10, Cdm = 0.012, and N = 1.5;  the second time with Cab = 20, Cdm = 0.012, and N = 1.5;  a third time with Cab = 30, Cdm = 0.012, and N = 1.5;  then fourth with Cab = 10, Cdm = 0.012 and N = 2.0, etc, cycling through values for Cab again.

7.	To execute the model
----------------------------------
The model can be executed by calling ‘SCOPE’ in the command window of Matlab. Alternatively, separate modules can be called, provided that the required input is given. The modules have a help text describing how to do this, which can be called by typing ‘help modulename’, for example: ‘help ebal’. It is however more difficult, because the structures need to be provided.
The output of each simulation is automatically saved in an output directory, together with files documenting the parameters used for this simulation, and the spreadsheet in directory ‘Parameters’.
It is also possible to use the executable SCOPE.exe. In that case you will first need to install the Matlab Runtime Compiler for Matlab 2015b, which can be found on the Mathworks web site.

8.	To plot the output
-------------------------------
An example of a module which creates graphs is provided with the model (plots.m). This function browses through the latest output directory, and plots all data present there in graphs. The titles of the graphs are the headings found in the output files.
