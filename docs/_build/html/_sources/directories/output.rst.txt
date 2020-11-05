Output
========

In this directory all output will be saved.

It has also ``verificationdata`` example output, against which the very first run is compared (``options.verification == 1``)


The function :func:`.output_data_binary` saves the output of SCOPE in an output directory.

In SCOPE2, output_data_binary is called after each calculation.

When all computations are completed binary files are converted to csv :func:`.bin_to_csv` if ``options.save_csv == 1``

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
