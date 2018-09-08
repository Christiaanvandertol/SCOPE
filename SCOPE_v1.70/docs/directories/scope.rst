SCOPE_v1.70
============

.. contents::

docs
-----

This docs


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

for files see :ref:`output_files`

src
----

.m files with the code.

* :ref:`+equations`
* :ref:`+helpers`
* :ref:`+io (input output)`
* :ref:`+plot`
* :ref:`not_used`

