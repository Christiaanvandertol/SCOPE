Tricks
========

Prior value selection
-----------------------

If you know your research area well, you might want to "recommend" the retrieval algorithm to stay close to, say, mean values of the parameters.

For this purpose on the ``Input`` tab:

1. Set **value** column to the mean value of your parameter.
2. Set **uncertainty** as the standard deviation of your parameter.
3. Uncomment the line in :func:`COST4SAIL` (currently 88)

.. code-block:: matlab
    :lineno-start: 87
    :caption: COST4SAIL.m

    er2 = 0;
    er2 = (p - prior.Apm) ./ prior.Aps;

    %% total error
    er = [er1 ; 3E-2* er2];  % change value of 3E-2 to higher / lower


Speed up optimization (not recommended)
-----------------------------------------

If you prefer quantity over quality you may provide additional :func:`lsqnonlin` parameters in the :func:`fit_spectra`.

With **stoptol** == 1E-6 one spectra takes ~12 seconds, with **stoptol** == 1E-6 ~3 seconds (currently line 15).

.. code-block:: matlab
    :lineno-start: 15
    :caption: fit_spectra.m

     stoptol = 1E-3;  % we recommend e-6

.. Warning::
    With stoptol 1E-3 we were not able to reproduce with high quality even SCOPE own spectra without any additional noise.

    stoptol = 1E-6:

    .. figure:: ../images/synthetic_ASD_1e-6.png

    stoptol = 1E-3:

    Clearly not all parameters lie on 1 : 1 line anymore.

    .. figure:: ../images/synthetic_ASD_1e-3.png

Synthetic data generation
-----------------------------

Having changed code, input data, introduced new sesnors you might want to know if the retrieval still works properly.

You can generate spectra with known characteristics (forward SCOPE run) with :func:`helpers.generate_synthetic`.

There are already ASD and MSI sensors simulated in folder ``./measured/synthetic``. The new sensors will be saved there as well.

