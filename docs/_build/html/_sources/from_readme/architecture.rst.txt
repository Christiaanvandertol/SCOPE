Model architecture
=======================

After reading in the main input file (a spreadsheet with tabs or alternatively, three text files), the supporting data are loaded. These consist of soil and leaf optical coefficients, and (optional) meteorological time series. The model performs a number of simulations. Each simulation starts with the leaf optical model FLUSPECT and the ratiative transfer model RTMo for scaling from leaf to canopy. Next, the energy balance can be calculated. The energy balance routine calculates the turbulent heat fluxes, photosynthesis, and the leaf and soil temperatures (for each leaf orientation, leaf layer, and for the sunlit and shaded fraction separately). The temperatures are solved by iteration until the energy balance closure error is small enough. The optional calculations that may follow include the radiative transfer of fluorescence, the thermal emission spectra (2.5- 50 μm), and the calculation of a complete BRDF using simulations of many observation angles.

Pseudocode:

.. code-block:: python

    Load input files: filenames, options and input data
    Loop over number of simulations
    Simulate leaf spectra (Fluspect)
    Simulate radiative transfer of incident light (RTMo)
        If option 'calculate energy balance' is on
            Repeat
                Simulate radiative transfer of emitted thermal radiation (RMTt)
                Simulate non-radiative energy balance fluxes (heatfluxes)
            Until the energy balance closes
        end
        If option 'calculate fluorescence' is on
            Simulate radiative transfer of emitted fluorescence
        end
        If option 'calculate Xanthophyll reflectance' is on
                Simulate radiative transfer of incident light due to change in reflectance per leaf (RTMz)
        end
        If option 'calculate Directional' is on
            Repeat the calculations above for many observation zenith and azimuth angles
        end
    Store the output
    end
    Optional: make plots of the output
    Optional: compare the output to a reference simulation

