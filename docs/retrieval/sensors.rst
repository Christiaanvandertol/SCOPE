Custom sensors
================

Sensor collects radiance with a certain spectral sampling interval (SSI) and defined spectral response function (SRF).
In the simples case SRF is just a gaussian curve with the centre at the band wavelength and width defined with full width half maximum (FWHM).

We convolve irradiance from **atmfile** (resolution 0.01 nm) (or **Esun, Esky** (resolution 1 nm)) and simulated reflected radiance to the sensor parameters.

Hyperspectral
---------------

Out of the box hyperspectral sensors are:

1. ASD
2. HyPlant2018
3. FLOX_SPEC2
4. A sensor with constant FWHM over whole measurement wavelength

Hyperspectral instruments typically have individual sensors for visible and SWIR regions with different FWHM.
However, taken that SSI is < 3 nm and FWHM is relatively small, we do not think it is necessary to reconstruct the gaussian curve for each spectral region.

.. figure:: ../images/FWHM.png

    From this figure you can see that FWHM is important for accurate radiance simulation, but not so important for reflectance.

Multispectral
---------------

Out of the box multispectral sensors are:

1. Sentinel-3:
    * S3A_OLCI
    * SLSTR
2. Sentinel-2:
    * MSI


.. warning::
    After adding SRF of sensor to ``input/sensors.xlsx`` add its name into :func:`io.read_fixed_input` (line 22)

    .. code-block:: matlab
        :caption: +io/read_fixed_input.m
        :lineno-start: 16
        :emphasize-lines: 6

            %% collect

            fixed.spectral = spectral;
            fixed.pcf = pcf;
            fixed.optipar = optipar.optipar;
            fixed.srf_sensors = {'S3A_OLCI', 'MSI', 'YOUR_NEW_SENSOR_NAME_AS_SHEET_NAME_IN_SENSORS_XLSX'};
        end


For satellites `SRFs are known`_ and can be added in ``input/sensors.xlsx``.

.. _`SRFs are known`: https://www.nwpsaf.eu/site/software/rttov/download/coefficients/spectral-response-functions/

For other multispectral sensors we recommend reconstructing the SRF as a gaussian curve from FWHM with :func:`helpers.create_sensor_from_fwhm(input_path)`.

input_path - path to an excel file with sheet 'fwhm' with two columns - band centre, fwhm.

The function will write the second sheet 'sensor' that can be copied to ``input/sensors.xlsx``.
The function will also plot the resulting SRFs.

The example of input_path is in ``input/bands_from_fwhm.xlsx``.

.. Note::
        Having introduced new sensor, you might want to perform a sensitivity analysis with :ref:`retrieval/tricks:Synthetic data generation`
        to find out which parameters can be retrieved with your sensor.