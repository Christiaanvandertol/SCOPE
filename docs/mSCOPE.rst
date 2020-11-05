mSCOPE
############

mSCOPE is an integrated **multi-layer** model of vegetation reflectance, photosynthesis, fluorescence, temperature and energy balance.

The code and detailed information are available at: https://github.com/peiqiyang/mSCOPE

Brief introduction
''''''''''''''''''''
The original SCOPE model assumes that vegetation canopies are vertically homogeneous and horizontally infinite, as its radiative transfer routines are based on the classical 1-D SAIL model. However, in reality, **canopies generally exhibit large vertical heterogeneity of both biophysical and biochemical properties**. The development of mSCOPE is to incorporate the vertical variations of vegetation properties. Therefore, the model can be considered as a 2-D model since it does not consider the horizontal variations. It is noted that mSCOPE works for homogeneous canopies as well by just setting all the layers identical or using one layer. The original paper on mSCOPE :cite:`Yang2017`.

This version of the mSCOPE model (mSCOPE_v1_beta) is based upon the SCOPE model (v1.61). It simulates the light interaction and energy balance of vertically heterogeneous canopies.
mSCOPE keeps the same model structure and output of SCOPE, but uses a different solution for the radiative transfer of incident and emitted radiation in vegetation canopies.
Questions related to mSCOPE, please contact p.yang@utwente.nl; or peiqiyangweb@gmail.com (Peiqi Yang).

Prerequisites
''''''''''''''''

.. Warning::
    This model works only with **Matlab2017a** and newer due to matrix multiplication issues.

    The input data are structured in a excel spreadsheet. The users are expected to have the Microsoft installed.

If there is not, please contact p.yang@utwente.nl. This issue can be easily fixed by giving alternative  input options.

What will you receive when you download the model
''''''''''''''''''''''''''''''''''''''''''''''''''''

In the main folder, you will find

- manuals
    - SCOPE manuals
    - SCOPE and mSCOPE presentations

- inputdata
    - input is the same as the original SCOPE: :ref:`directories/input:input`

- output
    .. Note:: the output is saved in the directory

- mSCOPE_code_v1 (similar to :ref:`api:api`)
    -  ``Input_data.xlsx``
        You change the input parameters in this file.

        .. Warning::
            Time-series simulations and LUT simulations are not possible yet (see :ref:`options:Rules of input reading`)

    -  mSCOPE.m
        The main function

    - RTMs - leaf and canopy radiative transfer models
        -	**fluspect_b**  simulating leaf reflectance, transmittance and fluorescence emission matrices
        -	**RTMo_m** canopy radiative transfer in the solar domain.
        - 	**RTMf_m** canopy radiative transfer model for fluorescence
        - 	**RTMt** canopy radiative transfer for emitted thermal radiation

    - Supporting
        -	**Brightness_T** converting radiant emittance(energy per time per area) to temperature for blackbody by inverting Stefan–Boltzmann law
        -  **calczenithangle**	computing solar position based on the location and time
        - 	**e2phot**	calculating the number of moles of photons corresponding to E Joules of energy of wavelength lambda
        -  **ephoton** calculating the energy content (J) of 1 photon of wavelength lambda (m)

    - IO
        - reading input data
        - exporting and plotting (optional) simulation results

    - Fluxes
        - Computing photosynthesis, latent and sensible heat, leaf temperature


Summary of the main changes in mSCOPE
''''''''''''''''''''''''''''''''''''''

1. input_mSCOPE.m
    It reads the vertical profiles of leaf optical properties (e.g. Cab, Cw) from input_data.xlsx in spratsheet 'mSCOPE'
    The input_mSCOPE is called in the main function mSCOPE.m L61 before executing fluspect_mSCOPE and canopy RTMs.
2. fluspect_mSCOPE.m
    It runs fluspect_b for different layers to obtain leaf reflectance, transmittance, Mb and Mf fluspect_mSCOPE is called in the main function mSCOPE.m L249
3. RTMo_m.m
    It is a replacement of the RTMo.m in SCOPE. Many changes have been made here.
    RTMo_m.m is called in the main function mSCOPE.m L279
4. RTMf_m.m
    It is a replacement of the RTMf.m in SCOPE. Many changes have been made here.
    RTMf_m is called in the main function mSCOPE.m L279

References
''''''''''''

:cite:`Yang2017` Yang, P., Verhoef, W., & Van Der Tol, C. (2017). The mSCOPE model: A simple adaptation to the SCOPE model to describe reflectance, fluorescence and photosynthesis of vertically heterogeneous canopies. Remote sensing of environment, 201, 1-11.

Authors
'''''''''

Peiqi Yang (p.yang@utwente.nl; peiqiyangweb@gmail.com)

Wout Verhoef  (w.verhoef@utwente.nl)

Christiaan van der Tol (c.vandertol@utwente.nl)

License
''''''''''
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

