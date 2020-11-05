Fluorescence
==============

``options.calc_fluor``

Definition
''''''''''''

Light reaching a leaf (pretty much any object) can be reflected, transmitted or absorbed. In plants absorbed light can be spent on three different processes:

#. photochemistry (assimilation CO2)
#. non-photochemical quenching (NPQ): heat dissipation
#. chlorophyll fluorescence excitation

.. figure:: ../images/aPAR_destiny.png
    :scale: 50 %

    Source: http://spie.org/newsroom/4725-remote-sensing-of-terrestrial-chlorophyll-fluorescence-from-space?SSO=1

>>> Fluorescence is light emitted by chlorophyll molecules in the range 640-800 nm.

SCOPE
'''''''

Fluorescence light (pretty much as any light) can be absorbed and scattered on its way to a sensor.

SCOPE model simulates 3 hemispherical fluorescence quantiles:

 #. fluorescence emitted by all photosystems without any scattering / re-absorption neither within leaf, nor withing canopy
 #. fluorescence emitted by all leaves without any scattering / re-absorption within canopy or from soil
 #. fluorescence emitted by canopy (all leaf layers) accounting for all scattering / re-absorption events

.. figure:: ../images/fluorescence.png

.. Note:: Notice the difference in ranges and units between directional and hemispherical fluorescence.

SCOPE model simulates directional fluorescence (the one that actually reaches a sensor) and its components coming from:

 #. sunlit leaves
 #. shaded leaves
 #. scattered by leaves and soil

It is also possible to partition directional fruorescence between photosystem I and II (PSI, PSII) with ``options.calc_PSI``, not recommended though.

.. figure:: ../images/fluorescence_contributors.png

.. Note:: There are much more outputs of :func:`.biochemical` related to Pulse-Amplitude-Modulation (PAM) Fluorometry quantiles that are stored in internal structure ``biochem_out``.
