BRDF
=====

``options.calc_directional`` & ``options.calc_ebal``

.. Warning:: This is an advanced topic, please refer to Schaepman-Strub et al. 2006 :cite:`Schaepman-Strub2006` for further explanations.

Definition
''''''''''''

Light consists of two components **direct** (aka specular) and **diffuse** (aka hemispherical).

To explain reflectance of each light component individually, different reflectance factors are used.

SCOPE model simulates the following reflectance factors:

* Incoming light is directional
    * CASE 1: bidirectional (BRF)
    * CASE 2: directional-hemispherical (DHRF)
* Incoming light is hemispherical
    * CASE 7: hemispherical-directional (HDRF)
    * CASE 9: bihemispherical (HRF)

After reflectance from a material *direct component* of incoming light contributes to both directional and hemispherical component of reflected light.

After reflectance from a material *diffuse component* of incoming light also contributes to both directional and hemispherical component of reflected light.

.. figure:: ../images/reflectance_factors.png

    From Schaepman-Strub et al. 2006 :cite:`Schaepman-Strub2006`.

.. Note::
    **Bidirectional Reflectance Distribution Function** is a function describing bidirectional reflectance from a material.

    - "input" of BRDF are four angles (solar zenith and azimuth angle (direction of incoming light); viewing zenith and azimuth angle (direction of observation))
    - "output" of BRDF is reflectance (BRF)

SCOPE
'''''''

To simulate BRDF enable ``options.calc_directional``.

SCOPE calculates BRDF itself and also directional fluorescence radiance and directional thermal radiance (or brightness temperature).

Directional plots have 3 components:

    * viewing zenith angle (towards the centre of the circle)
    * viewing azimuth angle (around the circle)
    * measured quantity (color)

On all graphs you can see a **hot spot** (red dot) where viewing azimuth angle is 0º and *viewing* zenith angle is 30º.

**Hot spot** occurs when the observation direction coincides with the illumination direction. Indeed, for this example *solar* zenith angle of 30º was used.

Directional plots are made per wavelength.

.. figure:: ../images/BRDF.png

    Courtesy of Peiqi Yang
