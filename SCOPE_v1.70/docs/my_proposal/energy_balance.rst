Energy balance
================

Definition
'''''''''''''

Net radiation (Rn):
.. http://www.indiana.edu/~geog109/topics/04_radiation/4c-RadiationBalance_nf.pdf

>>> Rn = (SW_in - SW_out) + (LW_in - LW_out)

    | SW - shortwave radiation (400-2400 nm)
    | LW - longwave (thermal) radiation

Net radiation consists of 3 (4) heat fluxes:

>>> Rn = H + lE + G

    | H - sensible heat
    | lE - latent heat
    | G - ground heat flux


SCOPE
'''''''

With ``options.calc_ebal`` energy balance loop is started until the energy balance is closed (net radiation become equal to heat fluxes).

To close energy balance leaf temperatures and Monin-Obukhov length ``L`` are iteratively adjusted.

This how it looks like:

.. figure:: ../images/energy_balance.png


Leaf temperatures are calculated by :func:`.biochemical` or :func:`.biochemical_MD12`.

Initial values of soil and leaf temperatures are equal to ambient temperature (``Ta``).

Monin-Obukhov length influences on aerodynamic resistances values.
