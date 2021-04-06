function E = ephoton(lambda,constants)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
%wavelength lambda (m)

h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;