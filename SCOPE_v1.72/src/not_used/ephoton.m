function E = ephoton(lambda)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of 
%wavelength lambda (m)

h         = 6.6262E-34; % [J s]         Planck's constant
c         = 299792458;  % [m s-1]       Speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;
