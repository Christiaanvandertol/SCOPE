function molphotons = e2phot(lambda,E,constants)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

e           = ephoton(lambda,constants);
photons     = E./e;
molphotons  = photons./constants.A;
return;