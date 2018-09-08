function molphotons = e2phot(lambda,E)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

A         = 6.02214E23; % [mol-1]       Constant of Avogadro
e           = ephoton(lambda);
photons     = E./e;
molphotons  = photons./A;
return;