function rwet = BSM(soilpar,spec,emp)

% Spectral parameters

%wl  = spec.wl;          % wavelengths
GSV = spec.GSV;         % Global Soil Vectors spectra (nwl * 3)
kw  = spec.Kw;          % water absorption spectrum
nw  = spec.nw;          % water refraction index spectrum

% Soil parameters

B   = soilpar.BSMBrightness;        % soil brightness
lat = soilpar.BSMlat;      % spectral shape latitude (range = 20 - 40 deg)
lon = soilpar.BSMlon;      % spectral shape longitude (range = 45 - 65 deg)
SMp = soilpar.SMC;      % soil moisture volume percentage (5 - 55)

% Empirical parameters

SMC  = emp.SMC;         % soil moisture capacity parameter
film = emp.film;        % single water film optical thickness

f1 = B * sind(lat);
f2 = B * cosd(lat) * sind(lon);
f3 = B * cosd(lat) * cosd(lon);

rdry = f1 * GSV(:,1) + f2 * GSV(:,2) + f3 * GSV(:,3);

% Soil moisture effect

rwet = soilwat(rdry,nw,kw,SMp,SMC,film);




function rwet = soilwat(rdry,nw,kw,SMp,SMC,deleff)

    % In this model it is assumed that the water film area is built up  
    % according to a Poisson process. The fractional areas are as follows:
    
    % P(0)   = dry soil area
    % P(1)   = single water film area
    % P(2)   = double water film area
    % ...
    % et cetera
    
    % The fractional areas are given by P(k) = mu^k * exp(-mu) / k! 
    
    % For water films of multiple thickness only the transmission loss due
    % to water absorption is modified, since surface reflectance effects 
    % are not influenced by the thickness of the film
    
    % Input parameters:
    
    % rdry   = dry soil reflectance                             [NW,1]
    % nw     = refraction index of water                        [NW,1]
    % kw     = absorption coefficient of water                  [NW,1]
    % SMp    = soil moisture volume percentage                  [1,NS]
    % SMC    = soil moisture capacity (recommended 0.25)        [1,1]
    % deleff = effective optical thickness of single water film [1,1]
    %          (recommended 0.015)
    
    % Output
    
    % rwet   = wet soil spectra                                 [NW,NS]
    
    % If SMp is given as a row-vector and rdry, nw, kw as column vectors 
    % of the same size (NW, # of wavelengths), then the output is a matrix 
    % of spectra for the different SMp, where each column is a spectrum
    
    % Wout Verhoef
    % Version 1.0
    % September 2012
    
    %---------------------------------------------------------------------%
    
    % two-way transmittance of elementary water film
    
    tw = exp(-kw * deleff);
    
    % Lekner & Dorf (1988) modified soil background reflectance
    % for soil refraction index = 2.0; uses the tav-function of PROSPECT
    
    rbac = 1 - (1-rdry) .* (rdry .* tav(90,2.0./nw) / tav(90,2.0) + 1-rdry);
    
    % total reflectance at bottom of water film surface
    
    p    = 1 - tav(90,nw) ./ nw.^2;
    
    % reflectance of water film top surface, use 40 degrees incidence angle, 
    % like in PROSPECT
    
    Rw  = 1 - tav(40,nw);
    
    % additional reflectance of single water film (Lekner & Dorf, 1988)
    % two-way transmission loss by water absorption is not included here
    % yet
    
    Radd   = (1-Rw) .* (1-p) .* rbac ./(1 - p .* rbac);
    
    % Mu-parameter of Poisson distribution
    
    mu  = (SMp - 5)/ SMC;
    
    % fraction of dry soil area
    
    fdry = exp(-mu);
    
    % contribution due to total water film area of single or 
    % multiple thickness
    
    fmul = (exp(tw * mu) - 1) * diag(fdry);
    
    % reflectance spectra of wet soil

    rwet = rdry * fdry + Rw * (1 - fdry) + Radd * ones(size(mu)) .* fmul;
    
return


