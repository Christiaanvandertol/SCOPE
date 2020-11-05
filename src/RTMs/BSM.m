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
    SMp = soilpar.SMC*1E2;      % soil moisture volume percentage (5 - 55)
    
    % Empirical parameters
    
    SMC  = emp.SMC;         % soil moisture capacity parameter
    film = emp.film;        % single water film optical thickness
    
    f1 = B * sind(lat);
    f2 = B * cosd(lat) * sind(lon);
    f3 = B * cosd(lat) * cosd(lon);
    
    rdry = f1 * GSV(:,1) + f2 * GSV(:,2) + f3 * GSV(:,3);
    
    % Soil moisture effect
    
    rwet = soilwat(rdry,nw,kw,SMp,SMC,film);
    
return 

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
    
% Peiqi Yang
% Version 1.1
% Jan 2020
% Note by PY (p.yang@utwente.nl)
% because mu<=2, P(k>6) is negligible

%---------------------------------------------------------------------%
k       =   0:6;                    % number of water film, '0' refers to dry soil
nk      =   length(k);              % the number of occurrences
mu      =   (SMp - 5)/ SMC;         % Mu-parameter of Poisson distribution
if mu   <=  0                       % the reason for adding this: if mu<0, fry>1.
    rwet = rdry;                    % we need to check SMC in other parts of SCOPE. soil fluxes routine.
else
    
    % Lekner & Dorf (1988) modified soil background reflectance
    % for soil refraction index = 2.0; uses the tav-function of PROSPECT
    rbac = 1 - (1-rdry) .* (rdry .* tav(90,2.0./nw) / tav(90,2.0) + 1-rdry); % Rbac
    
    % total reflectance at bottom of water film surface
    p    = 1 - tav(90,nw) ./ nw.^2;   % rho21, water to air, diffuse
    
    % reflectance of water film top surface, use 40 degrees incidence angle,
    % like in PROSPECT
    Rw  = 1 - tav(40,nw);             % rho12, air to water, direct
  
    
    % fraction of areas
    % P(0)   = dry soil area            fmul(1)
    % P(1)   = single water film area   fmul(2)
    % P(2)   = double water film area   fmul(3)
    % without loop
    fmul    =   poisspdf(k,mu)';                          % Pobability 
    tw      =   exp(-2*kw * deleff.*k);                   % two-way transmittance,exp(-2*kw*k Delta)
    Rwet_k  =   Rw + (1-Rw) .* (1-p) .*tw.* rbac ./(1 - p .*tw.* rbac);
    rwet   =   rdry * fmul(1) + Rwet_k(:,2:nk)*fmul(2:nk);
end     
return

function Tav = tav(alfa,nr)
n2                                  =   nr.^2;
np                                  =   n2 + 1;
nm                                  =   n2 - 1;

% Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees

% y1 = (3*n2+2*nr+1)./(3*(nr+1).^2);
% y2 = 2*nr.^3.*(nr.^2+2*nr-1)./(np.^2.*nm);
% y3 = n2.*np.*log(nr)./nm.^2;
% y4 = n2.*nm.^2./np.^3.*log((nr.*(nr+1)./(nr-1)));

% st = y1-y2+y3-y4;

a                                   =   +((nr+1).^2)/2;
k                                   =   -((n2-1).^2)/4;
sin_a                               =   sind(alfa);
%
if alfa~=0    
    B2                              =   sin_a^2 - np/2;
    B1                              =   (alfa~=90) * sqrt( B2.^2 + k );
   
    b                               =   B1 - B2;
    b3                              =   b.^3;
    a3                              =   a.^3;
    
    ts                              =   (k.^2./(6*b3) + k./b - b./2) - ...
                                        (k.^2./(6*a3) + k./a - a./2);
                                    
    tp1                             =   -2*n2.*    (   b  -  a   ) ./ (np.^2);
    tp2                             =   -2*n2.*np.*(  log(b./a)  ) ./ (nm.^2);
    tp3                             =      n2.*    ( 1./b - 1./a ) ./ 2; 
    
%     tp4                             =   16*n2.^2.* (n2.^2+1) .* ( log(2*np.*b - nm.^2) - log(2*np.*a - nm.^2) ) ./ (np.^3.*nm.^2);    
%     tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./ (np.^3       );

    tp4                             =	16*n2.^2.* (n2.^2+1) .* ( log((2*np.*b - nm.^2)./(2*np.*a - nm.^2))  ) ./(np.^3.*nm.^2);
    tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./(np.^3);							 
    tp                              =   tp1 + tp2 + tp3 + tp4 + tp5;
    Tav                             =   (ts + tp) ./ (2*sin_a.^2);
else
    Tav                             =   4 *nr/((nr+1)*(nr+1));
end
return

