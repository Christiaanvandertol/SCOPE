function [rad] = RTMt_sb(spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Tsu,Tsh,obsdir)

% function 'RTMt_sb' calculates total outgoing radiation in hemispherical
% direction and total absorbed radiation per leaf and soil component.
% Radiation is integrated over the whole thermal spectrum with
% Stefan-Boltzman's equation. This function is a simplified version of
% 'RTMt_planck', and is less time consuming since it does not do the
% calculation for each wavelength separately.
%
% Authors: Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
% date:     5  Nov 2007
% update:   13 Nov 2007
%           16 Nov 2007 CvdT    improved calculation of net radiation
%           27 Mar 2008 JT      added directional calculation of radiation
%           24 Apr 2008 JT      Introduced dx as thickness of layer (see parameters)
%           31 Oct 2008 JT      introduced optional directional calculation
%           31 Oct 2008 JT      changed initialisation of F1 and F2 -> zeros
%           07 Nov 2008 CvdT    changed layout
%           16 Mar 2009 CvdT    removed Tbright calculation
%              Feb 2013 WV      introduces structures for version 1.40
%
% Table of contents of the function
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     parameters       
%       0.3     geometric factors of Observer
%       0.4     geometric factors associated with extinction and scattering
%       0.5     geometric factors to be used later with rho and tau
%       0.6     fo for all leaf angle/azumith classes
%   1       calculation of upward and downward fluxes          
%   2       total net fluxes
%   Appendix A. Stefan-Boltzmann
%
% usage:
% [rad] = RTMt_sb(options,spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Tsu,Tsh)
%
% Most input and output are structures. These structures are further
% specified in a readme file. The temperatures Tcu, Tch, Tsu and Tsh are
% variables.
%
% Input:
%   options       calculation options
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%   Tcu         Temperature of sunlit leaves    (oC), [13x36x60]
%   Tch         Temperature of shaded leaves    (oC), [13x36x60]
%   Tsu         Temperature of sunlit soil      (oC), [1]
%   Tsh         Temperature of shaded soil      (oC), [1]
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients.
%               Here, thermal fluxes are added
%% 0.0 globals
global constants

%% 0.1 parameters

IT          = find(spectral.wlS == 10000);   % Take 10 microns as representative wavelength for the thermal

deg2rad     = constants.deg2rad;
nl          = canopy.nlayers;
lidf        = canopy.lidf;
litab       = canopy.litab;
lazitab     = canopy.lazitab;
nlazi       = length(lazitab);
tto         = angles.tto;
psi         = angles.psi;
Ps          = gap.Ps;
K           = gap.K;

rho         = leafopt.refl(IT);       % [nwl]               Leaf/needle reflection
tau         = leafopt.tran(IT);       % [nwl]               Leaf/needle transmission
rs          = soil.refl(IT);          % [nwl]               Soil reflectance
epsc        = 1-rho-tau;              % [nwl]               Emissivity vegetation
epss        = 1-rs;                   % [nwl]               Emissivity soil
crit        = max(1E-2);              % [1]                 Desired minimum accuracy
LAI         = canopy.LAI;
dx          = 1/nl;
iLAI        = LAI*dx;

%% 0.2 initialiations
Rnhc        = zeros(nl,1);  % [nl]
Rnuc        = zeros(size(Tcu));   % [13,36,nl]

%% 0.3 geometric factors of observer
if obsdir
    cos_tto = cos(tto*deg2rad);   % [1]                 cos observation angle
    sin_tto = sin(tto*deg2rad);   % [1]                 sin observation angle
end

%% 0.4 geometric factors associated with extinction and scattering
cos_ttl     = cos(litab*deg2rad);
if obsdir
    sin_ttl = sin(litab*deg2rad);
    cos_ttlo= cos((lazitab-psi)*deg2rad);
end
bfli        = cos_ttl.^2;
bf          = bfli'*lidf;

%% 0.5 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
ddb         = 0.5*(1+bf);         %                     f1^2 + f2^2
ddf         = 0.5*(1-bf);         %                     2*f1*f2
if obsdir
    dob     = 0.5*(K+bf);         %                     fo*f1
    dof     = 0.5*(K-bf);         %                     fo*f1
end

%% 0.6 fo for all leaf angle/azumith classes
if obsdir
    Co      = cos_ttl*cos_tto;    % [nli]               pag 305
    So      = sin_ttl*sin_tto;    % [nli]               pag 305
    cos_deltao  = Co*ones(1,nlazi) + So*cos_ttlo;% [nli, nlazi] projection of leaves in  in direction of sun (pag 125-126)
    fo      = cos_deltao/abs(cos_tto);% [nli, nlazi]    leaf area projection factors in direction of observation
end

%% 1. calculation of upward and downward fluxes pag 305
sigb        = ddb*rho + ddf*tau;  % [nwlt]              Diffuse backscatter scattering coefficient
sigf        = ddf*rho + ddb*tau;  % [nwlt]              Diffuse forward     scattering coefficient
if obsdir
    vb      = dob*rho + dof*tau;  % [nwlt]              Directional backscatter scattering coefficient for diffuse  incidence
    vf      = dof*rho + dob*tau;  % [nwlt]              Directional forward     scattering coefficient for diffuse  incidence
end
a           = 1-sigf;             % [nwlt]              Attenuation
m           = sqrt(a*a-sigb*sigb);% [nwlt]
rinf        = (a-m)/sigb;         % [nwlt]              Reflection coefficient for infinite thick canopy    
rinf2       = rinf*rinf;          % [nwlt]

fHs         = (1-rinf2)*(1-rs)/(1-rinf*rs);
fHc         = iLAI*m*(1-rinf);
fbottom     = (rs-rinf)/(1-rs*rinf);

%1.1 radiance by components
Hcsu3       = Stefan_Boltzmann(Tcu);%                   Radiance by sunlit leaves
Hcsh        = Stefan_Boltzmann(Tch);%                   Radiance by shaded leaves
Hssu        = Stefan_Boltzmann(Tsu);%                   Radiance by sunlit soil
Hssh        = Stefan_Boltzmann(Tsh);%                   Radiance by shaded soil

% 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)
v1 = repmat( 1/size(Hcsu3, 2), 1, size(Hcsu3, 2)); % vector for computing the mean
Hcsu2 = reshape(Hcsu3, size(Hcsu3, 1), []);   % create a block matrix from the 3D array
Hcsu = (v1 * reshape(Hcsu2'*lidf, size(Hcsu3, 2), []))'; % compute column means for each level

Hc          = Hcsu.*Ps(1:nl) + Hcsh.*(1-Ps(1:nl));      % hemispherical emittance by leaf layers
Hs          = Hssu.*Ps(nl+1) + Hssh.*(1-Ps(nl+1));      % hemispherical emittance by soil surface

% 1.3 Diffuse radiation
cont        = 1;                                        % continue iteration (1:yes, 0:no)
counter     = 0;                                        % number of iterations
F1          = zeros(nl+1,1);
F2          = zeros(nl+1,1);
F1top       = 0;
while cont
    F1topn  = -rinf*F2(1);
    F1(1)   = F1topn;
    for j   = 1:nl
        F1(j+1) = F1(j)*(1-m*iLAI)+ fHc*Hc(j);
    end
    F2(nl+1) = fbottom*F1(nl+1) + fHs*Hs;
    for j   = nl:-1:1
        F2(j)   = F2(j+1)*(1-m*iLAI) + fHc*Hc(j);
    end
    cont    = abs(F1topn-F1top)>crit;
    F1top   = F1topn;
    counter = counter + 1;
end

Emin        = (F1+rinf*F2)/(1-rinf2);
Eplu        = (F2+rinf*F1)/(1-rinf2);

% 1.4 Directional radiation
if obsdir
    piLo1       = iLAI*epsc*K*Hcsh'*(gap.Po(1:nl)-gap.Pso(1:nl));                % directional   emitted     radation by shaded leaves
    % JAK 2015-01: replaced earlier loop by this: all-at-once with more efficient mean
    absfo_rep = repmat(abs(fo), 1, nl);
    piLo2 = iLAI*epsc*(v1 * reshape( (Hcsu2.*absfo_rep)'*lidf, size(Hcsu3, 2), []))'.*gap.Pso(1:nl); % compute column means for each level

    piLo3       = iLAI*((vb*Emin(1:nl) + vf*Eplu(1:nl))'*gap.Po(1:nl));      % directional   scattered   radiation by vegetation for diffuse incidence
    piLo4       = epss*Hssh*(gap.Po(nl+1)-gap.Pso(nl+1));                        % directional   emitted     radiation by shaded soil
    piLo5       = epss*Hssu*gap.Pso(nl+1);                                   % directional   emitted     radiation by sunlit soil
    piLo6       = rs*Emin(nl+1)*gap.Po(nl+1);                                % directional   scattered   radiation by soil       for diffuse incidence       [1]
    
    piLot       = piLo1 + sum(piLo2) + piLo3 + piLo4 + piLo5 + piLo6;
else
    piLot       = NaN;
end
Lot             = piLot/pi;

%% 2. total net fluxes
% net radiation per component, in W m-2 (leaf or soil surface)
for j = 1:nl
    Rnuc(:,:,j) = (Emin(j) + Eplu(j+1) - 2*Hcsu3(:,:,j))*epsc;    % sunlit leaf
    Rnhc(j)     = (Emin(j) + Eplu(j+1) - 2*Hcsh(j))*epsc;         % shaded leaf
end
Rnus            = (Emin(nl+1) - Hssu)*epss;                       % sunlit soil
Rnhs            = (Emin(nl+1) - Hssh)*epss;                       % shaded soil

%% 3. Write the output to the rad structure

rad.Emint   = Emin;
rad.Eplut   = Eplu;
rad.Eoutte  = Eplu(1)-Emin(1); % 1)  
rad.Lot     = Lot;
rad.Rnuct   = Rnuc;
rad.Rnhct   = Rnhc;
rad.Rnust   = Rnus;
rad.Rnhst   = Rnhs;
return

% 1) CvdT, 11 December 2015. 
% We subtract Emin(1), because ALL incident (thermal) radiation from Modtran 
% has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!

%% Appendix A. Stefan-Boltzmann
function H      =   Stefan_Boltzmann(T_C)

global constants;
C2K     = constants.C2K;
sigmaSB = constants.sigmaSB;

H       = sigmaSB*(T_C + C2K).^4;
return