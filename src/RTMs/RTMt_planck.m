function [rad] = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,Tcu,Tch,Tsu,Tsh)

% function 'RTMt_sb' calculates total outgoing radiation in hemispherical
% direction and total absorbed radiation per leaf and soil component.
% Radiation is integrated over the whole thermal spectrum with
% Stefan-Boltzman's equation. This function is a simplified version of
% 'RTMt_planck', and is less time consuming since it does not do the
% calculation for each wavelength separately.
%
% Authors: Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl)
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
%           04 Dec 2019 CvdT    adapted for SCOPE-lite
%           17 Mar 2020 CvdT    mSCOPE representation
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
% [rad] = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Tsu,Tsh,obsdir,spectral)
%
% Most input and output are structures. These structures are further
% specified in a readme file. The temperatures Tcu, Tch, Tsu and Tsh are
% variables.
%
% Input:
%   constants   physical constants
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   Tcu         Temperature of sunlit leaves    (oC), [13x36x60]
%   Tch         Temperature of shaded leaves    (oC), [13x36x60]
%   Tsu         Temperature of sunlit soil      (oC), [1]
%   Tsh         Temperature of shaded soil      (oC), [1]
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, thermal fluxes are added

%% 0.1 parameters
IT          = spectral.IwlT;   %
wlt         = spectral.wlT;
%deg2rad     = constants.deg2rad;

nl          = canopy.nlayers;
lidf        = canopy.lidf;
Ps          = gap.Ps;
%
rho         = leafopt.refl(:, IT)';    % [1]               Leaf/needle reflection
tau         = leafopt.tran(:, IT)';    % [1]               Leaf/needle transmission
rs          = soil.refl(IT);        % [1]               Soil reflectance
epsc        = 1-rho-tau;              % [nwl]               Emissivity vegetation
epss        = 1-rs;                   % [nwl]               Emissivity soil
LAI         = canopy.LAI;
dx          = 1/nl;
iLAI        = LAI*dx;

Xdd         = rad.Xdd(:,IT);
Xsd         = rad.Xsd(:,IT);
Xss         = repmat(rad.Xss,canopy.nlayers,1);
R_dd        = rad.R_dd(:,IT);
R_sd        = rad.R_sd(:,IT);
rho_dd      = rad.rho_dd(:,IT);
tau_dd      = rad.tau_dd(:,IT);

%% 0.2  initialization of output variables
[piLot_,Eoutte_]    = deal(zeros(1,length(IT))); %          [1,nwlt]
[Emin_,Eplu_]       = deal(zeros(nl+1,length(IT)));       % [nl+1,nwlt]

%% 1. calculation of upward and downward fluxes pag 305

for i = 1:length(IT)
    % 1.1 radiance by components
    Hcsu3          = pi*Planck(wlt(i),Tcu+273.15,epsc(i));
    Hcsh           = pi*Planck(wlt(i),Tch+273.15,epsc(i));
    Hssu           = pi*Planck(wlt(i),Tsu+273.15,epss(i));
    Hssh           = pi*Planck(wlt(i),Tsh+273.15,epss(i));
    % 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)
    if size(Hcsu3,2)>1
        v1 = repmat( 1/size(Hcsu3, 2), 1, size(Hcsu3, 2)); % vector for computing the mean
        Hcsu2 = reshape(Hcsu3, size(Hcsu3, 1), []);   % create a block matrix from the 3D array
        Hcsu = (v1 * reshape(Hcsu2'*lidf, size(Hcsu3, 2), []))'; % compute column means for each level
    else
        Hcsu = Hcsu3;
    end
    Hc          = Hcsu.*Ps(1:nl) + Hcsh.*(1-Ps(1:nl));      % hemispherical emittance by leaf layers
    Hs          = Hssu.*Ps(nl+1) + Hssh.*(1-Ps(nl+1));      % hemispherical emittance by soil surface
    
    % 1.3 Diffuse radiation
    [U,Es_,Emin,Eplu]           = deal(zeros(nl+1,1));       % [nl+1,nwl]     direct, up and down diff. rad.
    
    U(nl+1)               =   Hs;
    Es_(1)              =   0;
    Emin(1)            =   0;
    
    for j=nl:-1:1      % from bottom to top
        Y(j)  =   (rho_dd(j).*U(j+1)+Hc(j)*iLAI)./(1-rho_dd(j).*R_dd(j+1));
        U(j)  =   tau_dd(j)*(R_dd(j+1).*Y(j)+U(j+1))+Hc(j)*iLAI;
    end
    for j=1:nl       % from top to bottom
        Es_(j+1)    =   Xss(j).*Es_(j);
        Emin(j+1)  =   Xsd(j).*Es_(j)+Xdd(j).*Emin(j)+Y(j);
        Eplu(j)    =   R_sd(j).*Es_(j)+R_dd(j).*Emin(j)+U(j);
    end
    Eplu(nl+1)    =   R_sd(nl).*Es_(nl)+R_dd(nl).*Emin(nl)+Hs;
    
    Emin_(:,i)      = Emin;                        %               downwelling diffuse radiance per layer
    Eplu_(:,i)      = Eplu;                        %               upwelling   diffuse radiance
    
    Eoutte_(i)      = Eplu(1);
    
    % 1.4 Directional radiation and brightness temperature
    K           = gap.K;
    vb          = rad.vb(end);
    vf          = rad.vf(end);
    piLov       = iLAI*...
       (K*Hcsh'*(gap.Po(1:nl)-gap.Pso(1:nl))+  ...              % directional   emitted     radation by shaded leaves
       K*Hcsu'*gap.Pso(1:nl)+ ... % compute column means for each level
       (vb*Emin(1:nl) + vf*Eplu(1:nl))'*gap.Po(1:nl));      % directional   scattered   radiation by vegetation for diffuse incidence   
   
    piLos       = (Hssh*(gap.Po(nl+1)-gap.Pso(nl+1))+ Hssu*gap.Pso(nl+1)); % directional   emitted     radiation by soil
        
    
    piLot_(i)   = piLov + piLos;
    
end
Lot_            = piLot_/pi;

%% 3. Write the output to the rad structure
[rad.Lot_,rad.Eoutte_] = deal(zeros(length(spectral.wlS),1));
rad.Lot_(IT)        = Lot_;
rad.Eoutte_(IT)     = Eoutte_;    %               emitted     diffuse radiance at top
rad.Eplut_          = Eplu_;
rad.Emint_          = Emin_;
return

% 1) CvdT, 11 December 2015.
% We subtract Emin(1), because ALL incident (thermal) radiation from Modtran
% has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!
