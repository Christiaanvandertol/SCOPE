function rad = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Tsu,Tsh,obsdir)
% function 'RTMt_planck' calculates the spectrum of outgoing thermal
% radiation in hemispherical and viewing direction
%
% Authors: Wout Verhoef and Christiaan van der Tol (tol@itc.nl)
% Date: 5 November 2007
% Update:   14 Nov 2007
%           16 Nov 2007 CvdT    improved calculation of net radiation
%           17 Dec 2007 JT      simplified, removed net radiation
%           07 Nov 2008 CvdT    changed layout
%           16 Mar 2009 CvdT    removed calculation of Tbright
%           12 Apr 2013 CvdT    introduced structures
%
% Table of contents of the function:
%   0.      preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     parameters
%       0.3     geometric factors of Observer
%       0.4     geometric factors associated with extinction and scattering
%       0.5     geometric factors to be used later with rho and tau
%       0.6     fo for all leaf angle/azumith classes     
%   1       calculation of upward and downward fluxes
%   2       outgoing fluxes, hemispherical and in viewing direction
%   A1      function planck (external function is now used)
%
% Usage:
%   function rad = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Tsu,Tsh,obsdir)
%
% Input:
%   Symbol      Description                 Unit            Dimension
%   ------      -----------                 ----            ---------
%   Tcu         temperature sunlit leaves   C               [13,36,nl]
%   Tch         temperature shaded leaves   C               [nl]
%   Tsu         temperature sunlit soil     C               [1]
%   Tsu         temperature shaded soil     C               [1]
%   rad         a structure containing 
%   soil        a structure containing soil reflectance
%   canopy      a structure containing LAI and leaf inclination

%   Ps          probability of sunlit leaves                [nl+1]
%   Po          probability of viewing a leaf or soil       [nl+1]
%   Pso         probability of viewing a sunlit leaf/soil   [nl+1]
%   K           extinction coefficient in viewing dir       [1]
%   tto         viewing angle               (degrees)       [1]
%   psi         azimuth angle difference between solar and viewing position
%
% Output
%   Symbol      Description                             Unit        Dimension
%   ------      -----------                             ----        ---------
%   Loutt_      Spectrum of outgoing hemispherical rad  (W m-2 um-1)[nwl]
%   Lot_        Spectrum of outgoing rad in viewing dir (W m-2 um-1)[nwl]
%   Eplu        Total downward diffuse radiation        (W m-2)     [nl+1]
%   Emin        Total downward diffuse radiation        (W m-2)     [nl+1]
%
% Notes:
%   nl      number of layers
%   nwl     number of wavelengths of input (net PAR)
% '_'means: a flux at different wavelengths (a vertically oriented vector)

%% 0.0 globals
global constants

%% 0.1 parameters

%for speed-up the calculation only uses thermal part of the spectrum
IT          = spectral.IwlT;   % 
wlt         = spectral.wlT;

deg2rad     = constants.deg2rad;
nl          = canopy.nlayers;
lidf        = canopy.lidf;
litab       = canopy.litab;
lazitab     = canopy.lazitab;
nlazi       = length(lazitab);
tto         = angles.tto;
psi         = angles.psi;
Ps          = gap.Ps;
Po          = gap.Po;
Pso         = gap.Pso;
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

%% 0.2  initialization of output variables
Hcsui               = zeros(nl,1);      %             [nl]      
piLot_              = zeros(1,length(IT)); %          [1,nwlt]    
[Emin_,Eplu_]       = deal(zeros(nl+1,length(IT)));       % [nl+1,nwlt]

%% 0.3 geometric factors of Observer
if obsdir
    cos_tto             = cos(tto*deg2rad);%            [1]         cos observation angle
    sin_tto             = sin(tto*deg2rad);%            [1]         sin observation angle
end

%% 0.4 geometric factors associated with extinction and scattering
cos_ttl             = cos(litab*deg2rad);           %   [nli]       cos leaf inclination angles
if obsdir
    sin_ttl         = sin(litab*deg2rad);           %   [nli]       sin leaf inclination angles
    cos_ttlo        = cos((lazitab-psi)*deg2rad);     % [nlazi]     sin leaf orientation angles
end
bfli                = cos_ttl.^2;                  %    [nli]
bf                  = bfli'*lidf;                   %   [1]

%% 0.5 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
ddb                 = 0.5*(1+bf);                   %   [1]         W of dif2dif back    (f1^2 + f2^2)                
ddf                 = 0.5*(1-bf);                   %   [1]         W of dif2dif forward (2*f1*f2    )
if obsdir
    dob             = 0.5*(K+bf);                   %   [1]         W of dif2dir back    (fo*f1      )
    dof             = 0.5*(K-bf);                   %   [1]         W of dif2dir back    (fo*f1      ) 
end

%% 0.6 fo for all leaf angle/azumith classes
if obsdir
    Co              = cos_ttl*cos_tto;             %    [nli]       pag 305 modified by Joris
    So              = sin_ttl*sin_tto;             %    [nli]       pag 305 modified by Joris
    
    cos_deltao      = Co*ones(1,nlazi) + So*cos_ttlo; % [nli,nlazi]	projection of leaves in  in direction of sun %(pag 125/126) 
    fo              = cos_deltao/abs(cos_tto);      %   [nli,nlazi]
end

%% 1. calculation of upward and downward fluxes
sigb                = ddb*rho + ddf*tau;            %   [nwlt]      diffuse     backscatter scattering coefficient for diffuse  incidence  pag 305
sigf                = ddf*rho + ddb*tau;            %   [nwlt]      diffuse     forward     scattering coefficient for forward  incidence 
if obsdir
    vb              = dob*rho + dof*tau;            %   [nwlt]      directional backscatter scattering coefficient for diffuse  incidence 
    vf              = dof*rho + dob*tau;            %   [nwlt]      directional forward     scattering coefficient for diffuse  incidence 
end
a                   = 1-sigf;                       %   [nwlt]      attenuation
m                   = sqrt(a.^2-sigb.^2);           %   [nwlt]	
rinf                = (a-m)./sigb;                  %   [nwlt]      reflection coefficient for infinite thick canopy
rinf2               = rinf.*rinf;                   %   [nwlt]        

fHs                 = (1-rinf2).*(1-rs)./(1-rinf.*rs);
fHc                 = iLAI*m.*(1-rinf);
fbottom             = (rs-rinf)./(1-rinf.*rs);

for i = 1:length(IT);
    % 1.1 radiance by components
    Hcsui3          = pi*Planck(wlt(i),Tcu+273.15,epsc(i));
    Hcshi           = pi*Planck(wlt(i),Tch+273.15,epsc(i));
    Hssui           = pi*Planck(wlt(i),Tsu+273.15,epss(i));
    Hsshi           = pi*Planck(wlt(i),Tsh+273.15,epss(i));
    % 1.2 radiance by leaf layers Hc and by soil Hs
    for j = 1:nl
        Hcsui2      = Hcsui3(:,:,j);                 %  [nli,nlazi] sunlit  leaves
        Hcsui(j)    = mean(Hcsui2'*lidf);            %  [nl]        sunlit  vegetation radiance per layer
    end    
    Hci             = Hcsui.*Ps(1:nl) + Hcshi.*(1-Ps(1:nl));  %[nl] emitted vegetation radiance per layer
    Hsi             = Hssui.*Ps(nl+1) + Hsshi.*(1-Ps(nl+1));  %[1]  emitted soil       radiance per layer
    
    % 1.3 Diffuse radiation
    cont            = 1;                            %   [1]         continue iteration (1:yes, 0:no)
    counter         = 0;                            %   [1]         iteration counter
    F1              = zeros(nl+1,1);                %   [nl+1]
    F2              = zeros(nl+1,1);                %   [nl+1]  
    F1top           = 0;                            %   [1]
    while cont
        F1topn      = -rinf(i)*F2(1);
        F1(1)       = F1topn;
        for j       = 1:nl
            F1(j+1) = F1(j)*(1-m(i)*iLAI)+ fHc(i)*Hci(j);
        end        
        F2(nl+1)    = fbottom(i)*F1(nl+1) + fHs(i)*Hsi;          
        for j       = nl:-1:1
            F2(j)   = F2(j+1)*(1-m(i)*iLAI) + fHc(i)*Hci(j);
        end      
        cont        = abs(F1topn-F1top)>crit;     %     [1]         check to continue
        F1top       = F1topn;                     %     [1]         Reset F1topn
        counter     = counter + 1;                %     [1]
    end
    Emini           = (F1+rinf(i)*F2)/(1-rinf2(i)); %   [nl+1]
    Eplui           = (F2+rinf(i)*F1)/(1-rinf2(i)); %   [nl+1]
    
    Emin_(:,i)      = Emini;                        %               downwelling diffuse radiance per layer
    Eplu_(:,i)      = Eplui;                        %               upwelling   diffuse radiance  
   
    % 1.4 Directional radiation
    if obsdir
        for j = 1:nl        
            Hcsui2      = Hcsui3(:,:,j).*abs(fo);
            Hcsui(j)    = mean(Hcsui2'*lidf);            
        end

        piLo1           = iLAI*K*sum(Hcshi.*(Po(1:nl)-Pso(1:nl)));         % directional   emitted     radiation by shaded leaves        
        piLo2           = iLAI*sum(Hcsui.*(Pso(1:nl         )));         % directional   emitted     radiation by sunlit leaves
        piLo3           = iLAI*((vb(i)*Emini(1:nl) + vf(i)*Eplui(1:nl))'*Po(1:nl));% directional   scattered   radiation by vegetation for diffuse incidence 
        piLo4           = Hsshi*(Po(nl+1)-Pso(nl+1));
        piLo5           = Hssui* Pso(nl+1);                                % directional   emitted     radiation by sunlit/shaded Soil
        piLo6           = rs(i)*Emini(nl+1)*Po(nl+1);                              % directional   scattered   radiation by soil       for diffuse incidence
        
        piLot_(i)       = piLo1 + sum(piLo2) + piLo3 + piLo4 + piLo5 + piLo6;      % directional   total       radiation 
        Lot_            = piLot_/pi;
    end
end

%% 2. Write the output to structure rad
[rad.Lot_,rad.Eoutte_] = deal(zeros(length(spectral.wlS),1));
rad.Lot_(IT)        = Lot_;
rad.Eoutte_(IT)      = Eplu_(1,:);                     %               emitted     diffuse radiance at top
rad.Eplut_          = Eplu_;
rad.Emint_          = Emin_;
return
