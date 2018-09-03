function [iter,fluxes,rad,thermal,profiles,soil]             ...  
         = ebal(iter,options,spectral,rad,gap,leafopt,  ...
                angles,meteo,soil,canopy,leafbio,xyt,k,profiles)
% function ebal.m calculates the energy balance of a vegetated surface
%
% authors:      Christiaan van der Tol (tol@itc.nl)
%               Joris Timmermans (j_timmermans@itc.nl)
% date          26 Nov 2007 (CvdT)
% updates       29 Jan 2008 (JT & CvdT)     converted into a function
%               11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
%               14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
%               31 Jul 2008 (CvdT)          Included Pntot in output
%               19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
%               07 Nov 2008 (CvdT)          Changed layout
%               18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
%                  Feb 2012 (WV)            introduced structures for variables
%                  Sep 2013 (JV, CvT)       introduced additional biochemical model
%
% parent: master.m (script)
% uses:
%       RTMt_sb.m, RTMt_planck.m (optional), RTMf.m (optional)
%       resistances.m
%       heatfluxes.m
%       biochemical.m
%       soil_respiration.m
%
% Table of contents of the function
%   
%   1. Initialisations for the iteration loop
%           intial values are attributed to variables
%   2. Energy balance iteration loop
%           iteration between thermal RTM and surface fluxes
%   3. Write warnings whenever the energy balance did not close
%   4. Calculate vertical profiles (optional)
%   5. Calculate spectrally integrated energy, water and CO2 fluxes
%
% The energy balance iteration loop works as follows:
%
% RTMo              More or less the classic SAIL model for Radiative 
%                   Transfer of sun and sky light (no emission by the vegetation)
% While continue	Here an iteration loop starts to close the energy 
%                   balance, i.e. to match the micro-meteorological model 
%                   and the radiative transfer model
% 	RTMt_sb         A numerical Radiative Transfer Model for thermal 
%                   radiation emitted by the vegetation
% 	resistances     Calculates aerodynamic and boundary layer resistances 
%                   of vegetation and soil (the micro-meteorological model)
% 	biochemical     Calculates photosynthesis, fluorescence and stomatal 
%                   resistance of leaves (or biochemical_MD12: alternative)
% 	heatfluxes      Calculates sensible and latent heat flux of soil and 
%                   vegetation
%                   Next soil heat flux is calculated, the energy balance 
%                   is evaluated, and soil and leaf temperatures adjusted 
%                   to force energy balance closure
% end {while continue}
% 
% meanleaf          Integrates the fluxes over all leaf inclinations 
%                   azimuth angles and layers, integrates over the spectrum
%
% usage:
%[iter,fluxes,rad,profiles,thermal]             ...  
%         = ebal(iter,options,spectral,rad,gap,leafopt,  ...
%                angles,meteo,soil,canopy,leafbio)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
% 
%   iter        numerical parameters used in the iteration for energy balance closure
%   options     calculation options
%   spectral    spectral resolutions and wavelengths
%   rad         incident radiation
%   gap         probabilities of direct light penetration and viewing
%   leafopt     leaf optical properties
%   angles      viewing and observation angles
%   soil        soil properties
%   canopy      canopy properties
%   leafbio     leaf biochemical parameters
%
% Output:
% 
%   iter        numerical parameters used in the iteration for energy balance closure
%   fluxes      energy balance, turbulent, and CO2 fluxes
%   rad         radiation spectra
%   profiles    vertical profiles of fluxes
%   thermal     temperatures, aerodynamic resistances and friction velocity

%% 1. initialisations and other preparations for the iteration loop
% initialisations
global constants

counter         = 0;              %           Iteration counter of ebal
maxit           = iter.maxit;
maxEBer         = iter.maxEBer;
Wc              = iter.Wc;

CONT            = 1;              %           is 0 when the calculation has finished

t               = xyt.t(k);   
Ta              = meteo.Ta;
ea              = meteo.ea;
Ca              = meteo.Ca;
Ts              = soil.Ts;
p               = meteo.p;
if options.soil_heat_method < 2 && options.simulation ==1
    if k>1
        Deltat          = (t-xyt.t(k-1))*86400;           %           Duration of the time interval (s)
    else
        Deltat          = 1/48*86400;
    end
    x 		= [1:12;1:12]'*Deltat;
    Tsold = soil.Tsold;
end

nl = canopy.nlayers;

Rnuc  = rad.Rnuc;
GAM   = soil.GAM;
Tch   = (Ta+.1)*ones(nl,1);       %           Leaf temperature (shaded leaves)
Tcu   = (Ta+.3)*ones(size(Rnuc)); %           Leaf tempeFrature (sunlit leaves)
ech   = ea*ones(nl,1);            %           Leaf H2O (shaded leaves)
ecu   = ea*ones(size(Rnuc));      %           Leaf H2O (sunlit leaves)
Cch   = Ca*ones(nl,1);            %           Leaf CO2 (shaded leaves)
Ccu   = Ca*ones(size(Rnuc));      %           Leaf CO2 (sunlit leaves)
%Tsold = Ts;                       %           Soil temperature of the previous time step
L     = -1;                       %           Monin-Obukhov length


MH2O  = constants.MH2O;
Mair  = constants.Mair;
rhoa  = constants.rhoa;
cp    = constants.cp;
g     = constants.g;
kappa = constants.kappa;
sigmaSB = constants.sigmaSB;
Ps    = gap.Ps;
nl    = canopy.nlayers;

SoilHeatMethod = options.soil_heat_method;
if ~(options.simulation==1), SoilHeatMethod = 2; end

kV   = canopy.kV;
xl   = canopy.xl;

% other preparations
e_to_q          = MH2O/Mair./p;             %           Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
Fs              = [1-Ps(end),Ps(end)];      %           Matrix containing values for 1-Ps and Ps of soil
Fc              = (1-Ps(1:end-1))'/nl;      %           Matrix containing values for Ps of canopy

if ~exist('SMCsf','var'), SMCsf = 1; end    % HERE COULD BE A STRESS FACTOR FOR VCMAX AS A FUNCTION OF SMC DEFINED
% but this is at present not
% incorporated

fVh             = exp(kV*xl(1:end-1));
fVu             = ones(13,36,60);

for i = 1:60
    fVu(:,:,i) = fVh(i);
end

LAI = canopy.LAI;

%% 2. Energy balance iteration loop

%'Energy balance loop (Energy balance and radiative transfer)

while CONT                          % while energy balance does not close
    
    % 2.1. Net radiation of the components
    % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad  = RTMt_sb(spectral,rad,soil,leafopt,canopy,gap,angles,Tcu,Tch,Ts(2),Ts(1),1);
    % Add net radiation of (1) solar and sky and (2) thermal emission model
    
    Rnhct = rad.Rnhct;
    Rnuct = rad.Rnuct;
    Rnhst = rad.Rnhst;
    Rnust = rad.Rnust;
    
    Rnhc = rad.Rnhc;
    Rnuc = rad.Rnuc;
    Rnhs = rad.Rnhs;
    Rnus = rad.Rnus;
    
    Rnch        = Rnhc + Rnhct;             %           Canopy (shaded) net radiation
    Rncu        = Rnuc + Rnuct;             %           Canopy (sunlit) net radiation
    Rnsh        = Rnhs + Rnhst;             %           Soil   (shaded) net radiation
    Rnsu        = Rnus + Rnust;             %           Soil   (sunlit) net radiation
    Rns         = [Rnsh Rnsu]';             %           Soil   (sun+sh) net radiation
    
    % 2.2. Aerodynamic roughness
    % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    
    resist_in.u   = max(meteo.u,.2);
    resist_in.L   = L;
    resist_in.LAI = canopy.LAI;
    resist_in.rbs = soil.rbs;
    resist_in.rss = soil.rss;
    resist_in.rwc = canopy.rwc;
    resist_in.zo  = canopy.zo;
    resist_in.d   = canopy.d;
    resist_in.z   = meteo.z;
    resist_in.hc  = canopy.hc;
    resist_in.w   = canopy.leafwidth;
    resist_in.Cd  = canopy.Cd;
    
    [resist_out]  = resistances(resist_in);
    
    ustar = resist_out.ustar;
    raa   = resist_out.raa;
    rawc  = resist_out.rawc;
    raws  = resist_out.raws;
    
    % 2.3. Biochemical processes
    
    % photosynthesis (A), fluorescence factor (F), and stomatal resistance (rcw), for shaded (1) and sunlit (h) leaves
    biochem_in.Fluorescence_model = options.Fluorescence_model;
    biochem_in.Type         = leafbio.Type;
    biochem_in.p            = p;
    biochem_in.m            = leafbio.m;
    biochem_in.O            = meteo.Oa;
    biochem_in.Rdparam      = leafbio.Rdparam;
    
    if options.Fluorescence_model==2    % specific for the v.Caemmerer-Magnani model
        b                   = @biochemical_MD12;
        biochem_in.Tyear        = leafbio.Tyear;
        biochem_in.beta         = leafbio.beta;
        biochem_in.qLs          = leafbio.qLs;
        biochem_in.NPQs        = leafbio.kNPQs;
        biochem_in.stressfactor = leafbio.stressfactor;
    else
        b                   = @biochemical; % specific for Berry-v.d.Tol model
        biochem_in.tempcor      = options.apply_T_corr;
        biochem_in.Tparams      = leafbio.Tparam;
        biochem_in.stressfactor = SMCsf;    
    end

    % for shaded leaves
    biochem_in.T        = Tch;
    biochem_in.eb       = ech;
    biochem_in.Vcmo     = fVh.*leafbio.Vcmo;
    biochem_in.Cs       = Cch;
    biochem_in.Q        = rad.Pnh_Cab*1E6;
    if counter == 0
        biochem_in.A        = -999;    
    else 
        biochem_in.A    = Ah;
    end
    
    biochem_out         = b(biochem_in);
    Ah                  = biochem_out.A;
    Cih                 = biochem_out.Ci;
    Fh                  = biochem_out.eta;
    rcwh                = biochem_out.rcw;
    qEh                 = biochem_out.qE; % vCaemmerer- Magnani does not generate this parameter (dummy value)
    Knh                 = biochem_out.Kn;
    
    % for sunlit leaves
    biochem_in.T        = Tcu;
    biochem_in.eb       = ecu;
    biochem_in.Vcmo     = fVu.*leafbio.Vcmo;
    biochem_in.Cs       = Ccu;
    biochem_in.Q        = rad.Pnu_Cab*1E6;
    if counter == 0
        biochem_in.A        = -999;    
    else 
        biochem_in.A    = Au;
    end
    
    biochem_out         = b(biochem_in);
 
    Au                  = biochem_out.A;
    Ciu                 = biochem_out.Ci;
    Fu                  = biochem_out.eta;
    rcwu                = biochem_out.rcw;
    qEu                 = biochem_out.qE;
    Knu                 = biochem_out.Kn;
    
    Pinh                = rad.Pnh;
    Pinu                = rad.Pnu;
    Pinh_Cab            = rad.Pnh_Cab;
    Pinu_Cab            = rad.Pnu_Cab;
    Rnh_PAR             = rad.Rnh_PAR;        
    Rnu_PAR             = rad.Rnu_PAR;
    
    % 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    
    PSIs = 0;%soil.PSIs;
    rss  = soil.rss;
    
    [lEch,Hch,ech,Cch]     = heatfluxes((LAI+1)*(raa+rawc),rcwh,Tch,ea,Ta,e_to_q,0,Ca,Cih);
    [lEcu,Hcu,ecu,Ccu]     = heatfluxes((LAI+1)*(raa+rawc),rcwu,Tcu,ea,Ta,e_to_q,0,Ca,Ciu);
    [lEs,Hs]               = heatfluxes((LAI+1)*(raa+raws),rss ,Ts ,ea,Ta,e_to_q,PSIs,Ca,Ca);
    
    % integration over the layers and sunlit and shaded fractions
    Hstot       = Fs*Hs;
    Hctot       = LAI*(Fc*Hch + meanleaf(canopy,Hcu,'angles_and_layers',Ps));
    Htot        = Hstot + Hctot;
    if SoilHeatMethod==2;
       G = 0.35*Rns;
    else      
       G = GAM/sqrt(pi) * 2* sum(([Ts'; Tsold(1:end-1,:)] - Tsold)/Deltat .* (sqrt(x) - sqrt(x-Deltat)));
       G = G';
    end
    % 2.5. Monin-Obukhov length L
    L           = -rhoa*cp*ustar.^3.*(Ta+273.15)./(kappa*g*Htot);           % [1]
    L(L<-1E3)   = -1E3;                                                     % [1] 
    L(L>1E2)    =  1E2;                                                     % [1]      
    L(isnan(L)) = -1;                                                       % [1] 
    
    % 2.6. energy balance errors, continue criterion and iteration counter
    EBerch      = Rnch -lEch -Hch;
    EBercu      = Rncu -lEcu -Hcu;
    EBers       = Rns  -lEs  -Hs - G;
    
    counter     = counter+1;                   %        Number of iterations
    maxEBercu   = max(max(max(abs(EBercu))));
    maxEBerch   = max(abs(EBerch));
    maxEBers    = max(abs(EBers));
     
    CONT        = ( maxEBercu >   maxEBer    |...
                    maxEBerch >   maxEBer     |...
                    maxEBers  >   maxEBer)    &...
                    counter   <   maxit+1;%        Continue iteration?
                
    % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1) 
    %Tch         = Ta + update(Tch-Ta,Wc,(raa + rawc)/(rhoa*cp).*(Rnch - lEch));
    Tch         = Tch + Wc*(Rnch-lEch-Hch)./((rhoa*cp)./((LAI+1)*(raa + rawc)) + 4*sigmaSB*(Tch+273.15).^3);
    %Tcu         = Ta + update(Tcu-Ta,Wc,(raa + rawc)/(rhoa*cp).*(Rncu - lEcu));
    Tcu         = Tcu + Wc*(Rncu-lEcu-Hcu)./((rhoa*cp)./((LAI+1)*(raa + rawc)) + 4*sigmaSB*(Tcu+273.15).^3);
    
    
    
    Ts(abs(Ts)>100 ) = Ta;
    %Ts          = Ta + update(Ts-Ta,Wc, (raa + raws)/(rhoa*cp).*(Rns - lEs - G));     
    Ts         = Ts + Wc*(Rns-lEs-Hs-G)./((rhoa*cp)./(raa + rawc) + (rhoa*cp)./(raa + rawc) + 4*sigmaSB*(Ts+273.15).^3);

%     if mean(abs(Hs))>1E4,
%         Ts(:) = Ta-1; Tcu(:) = Ta-1; Tch(:) = Ta-1;
%     end
    
    
    %     if t==0 || SoilHeatMethod == 2,
%         Ts      = Ta + update(Ts-Ta,Wc, (raa + raws)/(rhoa*cp).*(Rns - lEs - G));
%     else
%         Ts      = Tsold + G/GAM*sqrt(Deltat/pi);
%     end  
    % 2.8. error check (uncomment for debugging mode)
   % if max(isnan(Tch))                ,fprintf('warning: canopy temperature gives NaNs \r'),     end
   % if max(isnan(Ts))                 ,fprintf('warning soil temperature \r'),     end
  if counter>50, Wc = 0.2;  end
  
end

iter.counter = counter;
profiles.etah = Fh;
profiles.etau = Fu; 

if SoilHeatMethod<2
    Tsold(2:end,:) = soil.Tsold(1:end-1,:);
    Tsold(1,:) 	= Ts(:);
    if isnan(Ts), Tsold(1,:) = Tsold(2,:); end
    soil.Tsold = Tsold;
end

Tbr         = (rad.Eoutte/constants.sigmaSB)^0.25;
Lot_        = Planck(spectral.wlS',Tbr);
rad.LotBB_  = Lot_;           % Note that this is the blackbody radiance!

%% 3. Print warnings whenever the energy balance could not be solved
if counter>=maxit,
    fprintf(1,'%s \n','warning: maximum number of iteratations exceeded');
    fprintf(1,'%s ',['Energy balance error sunlit vegetation = ',sprintf('%4.2f',maxEBercu),'W m-2 ']);
    fprintf(1,'%s ',['Energy balance error shaded vegetation = ',sprintf('%4.2f',maxEBerch),'W m-2 ']);
    fprintf(1,'%s ',['Energy balance error soil              = ',sprintf('%4.2f',maxEBers ),'W m-2 ']);
    fprintf(1,'\r');
end

%% 4. Calculate the output per layer
if options.calc_vert_profiles   
    [Hcu1d  ]           = meanleaf(canopy,Hcu,          'angles');   % [nli,nlo,nl]      mean sens heat sunlit leaves
    [lEcu1d ]           = meanleaf(canopy,lEcu,         'angles');   % [nli,nlo,nl]      mean latent sunlit leaves
    [Au1d   ]           = meanleaf(canopy,Au,           'angles');   % [nli,nlo,nl]      mean phots sunlit leaves
    [Fu_Pn1d]           = meanleaf(canopy,Fu.*Pinu_Cab, 'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
    [qEuL   ]           = meanleaf(canopy,qEu,          'angles');   % [nli,nlo,nl]      mean fluor sunlit leaves
    %[Pnu1d  ]           = meanleaf(canopy,Pinu,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    %[Pnu1d_Cab  ]       = meanleaf(canopy,Pinu_Cab,     'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    [Rnu1d  ]           = meanleaf(canopy,Rncu,         'angles');   % [nli,nlo,nl]      mean net PAR sunlit leaves
    [Tcu1d  ]           = meanleaf(canopy,Tcu,          'angles');   % [nli,nlo,nl]      mean temp sunlit leaves
    
    profiles.Tchave     = mean(Tch);                                           % [1]               mean temp shaded leaves
    profiles.Tch        = Tch;                                                 % [nl]
    profiles.Tcu1d      = Tcu1d;                                               % [nl]
    profiles.Tc1d       = (1-Ps(1:nl)).*Tch       + Ps(1:nl).*(Tcu1d);         % [nl]              mean temp leaves, per layer
    profiles.Hc1d       = (1-Ps(1:nl)).*Hch       + Ps(1:nl).*(Hcu1d);         % [nl]              mean sens heat leaves, per layer
    profiles.lEc1d      = (1-Ps(1:nl)).*lEch      + Ps(1:nl).*(lEcu1d);        % [nl]              mean latent heat leaves, per layer
    profiles.A1d        = (1-Ps(1:nl)).*Ah        + Ps(1:nl).*(Au1d);          % [nl]              mean photos leaves, per layer
    profiles.F_Pn1d     = ((1-Ps(1:nl)).*Fh.*Pinh_Cab + Ps(1:nl).*(Fu_Pn1d));  %[nl]           mean fluor leaves, per layer
    profiles.qE         = ((1-Ps(1:nl)).*qEh      + Ps(1:nl).*(qEuL));         %[nl]           mean fluor leaves, per layer
    %profiles.Pn1d       = ((1-Ps(1:nl)).*Pinh     + Ps(1:nl).*(Pnu1d));        %[nl]           mean photos leaves, per layer
    %profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pinh_Cab + Ps(1:nl).*(Pnu1d_Cab));        %[nl]           mean photos leaves, per layer
    profiles.Rn1d       = ((1-Ps(1:nl)).*Rnch     + Ps(1:nl).*(Rnu1d));        %[nl]
end
        

%% 5. Calculate spectrally integrated energy, water and CO2 fluxes
% sum of all leaves, and average leaf temperature 
%     (note that averaging temperature is physically not correct...)

Rnctot          = LAI*(Fc*Rnch + meanleaf(canopy,Rncu,'angles_and_layers',Ps)); % net radiation leaves
lEctot          = LAI*(Fc*lEch + meanleaf(canopy,lEcu,'angles_and_layers',Ps)); % latent heat leaves
Hctot           = LAI*(Fc*Hch  + meanleaf(canopy,Hcu ,'angles_and_layers',Ps)); % sensible heat leaves
Actot           = LAI*(Fc*Ah   + meanleaf(canopy,Au  ,'angles_and_layers',Ps)); % photosynthesis leaves
Tcave           =     (Fc*Tch  + meanleaf(canopy,Tcu ,'angles_and_layers',Ps)); % mean leaf temperature
Pntot           = LAI*(Fc*Pinh + meanleaf(canopy,Pinu,'angles_and_layers',Ps)); % net PAR leaves
Pntot_Cab       = LAI*(Fc*Pinh_Cab + meanleaf(canopy,Pinu_Cab,'angles_and_layers',Ps)); % net PAR leaves
Rntot_PAR       = LAI*(Fc*Rnh_PAR  + meanleaf(canopy,Rnu_PAR, 'angles_and_layers',Ps));% net PAR leaves
aPAR_Cab_eta        = LAI*(Fc*(profiles.etah .* Rnh_PAR) + meanleaf(canopy,(profiles.etau .* Rnu_PAR), 'angles_and_layers',Ps)); 
% ... green ePAR * relative fluorescence emission efficiency


% sum of soil fluxes and average temperature
%   (note that averaging temperature is physically not correct...)
Rnstot          = Fs*Rns;           %                   Net radiation soil
lEstot          = Fs*lEs;           %                   Latent heat soil
%Hstot          = Fs*Hs;            %                   Sensible heat soil
Gtot            = Fs*G;             %                   Soil heat flux
Tsave           = Fs*Ts;            %                   Soil temperature
Resp            = Fs*soil_respiration(Ts);%             Soil respiration

% total fluxes (except sensible heat), all leaves and soil
Atot            = Actot;            %                   GPP
Rntot           = Rnctot + Rnstot;  %                   Net radiation
lEtot           = lEctot + lEstot;  %                   Latent heat
%Htot           = Hctot  + Hstot;   %                   Sensible heat

fluxes.Rntot    = Rntot;  % [W m-2]             total net radiation (canopy + soil)
fluxes.lEtot    = lEtot;  % [W m-2]             total latent heat flux (canopy + soil)
fluxes.Htot     = Htot;   % [W m-2]             total sensible heat flux (canopy + soil)
fluxes.Atot     = Atot;   % [umol m-2 s-1]      total net CO2 uptake (canopy + soil)
fluxes.Rnctot   = Rnctot; % [W m-2]             canopy net radiation
fluxes.lEctot   = lEctot; % [W m-2]             canopy latent heat flux
fluxes.Hctot    = Hctot;  % [W m-2]             canopy sensible heat flux
fluxes.Actot    = Actot;  % [umol m-2 s-1]      canopy net CO2 uptake
fluxes.Rnstot   = Rnstot; % [W m-2]             soil net radiation
fluxes.lEstot   = lEstot; % [W m-2]             soil latent heat flux
fluxes.Hstot    = Hstot;  % [W m-2]             soil sensible heat flux
fluxes.Gtot     = Gtot;   % [W m-2]             soil heat flux
fluxes.Resp     = Resp;   % [umol m-2 s-1]      soil respiration
fluxes.aPAR     = Pntot;  % [umol m-2 s-1]      absorbed PAR
fluxes.aPAR_Cab = Pntot_Cab;% [umol m-2 s-1]      absorbed PAR
fluxes.aPAR_Wm2 = Rntot_PAR;% [W m-2]      absorbed PAR
fluxes.aPAR_Cab_eta = aPAR_Cab_eta;

thermal.Ta    = Ta;       % [oC]                air temperature (as in input)
thermal.Ts    = Ts;       % [oC]                soil temperature, sunlit and shaded [2x1]
thermal.Tcave = Tcave;    % [oC]                weighted average canopy temperature
thermal.Tsave = Tsave;    % [oC]                weighted average soil temperature
thermal.raa   = raa;      % [s m-1]             total aerodynamic resistance above canopy
thermal.rawc  = rawc;     % [s m-1]             aerodynamic resistance below canopy for canopy
thermal.raws  = raws;     % [s m-1]             aerodynamic resistance below canopy for soil 
thermal.ustar = ustar;    % [m s-1]             friction velocity
thermal.Tcu   = Tcu;
thermal.Tch   = Tch;

fluxes.Au     = Au;
fluxes.Ah     = Ah;

profiles.Knu     = Knu;
profiles.Knh     = Knh;
% function Tnew = update(Told, Wc, innovation)
%     Tnew        = Wc.*innovation + (1-Wc).*Told;
% return