function [iter,rad,thermal,soil,bcu,bch,fluxes,resist_out,meteo]             ...
    = ebal(constants,options,rad,gap,  ...
    meteo,soil,canopy,leafbio,k,xyt,integr)
% function ebal.m calculates the energy balance of a vegetated surface
%
% authors:      Christiaan van der Tol (c.vandertol@utwente.nl)
%               Joris Timmermans
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
%               10 Dec 2019 (CvdT)          made a light version (layer
%                                           averaged fluxes)
% parent: SCOPE.m (script)
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
%   thermal     temperatures, aerodynamic resistances and friction velocity
%   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
%               respectively

%% 1. initialisations and other preparations for the iteration loop
% parameters for the closure loop
counter     = 0;                % iteration counter of ebal
maxit       = 100;              % maximum number of iterations
maxEBer     = 1;                % maximum energy balance error (any leaf) [Wm-2]
Wc          = 1;                % update step (1 is nominal, [0,1] possible)
CONT        = 1;                % boolean indicating whether iteration continues

% constants
MH2O        = constants.MH2O;
Mair        = constants.Mair;
rhoa        = constants.rhoa;
cp          = constants.cp;
sigmaSB     = constants.sigmaSB;

% input preparation
nl          = canopy.nlayers;   
GAM         = soil.GAM;
Ps          = gap.Ps;
kV          = canopy.kV;
xl          = canopy.xl;
LAI         = canopy.LAI;
rss         = soil.rss;

% functions for saturated vapour pressure 
es_fun      = @(T)6.107*10.^(7.5.*T./(237.3+T));
s_fun       = @(es, T) es*2.3026*7.5*237.3./(237.3+T).^2;

SoilHeatMethod = options.soil_heat_method;
if ~(options.simulation==1), SoilHeatMethod = 2; end
if SoilHeatMethod < 2
    if k > 1
        Deltat          = (datenum(xyt.t(k))-datenum(xyt.t(k-1)))*86400;           %           Duration of the time interval (s)
    else
        Deltat          = 1/48*86400;
    end   
    x 		= [1:12;1:12]'*Deltat;
    Tsold   = soil.Tsold;
end

% meteo
Ta          = meteo.Ta;
ea          = meteo.ea;
Ca          = meteo.Ca;
p           = meteo.p;
Rnuc        = rad.Rnuc;
ech         = ea*ones(nl,1);          % Leaf boundary vapour pressure (shaded/sunlit leaves)
Cch         = Ca*ones(nl,1);
ecu         = ea+0*Rnuc;
Ccu         = Ca+0*Rnuc;          % Leaf boundary CO2 (shaded/sunlit leaves)

% other preparations
e_to_q      = MH2O/Mair./p;             % Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
Fc          = Ps(1:end-1);
Fs          = [1-Ps(end),Ps(end)];      % Matrix containing values for 1-Ps and Ps of soil
%Fc          = (1-Ps(1:end-1))'/nl;      % Matrix containing values for Ps of canopy
fV          = exp(kV*xl(1:end-1));      % Vertical profile of Vcmax

% initial values for the loop
Ts          = (Ta+3)*ones(2,1);         % soil temperature (+3 for a head start of the iteration) 
Tch         = (Ta+.1)*ones(nl,1);       % leaf temperature (shaded leaves)
Tcu         = (Ta+.3)*ones(size(Rnuc)); % leaf tempeFrature (sunlit leaves)
meteo.L     = -1E6;                     % Monin-Obukhov length
[meteo_h,meteo_u]  = deal(meteo);

% this for is the exponential decline of Vcmax25. If 'lite' the dimensions
% are [nl], otherwise [13,36,nl]: 
if size(Rnuc,2)>1
    fVu      = ones(13,36,nl);
    for i = 1:nl
        fVu(:,:,i) = fV(i);
    end
else
    fVu = fV;
end

%% 2.1 Energy balance iteration loop
%Energy balance loop (Energy balance and radiative transfer)

while CONT                          % while energy balance does not close
    % 2.1. Net radiation of the components
    % Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad     = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
    Rnhc    = rad.Rnhc + rad.Rnhct;     % Canopy (shaded) net radiation
    Rnuc    = rad.Rnuc + rad.Rnuct;     % Canopy (sunlit) net radiation
    Rnhs    = rad.Rnhs+rad.Rnhst;       % Soil   (sun+sh) net radiation
    Rnus    = rad.Rnus+rad.Rnust;
    Rns     = [Rnhs Rnus]';
    
    % 2.3. Biochemical processes
    meteo_h.T       = Tch;
    meteo_h.eb      = ech;
    meteo_h.Cs      = Cch;
    meteo_h.Q       = rad.Pnh_Cab;
    meteo_u.T       = Tcu;
    meteo_u.eb      = ecu;
    meteo_u.Cs      = Ccu;
    meteo_u.Q       = rad.Pnu_Cab;
    if options.Fluorescence_model == 1
        b       = @biochemical_MD12;
    else
        b       = @biochemical;
    end
    bch     = b(leafbio,meteo_h,options,constants,fV);
    bcu     = b(leafbio,meteo_u,options,constants,fVu);
     
    % Aerodynamic roughness
    % calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]  
    [resist_out]  = resistances(constants,soil,canopy,meteo);
    meteo.ustar = resist_out.ustar;
    raa     = resist_out.raa;
    rawc    = resist_out.rawc;
    raws    = resist_out.raws;  
    rac     = (LAI+1)*(raa+rawc);
    ras     = (LAI+1)*(raa+raws);
    
    % Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    % in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    [lEch,Hch,ech,Cch,lambdah,sh]     = heatfluxes(rac,bch.rcw,Tch,ea,Ta,e_to_q,Ca,bch.Ci,constants, es_fun, s_fun);
    [lEcu,Hcu,ecu,Ccu,lambdau,su]     = heatfluxes(rac,bcu.rcw,Tcu,ea,Ta,e_to_q,Ca,bcu.Ci,constants, es_fun, s_fun);
    [lEs,Hs,~,~,lambdas,ss]           = heatfluxes(ras,rss ,Ts ,ea,Ta,e_to_q,Ca,Ca,constants, es_fun, s_fun);
    
    % integration over the layers and sunlit and shaded fractions
    Hstot   = Fs*Hs;
    %Hctot   = LAI*(meanleaf(canopy,Hcu,integr,Ps(1:end-1)) + meanleaf(canopy,Hch,'layers',1-Ps(1:end-1))  ); 
    Hctot   = aggregator(LAI,Hcu, Hch, Ps(1:end-1), canopy,integr);
    
    Htot    = Hstot + Hctot;
    if options.MoninObukhov
        meteo.L     = Monin_Obukhov(constants,meteo,Htot);     
    end
    
    % ground heat flux
    if SoilHeatMethod == 2
        G = 0.35*Rns;
        dG = 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3 * 0.35;
    elseif SoilHeatMethod == 3  % daily average flux
        G = 0*Rns;
        dG = 0*Rns;
    else
        G = GAM/sqrt(pi) * 2* sum(([Ts'; Tsold(1:end-1,:)] - Tsold)/Deltat .* (sqrt(x) - sqrt(x-Deltat)));
        G = G';
        dG = GAM/sqrt(pi) * 2* ((sqrt(x(1)) - sqrt(x(1)-Deltat)))/Deltat * ones(2,1);
    end
    
    % energy balance errors, continue criterion and iteration counter
    EBerch  = Rnhc -lEch -Hch;
    EBercu  = Rnuc -lEcu -Hcu;
    EBers   = Rns  -lEs  -Hs - G;
    
    counter     = counter+1;                   %        Number of iterations
    maxEBercu   = max(max(max(abs(EBercu))));
    maxEBerch   = max(abs(EBerch));
    maxEBers    = max(abs(EBers));
    
    CONT        = ( maxEBercu >   maxEBer    |...
        maxEBerch >   maxEBer     |...
        maxEBers  >   maxEBer)    &...
        counter   <   maxit+1;%        Continue iteration?
    if ~CONT
        if any(isnan([maxEBercu, maxEBerch, maxEBers]))
            fprintf('WARNING: NaN in fluxes, counter = %i\n', counter)
        end
        break
    end
    if counter==10, Wc = 0.8;  end
    if counter==20; Wc = 0.6;  end

    % if counter>99, plot(EBercu(:)), hold on, end
    % 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch         = Tch + Wc*EBerch./((rhoa*cp)./rac + rhoa*lambdah*e_to_q.*sh./(rac+bch.rcw)+ 4*leafbio.emis*sigmaSB*(Tch+273.15).^3);
    Tcu         = Tcu + Wc*EBercu./((rhoa*cp)./rac + rhoa*lambdau*e_to_q.*su./(rac+bcu.rcw)+ 4*leafbio.emis*sigmaSB*(Tcu+273.15).^3);
    Ts          = Ts + Wc*EBers./(rhoa*cp./ras + rhoa*lambdas*e_to_q.*ss/(ras+rss)+ 4*(1-soil.rs_thermal)*sigmaSB*(Ts+273.15).^3 + dG);
    Tch(abs(Tch)>100) = Ta;
    Tcu(abs(Tcu)>100) = Ta;
end

%% 2.2 emmissivity calculation
rad     = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
[blackleaf.tau_thermal,blackleaf.rho_thermal,blacksoil.rs_thermal] = deal(0);
rad0    = RTMt_sb(constants,rad,blacksoil,blackleaf,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);
rad.canopyemis = rad.Eoutte./rad0.Eoutte;

%% 3. Print warnings whenever the energy balance could not be solved
if counter>=maxit
    fprintf('WARNING: maximum number of iteratations exceeded\n');
    fprintf('Maximum energy balance error sunlit vegetation = %4.2f W m-2\n' ,maxEBercu);
    fprintf('Maximum energy balance error shaded vegetation = %4.2f W m-2\n' ,maxEBerch);
    fprintf('Energy balance error soil              = %4.2f W m-2\n' ,maxEBers);
    fprintf('Mean error sunlit vegetation = %4.2f W m-2\n' , mean(EBercu(:)));
end

%% 4. some more outputs
iter.counter    = counter;

thermal.Tcu     = Tcu;
thermal.Tch     = Tch;
thermal.Tsu     = Ts(2);
thermal.Tsh     = Ts(1);

fluxes.Rnctot = aggregator(LAI,Rnuc, Rnhc, Fc, canopy,integr);     % net radiation leaves
fluxes.lEctot = aggregator(LAI,lEcu, lEch, Fc, canopy,integr);     % latent heat leaves
fluxes.Hctot  = aggregator(LAI,Hcu, Hch, Fc,  canopy,integr);       % sensible heat leaves
fluxes.Actot  = aggregator(LAI,bcu.A, bch.A, Fc, canopy,integr);   % photosynthesis leaves
fluxes.Tcave  = aggregator(1,Tcu, Tch, Fc,  canopy,integr);         % mean leaf temperature
fluxes.Rnstot = Fs*Rns;           % Net radiation soil
fluxes.lEstot = Fs*lEs;           % Latent heat soil
fluxes.Hstot  = Fs*Hs;            % Sensible heat soil
fluxes.Gtot   = Fs*G;             % Soil heat flux
fluxes.Tsave  = Fs*Ts;            % Soil temperature
% fluxes.Resp   = Fs*equations.soil_respiration(Ts); %  Soil respiration = 0
fluxes.Rntot = fluxes.Rnctot + fluxes.Rnstot;
fluxes.lEtot = fluxes.lEctot + fluxes.lEstot;
fluxes.Htot = fluxes.Hctot + fluxes.Hstot;

resist_out.rss = rss; % this is simply a copy of the input rss

%% update soil temperatures history
if SoilHeatMethod < 2
    Tsold(2:end,:) = soil.Tsold(1:end-1,:);
    Tsold(1,:) 	= Ts(:);
    if isnan(Ts)
        Tsold(1,:) = Tsold(2,:);
    end
    soil.Tsold = Tsold;
end
return

function flux_tot = aggregator(LAI,sunlit_flux, shaded_flux, Fs, canopy,integr)
flux_tot = LAI*(meanleaf(canopy,sunlit_flux,integr,Fs) + meanleaf(canopy,shaded_flux,'layers',1-Fs));
return
