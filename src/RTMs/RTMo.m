function [rad,gap,canopy,profiles] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,constants,meteo,options)
% calculates the spectra of hemisperical and directional observed visible 
% and thermal radiation (fluxes E and radiances L), as well as the single 
% and bi-directional gap probabilities
%
% the function does not require any non-standard Matlab functions. No
% changes to the code have to be made to operate the function for a
% particular canopy. All necessary parameters and variables are input or
% global and need to be specified elsewhere.
%
% Authors:      Wout Verhoef            (w.verhoef@utwente.nl) 
%               Christiaan van der Tol  (c.vandertol@utwente.nl)
%               Joris Timmermans        ()
%
% updates:      10 Sep 2007 (CvdT)      - calculation of Rn
%                5 Nov 2007             - included observation direction
%               12 Nov 2007             - included abs. PAR spectrum output
%                                       - improved calculation efficiency
%               13 Nov 2007             - written readme lines
%               11 Feb 2008 (WV&JT)     - changed Volscat
%                           (JT)        - small change in calculation Po,Ps,Pso
%                                        - introduced parameter 'lazitab'
%                                       - changed nomenclature
%                                       - Appendix IV: cosine rule
%               04 Aug 2008 (JT)        - Corrections for Hotspot effect in the probabilities
%               05 Nov 2008 (CvdT)      - Changed layout
%               04 Jan 2011 (JT & CvdT) - Included Pso function (Appendix IV)
%                                       - removed the analytical function (for checking)  
%               02 Oct 2012 (CvdT)      - included incident PAR in output
%
%               Jan/Feb 2013 (WV)       - Major revision towards SCOPE version 1.40:
%                                       - Parameters passed using structures
%                                       - Improved interface with MODTRAN atmospheric data
%                                       - Now also calculates 4-stream
%                                         reflectances rso, rdo, rsd and rdd
%                                         analytically 
%               Apri 2013 (CvT)         - improvements in variable names
%                                           and descriptions
%                  Dec 2019 CvdT        mSCOPE representation, lite option
%
%                                      
% Table of contents of the function
%
%   0.      Preparations
%       0.1     parameters
%       0.2     initialisations
%   1.      Geometric quantities
%       1.1     general geometric quantities
%       1.2     geometric factors associated with extinction and scattering
%       1.3     geometric factors to be used later with rho and tau
%       1.4     solar irradiance factor for all leaf orientations
%       1.5     probabilities Ps, Po, Pso
%   2.      Calculation of upward and downward fluxes
%   3.      Outgoing fluxes, hemispherical and in viewing direction, spectrum
%   4.      Net fluxes, spectral and total, and incoming fluxes
%   A1      functions J1 and J2 (introduced for stable solutions)
%   A2      function volscat
%   A3      function e2phot
%   A4      function Pso
%
% references:
%{1} Verhoef (1998), 'Theory of radiative transfer models applied in
%    optical remote sensing of vegetation canopies'. PhD Thesis Univ. Wageninegn
%{2} Verhoef, W., Jia, L., Xiao, Q. and Su, Z. (2007) Unified optical -
%    thermal four - stream radiative transfer theory for homogeneous 
%    vegetation canopies. IEEE Transactions on geoscience and remote 
%    sensing, 45,6.
%{3} Verhoef (1985), 'Earth Observation Modeling based on Layer Scattering
%    Matrices', Remote sensing of Environment, 17:167-175 
%              
% Usage:
% function [rad,gap,profiles] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options)  
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   atmo        MODTRAN atmospheric parameters
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   angles      viewing and observation angles
%   meteo       has the meteorological variables. Is only used to correct
%               the total irradiance if a specific value is provided
%               instead of the usual Modtran output.
%   rad         initialization of the structure of the output 'rad'
%   options     simulation options. Here, the option
%               'calc_vert_profiles' is used, a boolean that tells whether 
%               or not to output data of 60 layers separately.
%
% Output:
%   gap         probabilities of direct light penetration and viewing
%   rad         a large number of radiative fluxes: spectrally distributed 
%               and integrated, and canopy radiative transfer coefficients.
%   profiles    vertical profiles of radiation variables such as absorbed
%               PAR.


%% 0. Preparations
deg2rad = constants.deg2rad;    % degree to rad
wl      = spectral.wlS;         % SCOPE wavelengths as a column-vector
nwl     = length(wl);           %
wlPAR   = spectral.wlPAR;       % PAR wavelength range
minPAR  = min(wlPAR);           % min PAR
maxPAR  = max(wlPAR);           % max PAR
Ipar    = find(wl>=minPAR & wl<=maxPAR); % Indices for PAR wavelenghts within wl
tts     = angles.tts;           % solar zenith angle
tto     = angles.tto;           % observer zenith angle
psi     = angles.psi;           % relative azimuth anglee

nl      = canopy.nlayers;       % number of canopy layers (nl)
litab   = canopy.litab;         % SAIL leaf inclibation angles % leaf inclination angles PY
lazitab = canopy.lazitab;       % leaf azimuth angles relative to the sun
nlazi   = canopy.nlazi;         % number of azimuth angles (36)
LAIeff  = canopy.LAI;           % leaf area index
lidf    = canopy.lidf;          % leaf inclination distribution function
xl      = canopy.xl;            % all levels except for the top
Cv      = canopy.Cv;
crownd  = canopy.crowndiameter;
theta   = crownd / (canopy.hc/2);   % crown shape factor
dx      = 1/nl;
LAI     = LAIeff/Cv;

rho = leafopt.refl;
tau = leafopt.tran;
kChlrel = leafopt.kChlrel;
rs      =   soil.refl;          % [nwl,nsoils] soil reflectance spectra
epsc    =   1-rho-tau;          % [nl,nwl]        emissivity of leaves
epss    =   1-rs;              % [nwl]        emissivity of soil
iLAI    =   LAI/nl;             % [1]          LAI of elementary layer

% initializations
Rndif                       = zeros(nl,1);                 % [nl+1]         abs. diffuse rad soil+veg
[Pdif,Pndif,Pndif_Cab,Rndif_Cab,Rndif_PAR]      = deal(zeros(nl,1));           % [nl]           incident and net PAR veg
%[Es_,Emin_,Eplu_]           = deal(zeros(nl+1,nwl));       % [nl+1,nwl]     direct, up and down diff. rad.
[Rndif_]                    = zeros(nl,nwl);               % [nl,nwl]       abs diff and PAR veg.
[Pndif_,Pndif_Cab_,Rndif_PAR_,Rndif_Cab_]         = deal(zeros(nl,length(Ipar)));
%[Puc,Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR]    = deal(zeros(nli,nlazi,nl));   %

%% 1. Geometric quantities
% 1.1 general geometric quantities. these variables are scalars
cos_tts     = cos(tts*deg2rad);             %           cos solar       angle
tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle
sin_tts     = sin(tts*deg2rad);             %           sin solar       angle
tan_tts     = tan(tts*deg2rad);             %           tan observation angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));

% 1.2 geometric factors associated with extinction and scattering
[chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab);   % volume scattering

cos_ttlo    = cos(lazitab*deg2rad);         % [1,36]    cos leaf azimuth angles

cos_ttli    = cos(litab*deg2rad);           % [13]      cos leaf angles
sin_ttli    = sin(litab*deg2rad);           % [13]      sinus leaf angles

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle
koli        = chi_o./cos_tto;               % [13]      p307{1} extinction coefficient in direction of observer   per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto);    % [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto);    % [13]      pag 309{1}
bfli        = cos_ttli.^2;                  % [13]

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.
K           = koli'*lidf;                   %           pag 307{1}    extinction coefficient in direction of observer
bf          = bfli'*lidf;                   %
sob         = sobli'*lidf;                  %           weight of specular2directional back    scatter coefficient
sof         = sofli'*lidf;                  %           weight of specular2directional forward scatter coefficient
% 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
% these variables are scalars
sdb         = 0.5*(k+bf);                   % fs*f1
sdf         = 0.5*(k-bf);                   % fs*f2     weight of specular2diffuse     foward  scatter coefficient
ddb         = 0.5*(1+bf);                   % f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient
ddf         = 0.5*(1-bf);                   % 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient
dob         = 0.5*(K+bf);                   % fo*f1     weight of diffuse2directional  back    scatter coefficient
dof         = 0.5*(K-bf);                   % fo*f2     weight of diffuse2directional  forward scatter coefficient

% 1.4 solar irradiance factor for all leaf orientations
Css          = cos_ttli*cos_tts;             % [nli]     pag 305 modified by Joris
Ss          = sin_ttli*sin_tts;             % [nli]     pag 305 modified by Joris
cos_deltas  = Css*ones(1,nlazi) + Ss*cos_ttlo;  % [nli,nlazi]
fs          = abs(cos_deltas/cos_tts);         % [nli,nlazi] pag 305

%% 2. Calculation of reflectance
% 2.1  reflectance, transmittance factors in a thin layer
% the following are vectors with lenght [nl,nwl]
sigb        = ddb*rho + ddf*tau;            % [nl,nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence
sigf        = ddf*rho + ddb*tau;            % [nl,nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence
sb          = sdb*rho + sdf*tau;            % [nl,nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence
sf          = sdf*rho + sdb*tau;            % [nl,nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence
vb          = dob*rho + dof*tau;            % [nl,nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence
vf          = dof*rho + dob*tau;            % [nl,nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence
w           = sob*rho + sof*tau;            % [nl,nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)
a           = 1-sigf;                       % [nl,nwl]     attenuation

% crown area projections
Cs          = 1-(1-Cv).^(1./cos_tts);       %           crown cover fraction projected in the direction of the sun
Co          = 1-(1-Cv).^(1./cos_tto);       %           crown cover fraction projected in the direction of viewing

%% 3. Flux calculation

% diffuse fluxes within the vegetation covered part
%iLAI    =   LAI/nl; 
tau_ss = repmat(1-k.*iLAI,nl,1);  %REPLACE when LIDF profile ready.
tau_dd = (1-a.*iLAI);
tau_sd = sf.*iLAI;
rho_sd = sb.*iLAI;
rho_dd = sigb.*iLAI;
[R_sd,R_dd,Xss,Xsd,Xdd] = calc_reflectances(tau_ss,tau_sd,tau_dd,rho_dd,rho_sd,rs,nl,nwl);
rdd     = R_dd(1,:)'* Cv + (1-Cv)*rs;
rsd     = R_sd(1,:)'* Cs + (1-Cs)*rs;
[Esun_,Esky_] = calcTOCirr(atmo,meteo,rdd,rsd,wl,nwl);

[Emins_,Eplus_] = calc_fluxprofile(Esun_,0*Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl);
[Emind_,Eplud_] = calc_fluxprofile(0*Esun_,Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl);
Emin_ = Emins_+Emind_;
Eplu_ = Eplus_+Eplud_;
%%
% 1.5 probabilities Ps, Po, Pso
d_h         = crownd / canopy.hc;
Ps0          =   exp(k*xl*LAI) ;                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po0          =   exp(K*xl*LAI) ;
Ps0(1:nl)    =   Ps0(1:nl) *(1-exp(-k*LAI*dx))/(k*LAI*dx);                                      % Correct Ps/Po for finite dx
Po0(1:nl)    =   Po0(1:nl) *(1-exp(-K*LAI*dx))/(K*LAI*dx);  % Correct Ps/Po for finite dx

xls         =   min(1,d_h*(1-Cv)/Cv / tan_tts);% theta is the ratio of d/h
xlo         =   min(1,d_h*(1-Cv)/Cv / tan_tto);

C1s         = exp(-k*iLAI); % extinction
C1o         = exp(-K*iLAI);
C2s         = tan_tts/theta/nl * exp(-k*iLAI/4); % extra term
C2o         = tan_tto/theta/nl * exp(-K*iLAI/4);
% coefficients in the extra term: 2 because it is a triangle
% 4 (need to check carefully) because we account for the extinction within the thin layer. The
% average path out is integral of all positions in the triangle to the exit
% on the left side.
%

[Ps1,Po1]     = deal(ones(nl+1,1));
Ps1(1)      = C2s;
Po1(1)      = C2o;
for j = 2:nl+1
    Ps1(j)   = C1s*Ps1(j-1)+C2s *double(((j-1)/nl)<=xls);
    Po1(j)   = C1o*Po1(j-1)+C2o *double(((j-1)/nl)<=xlo);
    
end


Ps          = Ps0 + Ps1; % for sure Ps0 has to be multiplied by a weight!
Po          = Po0 + Po1;
% [nl+1]  p154{1} probability of viewing a leaf in observation dir


q           =   canopy.hot;
Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
end
Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error
gap.Pso      = Pso;

%%
% 3.3 outgoing fluxes, hemispherical and in viewing direction, spectrum
% in viewing direction, spectral due to diffuse light

% vegetation contribution
piLocd_     = (sum(vb.*Po(1:nl).*Emind_(1:nl,:)) +...
              sum(vf.*Po(1:nl).*Eplud_(1:nl,:)))'*iLAI;
% soil contribution          
piLosd_     = rs.*(Emind_(end,:)'*Po(end)); 

% in viewing direction, spectral due to direct solar light
% vegetation contribution
piLocu_     = (sum(vb.*Po(1:nl).*Emins_(1:nl,:)) +...
              sum(vf.*Po(1:nl).*Eplus_(1:nl,:))+...
              sum(w.*Pso(1:nl).*Esun_'))'*iLAI;
% soil contribution
piLosu_     = rs.* (Emins_(end,:)'*Po(end) + Esun_*Pso(end)); 

piLod_      = piLocd_ + piLosd_;        % [nwl] piRad in obsdir from Esky
piLou_      = piLocu_ + piLosu_;        % [nwl] piRad in obsdir from Eskun
piLoc_      = piLocu_ + piLocd_;        % [nwl] piRad in obsdir from vegetation
piLos_      = piLosu_ + piLosd_;        % [nwl] piRad in obsdir from soil

piLo_       = piLoc_ + piLos_;          % [nwl] piRad in obsdir

Lo_         = piLo_/pi;                 % [nwl] Rad in obsdir
rso         = piLou_./Esun_;            % [nwl] obsdir reflectance of solar beam
rdo         = piLod_./Esky_;            % [nlw] obsir reflectance of sky irradiance
Refl        = piLo_./(Esky_+Esun_);     % [nwl] rso and rdo are not computed 

%% 4. net fluxes, spectral and total, and incoming fluxes
%4.1 incident PAR at the top of canopy, spectral and spectrally integrated
P_          = e2phot(wl(Ipar)*1E-9,(Esun_(Ipar)+Esky_(Ipar)),constants);      %
P           = .001 * Sint(P_,wlPAR);                                % mol m-2s-1
%Psun        = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar),constants),wlPAR);   % Incident solar PAR in PAR units
% Incident and absorbed solar radiation

%4.2 Absorbed radiation
%    absorbed radiation in Wm-2         (Asun)
%    absorbed PAR in mol m-2s-1         (Pnsun)
%    absorbed PAR in Wm-2               (Rnsun_PAR)
%    absorbed PAR by Chl in mol m-2s-1  (Pnsun_Cab)

[Asun,Pnsun,Rnsun_PAR,Pnsun_Cab,Rnsun_Cab]= deal(zeros(nl,1));
for j=1:nl
    Asun(j)        = 0.001 * Sint(Esun_.*epsc(j,:)',wl);                                 % Total absorbed solar radiation
    Pnsun(j)       = 0.001 * Sint(e2phot(wlPAR*1E-9,Esun_(Ipar).*epsc(j,Ipar)',constants),wlPAR);  % Absorbed solar radiation in PAR range in moles m-2 s-1
    Rnsun_Cab(j)   = 0.001 * Sint(Esun_(Ipar)'.*epsc(j,Ipar).*kChlrel(j,Ipar),wlPAR);
    Rnsun_PAR(j)   = 0.001 * Sint(Esun_(Ipar)'.*epsc(j,Ipar),wlPAR);
    Pnsun_Cab(j)   = 0.001 * Sint(e2phot(wlPAR*1E-9,kChlrel(j,Ipar)'.*Esun_(Ipar).*epsc(j,Ipar)',constants),wlPAR);
end

%4.3 total direct radiation (incident and net) per leaf area (W m-2 leaf)

% total direct radiation (incident and net) per leaf area (W m-2 leaf)
%Pdir        = fs * Psun;                        % [13 x 36]   incident
if options.lite
    fs      = lidf'*mean(fs,2);%
    Rndir       = fs * Asun(j);                        % [13 x 36 x nl]   net
    Pndir       = fs * Pnsun(j);                       % [13 x 36 x nl]   net PAR
    Pndir_Cab   = fs * Pnsun_Cab(j);                   % [13 x 36 x nl]   net PAR Cab
    Rndir_Cab   = fs * Rnsun_Cab(j);                   %   net PAR energy units
    Rndir_PAR   = fs * Rnsun_PAR(j);                   % [13 x 36 x nl]
else
    [Rndir,Pndir,Pndir_Cab,Rndir_Cab,Rndir_PAR]= deal(zeros(13,36,nl));
    for j=1:nl
        Rndir(:,:,j)       = fs * Asun(j);                        % [13 x 36 x nl]   net
        Pndir(:,:,j)       = fs * Pnsun(j);                       % [13 x 36 x nl]   net PAR
        Pndir_Cab(:,:,j)   = fs * Pnsun_Cab(j);                   % [13 x 36 x nl]   net PAR Cab
        Rndir_Cab(:,:,j)   = fs * Rnsun_Cab(j);                   %   net PAR energy units
        Rndir_PAR(:,:,j)   = fs * Rnsun_PAR(j);                   % [13 x 36 x nl]   net PAR energy units
    end
end

%4.4 total diffuse radiation (net) per leaf area (W m-2 leaf)
for j = 1:nl     % 1 top nl is bottom
    % diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
    E_         = .5*(Emin_(j,:) + Emin_(j+1,:)+ Eplu_(j,:)+ Eplu_(j+1,:));
   
    % incident PAR flux, integrated over all wavelengths (moles m-2 s-1)
    Pdif(j)    = .001 * Sint(e2phot(wlPAR'*1E-9,E_(Ipar),constants),wlPAR);  % [nl] , including conversion mW >> W
  
    % net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
    Rndif_(j,:)         = E_.*epsc(j,:);                                                    % [nl,nwl]  Net (absorbed) radiation by leaves
    Pndif_(j,:)         = .001 *(e2phot(wlPAR*1E-9, Rndif_(j,Ipar)',constants))';                     % [nl,nwl]  Net (absorbed) as PAR photons
    Rndif_Cab_(j,:)     = .001 *(kChlrel(j,Ipar).*Rndif_(j,Ipar))';    % [nl,nwl]  Net (absorbed) as PAR photons by Cab
    Pndif_Cab_(j,:)     = .001 *(e2phot(wlPAR*1E-9, (kChlrel(j,Ipar).*Rndif_(j,Ipar))',constants))';    % [nl,nwl]  Net (absorbed) as PAR photons by Cab
    Rndif_PAR_(j,:)     = Rndif_(j,Ipar);                                                   % [nl,nwlPAR]  Net (absorbed) as PAR energy 
   
    % net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
    Rndif(j)            = .001 * Sint(Rndif_(j,:),wl);              % [nl]  Full spectrum net diffuse flux
    Pndif(j)            =        Sint(Pndif_(j,Ipar),wlPAR);        % [nl]  Absorbed PAR
    Pndif_Cab(j)        =        Sint(Pndif_Cab_(j,Ipar),wlPAR);    % [nl]  Absorbed PAR by Cab integrated
    Rndif_Cab(j)        = .001 * Sint(Rndif_Cab_(j,Ipar),wlPAR);      % [nl]  Absorbed PAR by Cab integrated
    Rndif_PAR(j)        = .001 * Sint(Rndif_PAR_(j,Ipar),wlPAR);    % [nl]  Absorbed PAR by Cab integrated
end

% soil layer, direct and diffuse radiation
Rndirsoil   = .001 * Sint(Esun_.*epss,wl);          % [1] Absorbed solar flux by the soil
Rndifsoil   = .001 * Sint(Emin_(nl+1,:).*epss',wl); % [1] Absorbed diffuse downward flux by the soil (W m-2)

% net (n) radiation R and net PAR P per component: sunlit (u), shaded (h) soil(s) and canopy (c),
% [W m-2 leaf or soil surface um-1]
Rnhc        = Rndif;            % [nl] shaded leaves or needles
Pnhc        = Pndif;            % [nl] shaded leaves or needles
Pnhc_Cab    = Pndif_Cab;        % [nl] shaded leaves or needles
Rnhc_Cab    = Rndif_Cab;        % [nl] shaded leaves or needles
Rnhc_PAR    = Rndif_PAR;        % [nl] shaded leaves or needles

if ~options.lite
    [Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR,Rnuc_Cab] = deal(0*Rndir);
    for j = 1:nl
        %Puc(:,:,j)  = Pdir(:,:,j)      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc(:,:,j) = Rndir(:,:,j)     + Rndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc(:,:,j)  = Pndir(:,:,j)     + Pndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
        Pnuc_Cab(:,:,j)  = Pndir_Cab(:,:,j)  + Pndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_PAR(:,:,j)  = Rndir_PAR(:,:,j)  + Rndif_PAR(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
        Rnuc_Cab(:,:,j)  = Rndir_Cab(:,:,j)  + Rndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles   
    end
else   
    %    Puc(:,:,j)  = Pdir      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc = Rndir     + Rndif;     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc  = Pndir     + Pndif;     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc_Cab  = Pndir_Cab  + Pndif_Cab;% [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc_Cab  = Rndir_Cab  + Rndif_Cab;% [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc_PAR  = Rndir_PAR  + Rndif_PAR;% [13,36,nl] Total fluxes on sunlit leaves or needles       
end
Rnus        = Rndifsoil + Rndirsoil; % [1] sunlit soil 
Rnhs        = Rndifsoil;  % [1] shaded soil

if options.calc_vert_profiles   
    [Pnu1d  ]           = meanleaf(canopy,Pnuc,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    [Pnu1d_Cab  ]       = meanleaf(canopy,Pnuc_Cab,     'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    
    profiles.Pn1d       = ((1-Ps(1:nl)).*Pnhc     + Ps(1:nl).*(Pnu1d));        %[nl]           mean photos leaves, per layer
    profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pnhc_Cab + Ps(1:nl).*(Pnu1d_Cab));        %[nl]           mean photos leaves, per layer
else
    profiles = struct;
end

%% 5 Model output
% up and down and hemispherical out, cumulative over wavelenght
Eout_       = Eplu_(1,:)'*Cs + (1-Cs)*rs.*(Esun_+Esky_);
Eouto       = 0.001 * Sint(Eout_(spectral.IwlP),spectral.wlP);  %     [1] hemispherical out, in optical range (W m-2)
Eoutt       = 0.001 * Sint(Eout_(spectral.IwlT),spectral.wlT);  %     [1] hemispherical out, in thermal range (W m-2)
Lot         = 0.001 * Sint(Lo_(spectral.IwlT),spectral.wlT);    %     [1] hemispherical out, in thermal range (W m-2)

% place output in structure rad
gap.k       = k;        % extinction cofficient in the solar direction
gap.K       = K;        % extinction cofficient in the viewing direction
gap.Ps      = Ps;       % gap fraction in the solar direction
gap.Po      = Po;       % gap fraction in the viewing direction
rad.rsd     = rsd;      % TOC directional-hemispherical reflectance
rad.rdd     = rdd;      % TOC hemispherical-hemispherical reflectance
rad.rdo     = rdo;      % TOC hemispherical-directional reflectance 
rad.rso     = rso;      % TOC directional-directional reflectance   
rad.refl    = Refl;     % TOC reflectance
rad.rho_dd  = rho_dd;   % diffuse-diffuse reflectance for the thin layers
rad.tau_dd  = tau_dd;   % diffuse-diffuse transmittance for the thin layers
rad.rho_sd  = rho_sd;   % direct-diffuse reflectance for the thin layers
rad.tau_ss  = tau_ss;   % direct-direct transmittance for the thin layers
rad.tau_sd  = tau_sd;   % direct-diffuse transmittance for the thin layers
rad.R_sd    = R_sd;
rad.R_dd    = R_dd;
rad.vb      = vb;       % directional backscatter coefficient for diffuse incidence
rad.vf      = vf;       % directional forward scatter coefficient for diffuse incidence
rad.sigf    = sigf;     % forward scatter coefficient for specular flux
rad.sigb    = sigb;     % backscatter coefficient for specular flux
rad.Esun_   = Esun_;    % [nwlx1 double]   incident solar spectrum (mW m-2 um-1)
rad.Esky_   = Esky_;    % [nwlx1 double]   incident sky spectrum (mW m-2 um-1)
rad.PAR     = P;        % [1 double]        incident spectrally integrated PAR (moles m-2 s-1)
rad.Eplu_   = Eplu_;    % [nlxnwl double]  upward diffuse radiation in the canopy (mW m-2 um-1)
rad.Emin_   = Emin_;    % [nlxnwl double]  downward diffuse radiation in the canopy (mW m-2 um-1)
rad.Emins_  = Emins_;   % [nlxnwl double]  downward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
rad.Emind_  = Emind_;   % [nlxnwl double]  downward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)
rad.Eplus_  = Eplus_;   % [nlxnwl double]  upward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
rad.Eplud_  = Eplud_;   % [nlxnwl double]  upward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)
rad.Lo_     = Lo_;      % [nwlx1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)
rad.Eout_   = Eout_;    % [nwlx1 double]   TOC upward radiation (mW m-2 um-1)
rad.Eouto   = Eouto;    % [1 double]        TOC spectrally integrated upward optical ratiation (W m-2)
rad.Eoutt   = Eoutt;    % [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)
rad.Lot     = Lot;
rad.Rnhs    = Rnhs;     % [1 double]        net radiation (W m-2) of shaded soil
rad.Rnus    = Rnus;     % [1 double]        net radiation (W m-2) of sunlit soil
rad.Rnhc    = Rnhc;     % [60x1 double]     net radiation (W m-2) of shaded leaves
rad.Rnuc    = Rnuc;     % [13x36x60 double] net radiation (W m-2) of sunlit leaves
rad.Pnh     = 1E6*Pnhc;     % [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
rad.Pnu     = 1E6*Pnuc;     % [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
rad.Pnh_Cab = 1E6*Pnhc_Cab;% [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Pnu_Cab = 1E6*Pnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
rad.Rnh_Cab = Rnhc_Cab; % [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_Cab = Rnuc_Cab; % [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
rad.Rnh_PAR = Rnhc_PAR;     % [60x1 double]     net PAR absorbed by Cab (W m-2) of shaded leaves
rad.Rnu_PAR = Rnuc_PAR;     % [13x36x60 double] net PAR absorbed (W m-2) of sunlit
rad.Xdd     =   Xdd;
rad.Xsd     = Xsd;
rad.Xss     = Xss;

gap.LAI_Cv = LAI;

% Rn = canopy.LAI*(meanleaf(canopy,rad.Rnhc,'layers',(1-Ps(1:nl)))+meanleaf(canopy,rad.Rnuc,'angles_and_layers',Ps(1:nl)))
% %y1 = canopy.Cv*(rad.Eplu_(end,:)-rad.Eplu_(1,:) + rad.Emin_(1,:) - rad.Emin_(end,:));
% y1 = canopy.Cv*(rad.Eplu_(end,:)-rad.Eplu_(1,:) + rad.Emin_(1,:) - rad.Emin_(end-1,:));
% y2 = canopy.Cv*(1-gap.Ps(end))*rad.Esun_';
% 
% %y2 = (canopy.Cs-gap.Ps(end))*rad.Esun_';
% y = y1+ y2;
% Sint(y,spectral.wlS)*1E-3
% iLAI
%% APPENDIX I function volscat

function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
nli     = length(ttli);

psi_rad         = psi*deg2rad*ones(nli,1);

cos_psi         = cos(psi*deg2rad);                 %   cosine of relative azimuth angle

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

cos_tto         = cos(tto*deg2rad);                 %   cosine of observer zenith angle
sin_tto         = sin(tto*deg2rad);                 %   sine   of observer zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}

Co              = cos_ttli*cos_tto;                 %   p305{1}
So              = sin_ttli*sin_tto;                 %   p305{1}

As              = max([Ss,Cs],[],2);
Ao              = max([So,Co],[],2);

bts             = acos(-Cs./As);                    %   p305{1}
bto             = acos(-Co./Ao);                    %   p305{2}

chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);

delta1          = abs(bts-bto);                     %   p308{1}
delta2          = pi-abs(bts + bto - pi);           %   p308{1}

Tot             = psi_rad + delta1 + delta2;        %   pag 130{1}

bt1             = min([psi_rad,delta1],[],2);
bt3             = max([psi_rad,delta2],[],2);
bt2             = Tot - bt1 - bt3;

T1              = 2.*Cs.*Co + Ss.*So.*cos_psi;
T2              = sin(bt2).*(2*As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

Jmin            = (   bt2).*T1 - T2;
Jplus           = (pi-bt2).*T1 + T2;

frho            =  Jplus/(2*pi^2);
ftau            = -Jmin /(2*pi^2);

% pag.309 wl-> pag 135{1}
frho            = max([zeros(nli,1),frho],[],2);
ftau            = max([zeros(nli,1),ftau],[],2);
return

%% APPENDIX II function e2phot

function molphotons = e2phot(lambda,E,constants)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

e           = ephoton(lambda,constants);
photons     = E./e;
molphotons  = photons./constants.A;
return;

function E = ephoton(lambda,constants)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of
%wavelength lambda (m)

h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;

%% APPENDIX III function Pso

function pso    =   Psofunction(K,k,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
    
end
return;

function [R_sd,R_dd,Xss,Xsd,Xdd]  = calc_reflectances(tau_ss,tau_sd,tau_dd,rho_dd,rho_sd,rs,nl,nwl)
[R_sd,R_dd] =   deal(zeros(nl+1,nwl));      % surface reflectance
[Xsd,Xdd]   =   deal(zeros(nl,nwl));        % Effective transmittance
Xss         =   zeros(1,nl);                % Effective transmittance
R_sd(nl+1,:)   =   rs;
R_dd(nl+1,:)   =   rs;
for j=nl:-1:1 % from bottom to top. note nl+1 the background. 1 is the top of canopy.
    Xss      = tau_ss(j);
    dnorm       = 1-rho_dd(j,:).*R_dd(j+1,:);
    Xsd(j,:)    = (tau_sd(j,:)+tau_ss(j).*R_sd(j+1,:).*rho_dd(j,:))./dnorm;
    Xdd(j,:)    = tau_dd(j,:)./dnorm;
    R_sd(j,:)   = rho_sd(j,:)+tau_dd(j,:).*(R_sd(j+1,:).*Xss+R_dd(j+1,:).*Xsd(j,:));
    R_dd(j,:)   = rho_dd(j,:)+tau_dd(j,:).*R_dd(j+1,:).*Xdd(j,:);
end
return;

function [Emin_,Eplu_,Es_] = calc_fluxprofile(Esun_,Esky_,rs,Xss,Xsd,Xdd,R_sd,R_dd,nl,nwl)
[Es_,Emin_,Eplu_]           = deal(zeros(nl+1,nwl));       % [nl+1,nwl]     direct, up and down diff. rad.
Es_(1,:)       = Esun_;
Emin_(1,:)     = Esky_;

for j=1:nl % from top to bottom
    Es_(j+1,:)    =   Xss.*Es_(j,:);
    Emin_(j+1,:)  =   Xsd(j,:).*Es_(j,:)+Xdd(j,:).*Emin_(j,:);
    Eplu_(j,:)    =   R_sd(j,:).*Es_(j,:)+R_dd(j,:).*Emin_(j,:);
end
Eplu_(j+1,:)    =   rs'.*(Es_(j+1,:)+Emin_(j+1,:)); 
%CvdT added calculation of Eplu_(j+1,:)
return;


function [Esun_,Esky_] = calcTOCirr(atmo,meteo,rdd,rsd,wl,nwl)
% calculation of incoming light Esun_ and Esky_
% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths

Fd      = zeros(nwl,1);
Ls      = Planck(wl,meteo.Ta+273.15);

if ~isfield(atmo,'Esun_')
    
    t1  = atmo.M(:,1);
    t3  = atmo.M(:,2);
    t4  = atmo.M(:,3);
    t5  = atmo.M(:,4);
    t12 = atmo.M(:,5);
    t16 = atmo.M(:,6);
    
    % radiation fluxes, downward and upward (these all have dimenstion [nwl]
    % first calculate hemispherical reflectances rsd and rdd according to SAIL
    % these are assumed for the reflectance of the surroundings
    % rdo is computed with SAIL as well
    % assume Fd of surroundings = 0 for the momemnt
    % initial guess of temperature of surroundings from Ta;

    Esun_   = max(1E-6,pi*t1.*t4);
    Esky_   = max(1E-6,pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16));
    
    % fractional contributions of Esun and Esky to total incident radiation in
    % optical and thermal parts of the spectrum
    if meteo.Rin ~= -999
        % fractional contributions of Esun and Esky to total incident radiation in
        % optical and thermal parts of the spectrum
        
        [fEsuno,fEskyo,fEsunt,fEskyt]          = deal(0*Esun_);   %initialization
        
        J_o             = wl<3000;                          %find optical spectrum
        Esunto          = 0.001 * Sint(Esun_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration), including conversion mW >> W
        Eskyto          = 0.001 * Sint(Esky_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration)
        Etoto           = Esunto + Eskyto;                  %Calculate total fluxes
        fEsuno(J_o)     = Esun_(J_o)/Etoto;                 %fraction of contribution of Sun fluxes to total light
        fEskyo(J_o)     = Esky_(J_o)/Etoto;                 %fraction of contribution of Sky fluxes to total light
        
        J_t             = wl>=3000;                         %find thermal spectrum
        Esuntt          = 0.001 * Sint(Esun_(J_t),wl(J_t)); %Themal solar fluxes
        Eskytt          = 0.001 * Sint(Esky_(J_t),wl(J_t)); %Thermal Sky fluxes
        Etott           = Eskytt + Esuntt;                  %Total
        fEsunt(J_t)     = Esun_(J_t)/Etott;                 %fraction from Esun
        fEskyt(J_t)     = Esky_(J_t)/Etott;                 %fraction from Esky
        
        Esun_(J_o) = fEsuno(J_o)*meteo.Rin;
        Esky_(J_o) = fEskyo(J_o)*meteo.Rin;
        Esun_(J_t) = fEsunt(J_t)*meteo.Rli;
        Esky_(J_t) = fEskyt(J_t)*meteo.Rli;
    end
    
else
    Esun_ = atmo.Esun_;
    Esky_ = atmo.Esky_;
end
return;
