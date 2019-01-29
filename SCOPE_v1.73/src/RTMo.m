function [rad,gap,profiles] = RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options)    
      
%% function RTMo 
%
% calculates the spectra of hemisperical and directional observed visible 
% and thermal radiation (fluxes E and radiances L), as well as the single 
% and bi-directional gap probabilities
%
% the function does not require any non-standard Matlab functions. No
% changes to the code have to be made to operate the function for a
% particular canopy. All necessary parameters and variables are input or
% global and need to be specified elsewhere.
%
% Authors:      Wout Verhoef            (verhoef@nlr.nl) 
%               Christiaan van der Tol  (tol@itc.nl)
%               Joris Timmermans        (j_timmermans@itc.nl)
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
deg2rad = pi/180;
wl      = spectral.wlS';                % SCOPE wavelengths as a column-vector
nwl     = length(wl);

wlP     = spectral.wlP;
wlT     = spectral.wlT;
wlPAR   = spectral.wlPAR';               % PAR wavelength range
minPAR  = min(wlPAR);
maxPAR  = max(wlPAR);
Ipar    = find(wl>=minPAR & wl<=maxPAR); % Indices for PAR wavelenghts within wl

tts     = angles.tts;              % solar zenith angle
tto     = angles.tto;              % observer zenith angle
psi     = angles.psi;              % relative azimuth anglee

nl      = canopy.nlayers;       % number of canopy layers (60)
litab   = canopy.litab;         % SAIL leaf inclibation angles
lazitab = canopy.lazitab;       % leaf azimuth angles relative to the sun
nli     = canopy.nlincl;        % numler of leaf inclinations (13)
nlazi   = canopy.nlazi;         % number of azimuth angles (36)
LAI     = canopy.LAI;           % leaf area index
lidf    = canopy.lidf;          % leaf inclination distribution function
x       = canopy.x;             % all levels except for the top
dx      = 1/nl;

kChlrel = leafopt.kChlrel;
rho     = leafopt.refl;         % [nwl]        leaf/needle reflection  
tau     = leafopt.tran;         % [nwl]        leaf/needle transmission
rs      = soil.refl;            % [nwl,nsoils] soil reflectance spectra
epsc    = 1-rho-tau;            % [nwl]        emissivity of leaves
epss    = 1-rs;                 % [nwl]        emissivity of soil
iLAI    = LAI/nl;               % [1]          LAI of elementary layer
xl      = [0; x];               % [nl+1]       all levels + soil

% 0.2 initialisations (allocation of memory)
Rndif                       = zeros(nl,1);                  % [nl+1]         abs. diffuse rad soil+veg
[Pdif,Pndif,Pndif_Cab,Rndif_PAR] = deal(zeros(nl,1));       % [nl]           incident and net PAR veg
[Emin_,Eplu_]               = deal(zeros(nl+1,nwl));        % [nl+1,nwl]     up and down diff. rad.
[Rndif_]                    = zeros(nl,nwl);                % [nl,nwl]       abs diff and PAR veg.
[Pndif_,Pndif_Cab_,Rndif_PAR_]       = deal(zeros(nl,length(Ipar)));
[Puc,Rnuc,Pnuc,Pnuc_Cab,Rnuc_PAR]    = deal(zeros(nli,nlazi,nl));   % [nli,nlazi,nl] inc and net rad and PAR sunlit

%% 1.0 Geometric quantities
% 1.1 general geometric quantities
% these variables are scalars
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
Cs          = cos_ttli*cos_tts;             % [nli]     pag 305 modified by Joris
Ss          = sin_ttli*sin_tts;             % [nli]     pag 305 modified by Joris

cos_deltas  = Cs*ones(1,nlazi) + Ss*cos_ttlo;  % [nli,nlazi]
fs          = abs(cos_deltas/cos_tts);         % [nli,nlazi] pag 305

% 1.5 probabilities Ps, Po, Pso
Ps          =   exp(k*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po          =   exp(K*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in observation dir

Ps(1:nl)    =   Ps(1:nl) *(1-exp(-k*LAI*dx))/(k*LAI*dx);                                      % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-K*LAI*dx))/(K*LAI*dx);  % Correct Ps/Po for finite dx


q           =   canopy.hot;
Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
end

Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error

gap.Pso      = Pso;

%% 2. Calculation of upward and downward fluxes

% the following are vectors with lenght nwl
sigb        = ddb*rho + ddf*tau;            % [nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence 
sigf        = ddf*rho + ddb*tau;            % [nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence 
sb          = sdb*rho + sdf*tau;            % [nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence 
sf          = sdf*rho + sdb*tau;            % [nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence 
vb          = dob*rho + dof*tau;            % [nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence 
vf          = dof*rho + dob*tau;            % [nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence 
w           = sob*rho + sof*tau;            % [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)         
a           = 1-sigf;                       % [nwl]     attenuation
m           = sqrt(a.^2-sigb.^2);           % [nwl]
rinf        = (a-m)./sigb;                  % [nwl]
rinf2       = rinf.*rinf;                   % [nwl]

% direct solar radiation
J1k        = calcJ1(-1, m,k,LAI);          % [nwl]
J2k        = calcJ2( 0, m,k,LAI);          % [nwl]
J1K        = calcJ1(-1, m,K,LAI);          % [nwl]   % added for calculation of rdo
J2K        = calcJ2( 0, m,K,LAI);          % [nwl]   % added for calculation of rdo

e1          = exp(-m*LAI);                  % [nwl]
e2          = e1.^2;                        % [nwl]
re          = rinf.*e1;                     % [nwl]

denom       = 1-rinf2.*e2;                  % [nwl]

s1          = sf+rinf.*sb;
s2          = sf.*rinf+sb;
v1          = vf+rinf.*vb;
v2          = vf.*rinf+vb;

Pss         = s1.*J1k;          % [nwl]
Qss         = s2.*J2k;          % [nwl]

Poo         = v1.*J1K;          % (nwl)   % added for calculation of rdo
Qoo         = v2.*J2K;          % [nwl]   % added for calculation of rdo

tau_ss      = exp(-k*LAI);                  % [1]
tau_oo      = exp(-K*LAI);                  % [1]

Z           = (1 - tau_ss * tau_oo)/(K + k);  % needed for analytic rso

tau_dd      = (1-rinf2).*e1 ./denom;        % [nwl]
rho_dd      = rinf.*(1-e2)  ./denom;        % [nwl]
tau_sd      = (Pss-re.*Qss) ./denom;        % [nwl]
tau_do      = (Poo-re.*Qoo) ./denom;        % [nwl]
rho_sd      = (Qss-re.*Pss) ./denom;        % [nwl]
rho_do      = (Qoo-re.*Poo) ./denom;        % (nwl)

T1          = v2.*s1.*(Z-J1k*tau_oo)./(K+m)+v1.*s2.*(Z-J1K*tau_ss)./(k+m);
T2          = -(Qoo.*rho_sd+Poo.*tau_sd).*rinf;
rho_sod     = (T1+T2)./(1-rinf2);

rho_sos     = w * sum(Pso(1:nl))*iLAI;
rho_so      = rho_sod + rho_sos;

Pso2w       = Pso(nl+1);

% Analytical rso following SAIL

rso         = rho_so + rs * Pso2w                                        ...
              + ((tau_sd+tau_ss*rs.*rho_dd)*tau_oo+(tau_sd+tau_ss).*tau_do) ...
               .*rs./denom;
           
% Extract MODTRAN atmosphere parameters at the SCOPE wavelengths
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

denom       = 1-rs.*rho_dd;

% SAIL analytical reflectances

rsd     = rho_sd + (tau_ss + tau_sd).*rs.*tau_dd./denom;
rdd     = rho_dd + tau_dd.*rs.*tau_dd./denom;

rdo     = rho_do + (tau_oo + tau_do).*rs.*tau_dd./denom;



% assume Fd of surroundings = 0 for the momemnt
% initial guess of temperature of surroundings from Ta;

Fd      = zeros(nwl,1);
Ls      = equations.Planck(wl,atmo.Ta+273.15);

% Solar and sky irradiance using 6 atmosperic functions
%keyboard
Esun_   = pi*t1.*t4;
Esky_   = pi./(1-t3.*rdd).*(t1.*(t5+t12.*rsd)+Fd+(1-rdd).*Ls.*t3+t16);

% fractional contributions of Esun and Esky to total incident radiation in
% optical and thermal parts of the spectrum
[fEsuno,fEskyo,fEsunt,fEskyt]          = deal(0*Esun_);   %initialization

J_o             = wl<3000;                          %find optical spectrum
Esunto          = 0.001 * helpers.Sint(Esun_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration), including conversion mW >> W 
Eskyto          = 0.001 * helpers.Sint(Esky_(J_o),wl(J_o)); %Calculate optical sun fluxes (by Integration)
Etoto           = Esunto + Eskyto;                  %Calculate total fluxes
fEsuno(J_o)     = Esun_(J_o)/Etoto;                 %fraction of contribution of Sun fluxes to total light
fEskyo(J_o)     = Esky_(J_o)/Etoto;                 %fraction of contribution of Sky fluxes to total light

J_t             = wl>=3000;                         %find thermal spectrum
Esuntt          = 0.001 * helpers.Sint(Esun_(J_t),wl(J_t)); %Themal solar fluxes
Eskytt          = 0.001 * helpers.Sint(Esky_(J_t),wl(J_t)); %Thermal Sky fluxes
Etott           = Eskytt + Esuntt;                  %Total
fEsunt(J_t)     = Esun_(J_t)/Etott;                 %fraction from Esun 
fEskyt(J_t)     = Esky_(J_t)/Etott;                 %fraction from Esky

if meteo.Rin ~= -999
    Esun_(J_o) = fEsuno(J_o)*meteo.Rin;
    Esky_(J_o) = fEskyo(J_o)*meteo.Rin;
    Esun_(J_t) = fEsunt(J_t)*meteo.Rli;
    Esky_(J_t) = fEskyt(J_t)*meteo.Rli;
end

Eplu_1      = rs.*((tau_ss+tau_sd).*Esun_+tau_dd.*Esky_)./denom;
Eplu0       = rho_sd.*Esun_ + rho_dd.*Esky_ + tau_dd.*Eplu_1;
Emin_1      = tau_sd.*Esun_ + tau_dd.*Esky_ + rho_dd.*Eplu_1;
delta1      = Esky_  - rinf.*Eplu0;
delta2      = Eplu_1 - rinf.*Emin_1;

% calculation of the fluxes in the canopy
for i = 1:nwl
    J1kx            = calcJ1(xl,m(i),k,LAI);    %           [nl]
    J2kx            = calcJ2(xl,m(i),k,LAI);    %           [nl]
    F1              = Esun_(i)*J1kx*(sf(i)+rinf(i)*sb(i)) + delta1(i)*exp( m(i)*LAI*xl);      %[nl]
    F2              = Esun_(i)*J2kx*(sb(i)+rinf(i)*sf(i)) + delta2(i)*exp(-m(i)*LAI*(xl+1));  %[nl]
    Emin_(:,i)      = (F1+rinf(i)*F2)/(1-rinf2(i));%        [nl,nwl]
    Eplu_(:,i)      = (F2+rinf(i)*F1)/(1-rinf2(i));%        [nl,nwl]
end

% Incident and absorbed solar radiation
Psun        = 0.001 * helpers.Sint(e2phot(wlPAR*1E-9,Esun_(Ipar)),wlPAR);   % Incident solar PAR in PAR units
%Psky        = 0.001 * helpers.Sint(e2phot(wlPAR*1E-9,Esky_(Ipar)),wlPAR);
Asun        = 0.001 * helpers.Sint(Esun_.*epsc,wl);                         % Total absorbed solar radiation
Pnsun       = 0.001 * helpers.Sint(e2phot(wlPAR*1E-9,Esun_(Ipar).*epsc(Ipar)),wlPAR);  % Absorbed solar radiation  in PAR range in moles m-2 s-1
Rnsun_PAR   = 0.001 * helpers.Sint(Esun_(Ipar).*epsc(Ipar),wlPAR);  
Pnsun_Cab   = 0.001 * helpers.Sint(e2phot(wlPAR*1E-9,kChlrel(Ipar).*Esun_(Ipar).*epsc(Ipar)),wlPAR);  
                                                                    % Absorbed solar radiation by Cab in PAR range in moles m-2 s-1

%% 3. outgoing fluxes, hemispherical and in viewing direction, spectrum
% (CvdT 071105: compared with analytical solution: is OK)
% hemispherical, spectral
Eout_   = Eplu_(1,:)';                  %           [nwl]

% in viewing direction, spectral
piLoc_      = (vb.*(Emin_(1:nl,:)'*Po(1:nl)) +...
               vf.*(Eplu_(1:nl,:)'*Po(1:nl)) +...
               w.*Esun_*sum(Pso(1:nl)))*iLAI;
piLos_      = rs.*Emin_(nl+1,:)'*Po(nl+1) + rs.*Esun_*Pso(nl+1);
piLo_       = piLoc_ + piLos_;              %           [nwl]
Lo_         = piLo_/pi;

% up and down and hemispherical out, cumulative over wavelenght
IwlP        = spectral.IwlP;
IwlT        = spectral.IwlT;
Eouto       = 0.001 * helpers.Sint(Eout_(IwlP),wlP);   %     [1] hemispherical out, in optical range (W m-2)
Eoutt       = 0.001 * helpers.Sint(Eout_(IwlT),wlT);   %     [1] hemispherical out, in thermal range (W m-2)

%% 4. net fluxes, spectral and total, and incoming fluxes
% incident PAR at the top of canopy, spectral and spectrally integrated
P_          = e2phot(wl(Ipar)*1E-9,(Esun_(Ipar)+Esky_(Ipar)));  
P           = .001 * helpers.Sint(P_,wlPAR);

% total direct radiation (incident and net) per leaf area (W m-2 leaf)
Pdir        = fs * Psun;                        % [13 x 36]   incident    
Rndir       = fs * Asun;                        % [13 x 36]   net         
Pndir       = fs * Pnsun;                       % [13 x 36]   net PAR     
Pndir_Cab   = fs * Pnsun_Cab;                   % [13 x 36]   net PAR Cab 
Rndir_PAR   = fs * Rnsun_PAR;                   % [13 x 36]   net PAR energy units 


% canopy layers, diffuse radiation
for j = 1:nl                                    
    % diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
    E_         = .5*(Emin_(j,:) + Emin_(j+1,:)+ Eplu_(j,:)+ Eplu_(j+1,:));
    
    % incident PAR flux, integrated over all wavelengths (moles m-2 s-1) 
    Pdif(j)    = .001 * helpers.Sint(e2phot(wlPAR*1E-9,E_(Ipar)'),wlPAR);  % [nl] , including conversion mW >> W 
    
    % net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
    Rndif_(j,:)         = E_.*epsc';                                                   % [nl,nwl]  Net (absorbed) radiation by leaves
    Pndif_(j,:)         = .001 *(e2phot(wlPAR*1E-9, Rndif_(j,Ipar)'))';                % [nl,nwl]  Net (absorbed) as PAR photons
    Pndif_Cab_(j,:)     = .001 *(e2phot(wlPAR*1E-9, kChlrel(Ipar).*Rndif_(j,Ipar)'))';  % [nl,nwl]  Net (absorbed) as PAR photons by Cab
    Rndif_PAR_(j,:)     = Rndif_(j,Ipar);  % [nl,nwlPAR]  Net (absorbed) as PAR energy
    
    % net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
    Rndif(j)            = .001 * helpers.Sint(Rndif_(j,:),wl);              % [nl]  Full spectrum net diffuse flux
    Pndif(j)            =        helpers.Sint(Pndif_(j,Ipar),wlPAR);        % [nl]  Absorbed PAR 
    Pndif_Cab(j)        =        helpers.Sint(Pndif_Cab_(j,Ipar),wlPAR);    % [nl]  Absorbed PAR by Cab integrated
    Rndif_PAR(j)        = .001 * helpers.Sint(Rndif_PAR_(j,Ipar),wlPAR);    % [nl]  Absorbed PAR by Cab integrated
end

% soil layer, direct and diffuse radiation
Rndirsoil   = .001 * helpers.Sint(Esun_.*epss,wl);          % [1] Absorbed solar flux by the soil
Rndifsoil   = .001 * helpers.Sint(Emin_(nl+1,:).*epss',wl); % [1] Absorbed diffuse downward flux by the soil (W m-2)

% net (n) radiation R and net PAR P per component: sunlit (u), shaded (h) soil(s) and canopy (c),
% [W m-2 leaf or soil surface um-1]
Rnhc        = Rndif;            % [nl] shaded leaves or needles      
Pnhc        = Pndif;            % [nl] shaded leaves or needles 
Pnhc_Cab    = Pndif_Cab;        % [nl] shaded leaves or needles 
Rnhc_PAR    = Rndif_PAR;        % [nl] shaded leaves or needles 

for j = 1:nl
    Puc(:,:,j)  = Pdir      + Pdif(j);      % [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc(:,:,j) = Rndir     + Rndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc(:,:,j)  = Pndir     + Pndif(j);     % [13,36,nl] Total fluxes on sunlit leaves or needles
    Pnuc_Cab(:,:,j)  = Pndir_Cab  + Pndif_Cab(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
    Rnuc_PAR(:,:,j)  = Rndir_PAR  + Rndif_PAR(j);% [13,36,nl] Total fluxes on sunlit leaves or needles
end
Rnhs        = Rndifsoil;                  % [1] shaded soil
Rnus        = Rndifsoil + Rndirsoil;      % [1] sunlit soil

%% 
if options.calc_vert_profiles   
    [Pnu1d  ]           = equations.meanleaf(canopy,Pnuc,         'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    [Pnu1d_Cab  ]       = equations.meanleaf(canopy,Pnuc_Cab,     'angles');   % [nli,nlo,nl]      mean net radiation sunlit leaves
    
    profiles.Pn1d       = ((1-Ps(1:nl)).*Pnhc     + Ps(1:nl).*(Pnu1d));        %[nl]           mean photos leaves, per layer
    profiles.Pn1d_Cab   = ((1-Ps(1:nl)).*Pnhc_Cab + Ps(1:nl).*(Pnu1d_Cab));        %[nl]           mean photos leaves, per layer
else
    profiles = struct;
end

%% place output in structure rad
gap.k       = k;
gap.K       = K;
gap.Ps      = Ps;
gap.Po      = Po;

rad.rsd     = rsd;
rad.rdd     = rdd;
rad.rdo     = rdo;
rad.rso     = rso;

rad.vb      = vb;
rad.vf      = vf;

rad.Esun_   = Esun_;        % [2162x1 double]   incident solar spectrum (mW m-2 um-1)
rad.Esky_   = Esky_;        % [2162x1 double]   incident sky spectrum (mW m-2 um-1)
rad.PAR     = P;            % [1 double]        incident spectrally integrated PAR (moles m-2 s-1)

rad.fEsuno  = fEsuno;       % [2162x1 double]   normalized spectrum of direct light (optical)
rad.fEskyo  = fEskyo;       % [2162x1 double]   normalized spectrum of diffuse light (optical)
rad.fEsunt  = fEsunt;       % [2162x1 double]   normalized spectrum of direct light (thermal)
rad.fEskyt  = fEskyt;       % [2162x1 double]   normalized spectrum of diffuse light (thermal)

rad.Eplu_   = Eplu_;        % [61x2162 double]  upward diffuse radiation in the canopy (mW m-2 um-1) 
rad.Emin_   = Emin_;        % [61x2162 double]  downward diffuse radiation in the canopy (mW m-2 um-1) 

rad.Lo_     = Lo_;          % [2162x1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1) 
rad.Eout_   = Eout_;        % [2162x1 double]   TOC upward radiation (mW m-2 um-1) 
rad.Eouto   = Eouto;        % [1 double]        TOC spectrally integrated upward optical ratiation (W m-2) 
rad.Eoutt   = Eoutt;        % [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)

rad.Rnhs    = Rnhs;         % [1 double]        net radiation (W m-2) of shaded soil
rad.Rnus    = Rnus;         % [1 double]        net radiation (W m-2) of sunlit soil
rad.Rnhc    = Rnhc;         % [60x1 double]     net radiation (W m-2) of shaded leaves
rad.Rnuc    = Rnuc;         % [13x36x60 double] net radiation (W m-2) of sunlit leaves
rad.Pnh     = Pnhc;         % [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
rad.Pnu     = Pnuc;         % [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
rad.Pnh_Cab = Pnhc_Cab;     % [60x1 double]     net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Pnu_Cab = Pnuc_Cab;     % [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
rad.Rnh_PAR = Rnhc_PAR;     % [60x1 double]     net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
rad.Rnu_PAR = Rnuc_PAR;     % [13x36x60 double] net PAR absorbed (W m-2) of sunlit 
rad.Etoto   = Etoto;

%% APPENDIX I functions J1 and J2 (introduced for numerically stable solutions)

function J1 = calcJ1(x,m,k,LAI)
if abs(m-k)>1E-3;
    J1 = (exp(m*LAI*x)-exp(k*LAI*x))./(k-m);
else
    J1 = -.5*(exp(m*LAI*x)+exp(k*LAI*x))*LAI.*x.*(1-1/12*(k-m).^2*LAI^2.*x.^2);
end
return

function J2 = calcJ2(x,m,k,LAI)
J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))./(k+m);
return;

%% APPENDIX II function volscat

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

%% APPENDIX III function e2phot

function molphotons = e2phot(lambda,E)
%molphotons = e2phot(lambda,E) calculates the number of moles of photons
%corresponding to E Joules of energy of wavelength lambda (m)

global constants;
A = constants.A;

e           = ephoton(lambda);
photons     = E./e;
molphotons  = photons./A;
return;

function E = ephoton(lambda)
%E = phot2e(lambda) calculates the energy content (J) of 1 photon of 
%wavelength lambda (m)

global constants;
h       = constants.h;           % [J s]         Planck's constant
c       = constants.c;           % [m s-1]       speed of light
E       = h*c./lambda;           % [J]           energy of 1 photon
return;

%% APPENDIX IV function Pso

function pso    =   Psofunction(K,k,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
end