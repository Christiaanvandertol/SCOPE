function [rad,profiles] = RTMz(spectral,rad,soil,leafopt,canopy,gap,angles,profiles)

% function 'RTMz' calculates the small modification of TOC outgoing
% radiance due to the conversion of Violaxanthin into Zeaxanthin in leaves
%
% Author:  Christiaan van der Tol (c.vandertol@utwente.nl)
% Date:     08 Dec 2016     
%
% The inputs and outputs are structures. These structures are further
% specified in a readme file.
%
% Input:
%   spectral    information about wavelengths and resolutions
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%   soil        soil properties
%   leafopt     leaf optical properties
%   canopy      canopy properties (such as LAI and height)
%   gap         probabilities of direct light penetration and viewing
%   angles      viewing and observation angles
%   profiles    vertical profiles of fluxes
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluorescence fluxes are added
%% 0.0 globals
global constants

%% initialisations
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlZ          = spectral.wlZ';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlZ); %#ok<ASGLU>
nf           = length(iwlfi);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nlinc        = length(litab);
nlazi        = length(lazitab);
nlori        = nlinc * nlazi;           % total number of leaf orientations
layers       = 1:nl;

RZ          = leafopt.reflZ(iwlfi)-leafopt.refl(iwlfi);
TZ          = leafopt.tranZ(iwlfi)-leafopt.tran(iwlfi);

vb           = rad.vb(iwlfi);		    % added for rescattering of SIF fluxes
vf           = rad.vf(iwlfi);

Ps           = gap.Ps;                  % [nl+1]
Po           = gap.Po;
Pso          = gap.Pso;

Qso         = (Pso(layers) + Pso(layers+1))/2;
Qs          = (Ps(layers)  + Ps(layers+1))/2;
Qo          = (Po(layers)  + Po(layers+1))/2;
Qsho        = Qo - Qso;

etah         = zeros(nl,1);
etau         = zeros(nl,nlori);   % modified dimensions to facilitate vectorization

LoZ_         = zeros(nf,1);
Zhem_        = zeros(nf,1);

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum

Esunf_       = rad.Esun_(iwlfi);
Eminf_       = rad.Emin_(:,iwlfi)';          % transpose into [nf,nl+1] matrix
Epluf_       = rad.Eplu_(:,iwlfi)';
iLAI         = LAI/nl;                       % LAI of a layer

%% optical quantities

rho          = leafopt.refl(iwlfi);          % [nf]     leaf/needle reflectance
tau          = leafopt.tran(iwlfi);          % [nf]     leaf/needle transmittance
rs           = soil.refl(iwlfi);                  % [nf]     soil reflectance

% geometric factors

deg2rad      = constants.deg2rad;
tto          = angles.tto;
tts          = angles.tts;
psi          = angles.psi;

cos_tto      = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto      = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts      = cos(tts*deg2rad);             % cos solar angle
sin_tts      = sin(tts*deg2rad);             % sin solar angle

cos_ttli     = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli     = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils    = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo    = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% geometric factors for all leaf angle/azumith classes
cds          = cos_ttli*cos_tts*ones(1,nlazi) + sin_ttli*sin_tts*cos_phils;  % [nlinc,nlazi]
cdo          = cos_ttli*cos_tto*ones(1,nlazi) + sin_ttli*sin_tto*cos_philo;  % [nlinc,nlazi]
fs           = cds/cos_tts;                                                  % [nlinc,nlazi]
absfs        = abs(fs);                                                      % [nlinc,nlazi]
fo           = cdo/cos_tto;                                                  % [nlinc,nlazi]
absfo        = abs(fo);                                                      % [nlinc,nlazi]
fsfo         = fs.*fo;                                                       % [nlinc,nlazi]
absfsfo      = abs(fsfo);                                                    % [nlinc,nlazi]
foctl        = fo.*(cos_ttli*ones(1,nlazi));                                 % [nlinc,nlazi]
fsctl        = fs.*(cos_ttli*ones(1,nlazi));                                 % [nlinc,nlazi]
ctl2         = cos_ttli.^2*ones(1,nlazi);                                    % [nlinc,nlazi]

%% calculation of le in observation direction

% Cx as a function of Kn

etahi = Kn2Cx(profiles.Knh);
etaur = permute(Kn2Cx(profiles.Knu),[3 1 2]);    % make dimensions [nl,nlinc,nlazi]
etaui = reshape(etaur,nl,nlori);           % expand orientations in a vector >> [nl,nlori]

[Fmin_,Fplu_] = deal(zeros(nf,nl+1));
etah(:) = etahi(:);   etau(:) = etaui(:);

MpluEmin  = (RZ*ones(1,nl+1)) .* Eminf_;	    % [nf,nl+1]  = (nf,ne) * (ne,nl+1)
MpluEplu  = (RZ*ones(1,nl+1))  .* Epluf_;
MminEmin  = (TZ*ones(1,nl+1))  .* Eminf_;
MminEplu  = (TZ*ones(1,nl+1))  .* Epluf_;

MpluEsun  = RZ .* Esunf_;      %
MminEsun  = TZ .* Esunf_;

xdd2    = mean(ctl2' * lidf);                           % lidf-weighted cosine squared of leaf inclination
mn_etau = mean(reshape(etau,nl,nlinc,nlazi),3) * lidf;  % lidf-weighted mean of etau per layer [nl]

% we calculate the spectrum for all individual leaves, sunlit and
% shaded

[Fmin,Fplu]   = deal(zeros(nf,nl+1));
[G1,G2]       = deal(zeros(nl+1,1));
[Mplu_i,Mmin_i] = deal(zeros(nl,1));

for i = 1 : nf
    
    Qso_wfEs = Qso * reshape(absfsfo * MpluEsun(i) + fsfo * MminEsun(i),1,nlori);  % [1,nlori]
    Qs_sfEs  = Qs *reshape(absfs * MpluEsun(i) - fsctl * MminEsun(i),1,nlori);
    Qs_sbEs  = Qs * reshape(absfs * MpluEsun(i) + fsctl * MminEsun(i),1,nlori);
    
    Mplu_i(layers) = MpluEmin(i,layers) + MpluEplu(i,layers+1);  % [1,nl]
    Mmin_i(layers) = MminEmin(i,layers) - MminEplu(i,layers+1);
    
    sigfEmini_sigbEplui = Mplu_i - xdd2 * Mmin_i;			     % [nl]
    sigbEmini_sigfEplui = Mplu_i + xdd2 * Mmin_i;
    
    Qso_Mplu  = Qso .* Mplu_i;                          % [nl]
    Qso_Mmin  = Qso .* Mmin_i;
    Qsho_Mplu = Qsho .* Mplu_i;
    Qsho_Mmin = Qsho .* Mmin_i;
    
    Qso_vEd  = Qso_Mplu * reshape(absfo,1,nlori) + Qso_Mmin * reshape(foctl,1,nlori);
    Qsh_vEd  = Qsho_Mplu * reshape(absfo,1,nlori) + Qsho_Mmin * reshape(foctl,1,nlori);
    
    % Directly observed radiation contributions from sunlit and
    % shaded leaves
    
    piLs = mean(reshape(etau .*(Qso_wfEs + Qso_vEd),nl,nlinc,nlazi),3) * lidf;
    piLd = etah .* (mean(reshape(Qsh_vEd,nl,nlinc,nlazi),3) * lidf);
    
    piLo1 = iLAI * sum(piLs) ;
    piLo2 = iLAI * sum(piLd);
    
    Qs_Fsmin = mean(reshape(etau .* Qs_sfEs,nl,nlinc,nlazi),3) * lidf ...
        + Qs .* mn_etau .* sigfEmini_sigbEplui;
    Qs_Fsplu = mean(reshape(etau .* Qs_sbEs,nl,nlinc,nlazi),3) * lidf ...
        + Qs .* mn_etau .* sigbEmini_sigfEplui;
    
    Qd_Fdmin = (1-Qs) .* etah .* sigfEmini_sigbEplui;
    Qd_Fdplu = (1-Qs) .* etah .* sigbEmini_sigfEplui;
    
    Fmin(i,layers+1) = Qs_Fsmin + Qd_Fdmin;
    Fplu(i,layers)   = Qs_Fsplu + Qd_Fdplu;
    
    t2   = xdd2 * (rho(i)-tau(i))/2;
    att  = 1-(rho(i)+tau(i))/2+t2;
    sig  = (rho(i)+tau(i))/2+t2;
    m    = sqrt(att^2-sig^2);
    rinf = (att - m)/sig;
    fac  = 1 - m * iLAI;
    facs = (rs(i)-rinf)/(1-rs(i)*rinf);
    
    % Transformed radiance calculated numerically
    
    G1(1) = 2; Gnew = 0;    % (to ensure we will enter the loop the first time)
    
    dF1 = (Fmin(i,layers+1) + rinf * Fplu(i,layers)) * iLAI;    % Thanks to JAK
    dF2 = (rinf * Fmin(i,layers+1) + Fplu(i,layers)) * iLAI;    % These are the source functions
    
    while abs(Gnew-G1(1)) > 1e-3
        G1(1) = Gnew;
        for j=2:nl+1
            G1(j) = fac * G1(j-1) + dF1(j-1);
        end
        G2(nl+1) = G1(nl+1) * facs;
        for j=nl:-1:1
            G2(j) = fac * G2(j+1) + dF2(j);
        end
        Gnew = -rinf * G2(1);
    end
    
    % Inverse transformation to get back the hemispherical fluxes
    
    Fplu_(i,:) = (rinf*G1+G2)/(1-rinf^2);
    Fmin_(i,:) = (rinf*G2+G1)/(1-rinf^2);
    
    Zhem_(i) = Fplu_(i,1);
    
    % The following contributions are coming from:
    
    %   3) Rescattered radiation of  observed leaves
    %   4) radiation reflected by observed soil
    
    piLo3       = iLAI * ((vb(i)*Fmin_(i,layers) + vf(i)*Fplu_(i,layers+1)) * Qo);
    piLo4       = rs(i) * Fmin_(i,nl+1) * Po(nl+1);
    
    piLtoti     = piLo1 + piLo2 + piLo3 + piLo4;
    LoZ_(i)  = piLtoti/pi;
    
end

rad.Lo_(iwlfi)          = rad.Lo_(iwlfi)+LoZ_;
rad.Eout_(iwlfi)        = rad.Eout_(iwlfi) + Zhem_;

function Cx = Kn2Cx(Kn)
Cx = 0.70*Kn;  % empirical fit by N Vilfan
return