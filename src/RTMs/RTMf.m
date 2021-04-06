function rad = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,etau,etah)

% function 'RTMf' calculates the spectrum of fluorescent radiance in the
% observer's direction and also the TOC spectral hemispherical upward Fs flux
%
% Authors:  Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl.nl)
% Date:     12 Dec 2007
% Update:   26 Aug 2008 CvdT        small correction to matrices
%           07 Nov 2008 CvdT        changed layout
% Update:   19 Mar 2009 CvdT        major corrections: lines 95-96,
%                                   101-107, and 119-120.
% Update:    7 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
%                                   dependence was not there in previous verions (implicit assumption of
%                                   azimuth(solar-viewing) = 0). This has been corrected
% Update:   May-June 2012 WV & CvdT Add calculation of hemispherical Fs
%                                   fluxes
% Update:   Jan-Feb 2013 WV         Inputs and outputs via structures for
%                                   SCOPE Version 1.40
% Update:   Aug-Oct 2016 PY         Re-write the calculation of emitted SIF
%                                   of each layer. It doesnt use loop at
%                                   all. with the function bsxfun, the
%                                   calculation is much faster
% Update:   Oct 2017-Feb 2018 PY    Re-write the RTM of fluorescence
% Update:   Jan 2020 CvdT           Modified to include 'lite' option,
%                                   mSCOPE representation
% Update:   25 Jun 2020 PY          Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees

% Table of contents of the function:
%   0       preparations
%       0.0     globals
%       0.1     initialisations
%       0.2     geometric quantities
%       0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
%   1.0     calculation of fluorescence flux in observation direction
%
% Usage: [rad] = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,etau,etah)
%
% The input and output are structures. These structures are further
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
%   etau,etah   relative fluorescence emission efficiency for sunlit and
%               shaded leaves, respectively
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluorescence fluxes are added


%% 0.1 initialisations
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlF          = (640:4:850)';%spectral.wlF';       % Fluorescence wavelengths
wlE          =   (400:5:750)'; %spectral.wlE';    % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF); %#ok<ASGLU>
nf           = length(iwlfo);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nlazi        = length(lazitab);         % azumith angle
nlinc        = length(litab);           % inclination
nlori        = nlinc * nlazi;           % total number of leaf orientations
layers       = 1:nl;

Ps           = gap.Ps;
Po           = gap.Po;
Pso          = gap.Pso;

Qs          = Ps(1:end-1);

[MpluEmin   ,...
    MpluEplu   , ...
    MminEmin   , ...
    MminEplu   , ...
    MpluEsun   , ...
    MminEsun]       = deal(zeros(nf,nl));

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_             = rad.Esun_(iwlfi);
Eminf_             = rad.Emin_(:,iwlfi)';          % transpose into [nwlfo,nl] matrix
Epluf_             = rad.Eplu_(:,iwlfi)';
iLAI               = LAI/nl;                       % LAI of a layer        [1]

Xdd         = rad.Xdd(:,iwlfo);
rho_dd      = rad.rho_dd(:,iwlfo);
R_dd        = rad.R_dd(:,iwlfo);
tau_dd      = rad.tau_dd(:,iwlfo);
vb          = rad.vb(:,iwlfo);
vf          = rad.vf(:,iwlfo);

%% 0.2 geometric quantities
Mb                = leafopt.Mb;
Mf                = leafopt.Mf;

% geometric factors
deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;
rs                  = soil.refl(iwlfo,:);           % [nwlfo]     soil reflectance
cos_tto             = cos(tto*deg2rad);             % cos observation zenith angle
sin_tto             = sin(tto*deg2rad);             % sin observation zenith angle

cos_tts             = cos(tts*deg2rad);             % cos solar angle
sin_tts             = sin(tts*deg2rad);             % sin solar angle

cos_ttli            = cos(litab*deg2rad);           % cos leaf inclinaation angles
sin_ttli            = sin(litab*deg2rad);           % sin leaf inclinaation angles
cos_phils           = cos(lazitab*deg2rad);         % cos leaf azimuth angles rel. to sun azi
cos_philo           = cos((lazitab-psi)*deg2rad);   % cos leaf azimuth angles rel. to viewing azi

%% 0.3 geometric factors for all leaf angle/azumith classes
cds                 = cos_ttli*cos_tts*ones(1,36) + sin_ttli*sin_tts*cos_phils;  % [nli,nlazi]
cdo                 = cos_ttli*cos_tto*ones(1,36) + sin_ttli*sin_tto*cos_philo;  % [nli,nlazi]
fs                  = cds/cos_tts;                                               % [nli,nlazi]
absfs               = abs(fs);                                                   % [nli,nlazi]
fo                  = cdo/cos_tto;                                               % [nli,nlazi]
absfo               = abs(fo);                                                   % [nli,nlazi]
fsfo                = fs.*fo;                                                    % [nli,nlazi]
absfsfo             = abs(fsfo);                                                 % [nli,nlazi]
foctl               = fo.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
fsctl               = fs.*(cos_ttli*ones(1,36));                                 % [nli,nlazi]
ctl2                = cos_ttli.^2*ones(1,36);                                    % [nli,nlazi]

% reshape all the variables
absfs               = reshape(absfs,nlori,1);                                    % [nlori,1]
absfo               = reshape(absfo,nlori,1);                                    % [nlori,1]
fsfo                = reshape(fsfo,nlori,1);                                     % [nlori,1]
absfsfo             = reshape(absfsfo,nlori,1);                                  % [nlori,1]
foctl               = reshape(foctl,nlori,1);                                    % [nlori,1]
fsctl               = reshape(fsctl,nlori,1);                                    % [nlori,1]
ctl2                = reshape(ctl2,nlori,1);                                     % [nlori,1]

%% 1.0 calculation of fluorescence flux in observation direction

% fluorescence matrices and efficiencies
[U,Fmin_,Fplu_] =deal(zeros(nl+1,size(leafopt.Mb,1)));

Mplu = 0.5*(Mb+Mf);    % [nwlfo,nwlfi]
Mmin = 0.5*(Mb-Mf);    % [nwlfo,nwlfi]

% in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.
for j = 1:nl  
    ep = constants.A*ephoton(wlF*1E-9,constants);
    MpluEmin(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Eminf_(:,j),constants));          % [nwlfo,nl+1]
    MpluEplu(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Epluf_(:,j),constants));          % [nwlfo,nl+1]
    MminEmin(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Eminf_(:,j),constants));          % [nwlfo,nl+1]
    MminEplu(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Epluf_(:,j),constants));          % [nwlfo,nl+1]
    MpluEsun(:,j)   = ep.*(Mplu(:,:,j)* e2phot(wlE*1E-9,Esunf_,constants));               % integration by inproduct
    MminEsun(:,j)   = ep.*(Mmin(:,:,j)* e2phot(wlE*1E-9,Esunf_,constants));               % integration by inproduct
end

laz= 1/36;

if size(etau,2)<2
    etau = repmat(etau,1,13,36);
    etau = permute(etau,[2 3 1]);
end

etau_lidf = bsxfun(@times,reshape(etau,nlori,nl),repmat(lidf*laz,36,1));     %[nlori,nl]
etah_lidf = bsxfun(@times,repmat(etah,1,nlori)',repmat(lidf*laz,36,1));

wfEs      =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfsfo)),MpluEsun) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsfo)),MminEsun);
sfEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
sbEs     =  bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfs)),MpluEsun) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,fsctl)),MminEsun);
vfEplu_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEplu);
vfEplu_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEplu);
vbEmin_h  = bsxfun(@times,sum(bsxfun(@times,etah_lidf,absfo)),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,foctl)),MminEmin);
vbEmin_u  = bsxfun(@times,sum(bsxfun(@times,etau_lidf,absfo)),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,foctl)),MminEmin);
sigfEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
sigfEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
sigbEmin_h  = bsxfun(@times,sum(etah_lidf),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEmin);
sigbEmin_u  = bsxfun(@times,sum(etau_lidf),MpluEmin) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEmin);
sigfEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
sigfEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) -...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);
sigbEplu_h  = bsxfun(@times,sum(etah_lidf),MpluEplu) +...
    bsxfun(@times,sum(bsxfun(@times,etah_lidf,ctl2)),MminEplu);
sigbEplu_u  = bsxfun(@times,sum(etau_lidf),MpluEplu) +...
    bsxfun(@times,sum(bsxfun(@times,etau_lidf,ctl2)),MminEplu);

%   Emitted fluorescence
piLs        =   wfEs+vfEplu_u+vbEmin_u;         % sunlit for each layer
piLd        =   vbEmin_h+vfEplu_h;              % shade leaf for each layer
Fsmin       =   sfEs+sigfEmin_u+sigbEplu_u;     % Eq. 29a for sunlit leaf
Fsplu       =   sbEs+sigbEmin_u+sigfEplu_u;     % Eq. 29b for sunlit leaf
Fdmin       =   sigfEmin_h+sigbEplu_h;          % Eq. 29a for shade leaf
Fdplu       =   sigbEmin_h+sigfEplu_h;          % Eq. 29b for shade leaf
Femmin      =   iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin);
Femplu      =   iLAI*bsxfun(@times,Qs', Fsplu) +iLAI* bsxfun(@times,(1-Qs)',Fdplu);

for j=nl:-1:1      % from bottom to top
    Y(j,:)  =(rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
    U(j,:) =tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';
end

for j=1:nl          % from top to bottom
    Fmin_(j+1,:)  = Xdd(j,:).*Fmin_(j,:)+Y(j,:);
    Fplu_(j,:)  = R_dd(j,:).*Fmin_(j,:)+U(j,:);
end
piLo1     = iLAI*piLs*Pso(1:nl);
piLo2     = iLAI*piLd*(Po(1:nl)-Pso(1:nl));
piLo3     = iLAI*(vb.*Fmin_(layers,:)  + vf.*Fplu_(layers,:))'*Po(1:nl);
piLo4     = rs .* Fmin_(nl+1,:)' * Po(nl+1);

piLtot      = piLo1 + piLo2 + piLo3 + piLo4;
LoF_        = piLtot/pi;
Fhem_       = Fplu_(1,:)';

method = 'spline';  % M2020a name
if verLessThan('matlab', '9.8')
    method = 'splines';
end

rad.LoF_    = interp1(wlF,LoF_,spectral.wlF',method);
rad.EoutF_   = interp1(wlF,Fhem_,spectral.wlF',method);

rad.LoF_sunlit      = interp1(wlF,piLo1/pi,spectral.wlF',method);
rad.LoF_shaded      = interp1(wlF,piLo2/pi,spectral.wlF',method);
rad.LoF_scattered   = interp1(wlF,piLo3/pi,spectral.wlF',method);
rad.LoF_soil        = interp1(wlF,piLo4/pi,spectral.wlF',method);

rad.EoutF   = 0.001 * Sint(Fhem_,wlF);
rad.LoutF   = 0.001 * Sint(LoF_,wlF);

[rad.F685,iwl685]  = max(rad.LoF_(1:55));
rad.wl685 = spectral.wlF(iwl685);
if iwl685 == 55; [rad.F685,rad.wl685] = deal(NaN); end
[rad.F740,iwl740]  = max(rad.LoF_(70:end));
rad.wl740 = spectral.wlF(iwl740+69);
rad.F684  = rad.LoF_(685-spectral.wlF(1));
rad.F761  = rad.LoF_(762-spectral.wlF(1));
return

