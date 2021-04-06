function [rad] = RTMz(constants,spectral,rad,soil,leafopt,canopy,gap,angles,Knu,Knh)

% function 'RTMz' calculates the small modification of TOC outgoing
% radiance due to the conversion of Violaxanthin into Zeaxanthin in leaves
%
% Author:  Christiaan van der Tol (c.vandertol@utwente.nl)
% Date:     08 Dec 2016
%           17 Mar 2020     CvdT    added cluming, mSCOPE representation
%           25 Jun 2020     PY      Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees

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
%
% Output:
%   rad         a large number of radiative fluxes: spectrally distributed
%               and integrated, and canopy radiative transfer coefficients.
%               Here, fluxes are added

%% 0.1 initialisations
wlS         = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlZ         = spectral.wlZ';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlZ); %#ok<ASGLU>
nwlZ        = length(spectral.wlZ);
nl          = canopy.nlayers;
LAI         = canopy.LAI;
iLAI        = LAI/nl;                       % LAI of a layer        [1]
litab       = canopy.litab;
lazitab     = canopy.lazitab;
lidf        = canopy.lidf;
nlazi       = length(lazitab);         % azumith angle
nlinc       = length(litab);           % inclination
nlori       = nlinc * nlazi;           % total number of leaf orientations
layers      = 1:nl;

RZ          = (leafopt.reflZ(:, iwlfi)-leafopt.refl(:, iwlfi))';
TZ          = (leafopt.tranZ(:, iwlfi)-leafopt.tran(:, iwlfi))';

Ps          = gap.Ps;
Po          = gap.Po;
Pso         = gap.Pso;
Qs          = Ps(1:end-1);

% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum
Esunf_             = rad.Esun_(iwlfi);
[Eminf_,Epluf_]    = deal(zeros(nwlZ,nl+1,2));
Eminf_(:,:,1)      = rad.Emins_(:,iwlfi)';
Eminf_(:,:,2)      = rad.Emind_(:,iwlfi)';
Epluf_(:,:,1)      = rad.Eplus_(:,iwlfi)';
Epluf_(:,:,2)      = rad.Eplud_(:,iwlfi)';

Xdd         = rad.Xdd(:,iwlfi);
rho_dd      = rad.rho_dd(:,iwlfi);
R_dd        = rad.R_dd(:,iwlfi);
tau_dd      = rad.tau_dd(:,iwlfi);
vb          = rad.vb(:,iwlfi);
vf          = rad.vf(:,iwlfi);

%% 0.2 geometric quantities

% geometric factors
deg2rad             = constants.deg2rad;
tto                 = angles.tto;
tts                 = angles.tts;
psi                 = angles.psi;
rs                  = soil.refl(iwlfi,:);           % [nwlfo]     soil reflectance
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

%% 1.0 calculation of 'flux' in observation direction
[Fmin_,Fplu_]       = deal(zeros(nl+1,nwlZ,2));
LoF_                = zeros(nwlZ,2);
laz= 1/36;
etah = Kn2Cx(Knh);

if size(Knu,2)==1
    etau = repmat(Kn2Cx(Knu),1,13,36);
else
    etau           = permute(Kn2Cx(Knu),[3 1 2]);    % make dimensions [nl,nlinc,nlazi]
    etau           = reshape(etau,nl,468);    % expand orientations in a vector >> [nl,nlori]
end
etau_lidf = bsxfun(@times,reshape(etau,nlori,nl),repmat(lidf*laz,36,1));     %[nlori,nl]
etah_lidf = bsxfun(@times,repmat(etah,1,nlori)',repmat(lidf*laz,36,1));

for k = 1:2
    [U,Y]     = deal(zeros(nl+1,nwlZ)); 
    MpluEsun  = RZ .* Esunf_*(k<2);      %
    MminEsun  = TZ .* Esunf_*(k<2);
    
    MpluEmin  = RZ  .* Eminf_(:,1:nl,k);	    % [nf,nl+1]  = (nf,ne) * (ne,nl+1)
    MpluEplu  = RZ  .* Epluf_(:,1:nl,k);
    MminEmin  = TZ  .* Eminf_(:,1:nl,k);
    MminEplu  = TZ  .* Epluf_(:,1:nl,k);
    
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
    
    %   Emitted 'flux'
    
    piLs        =   wfEs+vfEplu_u+vbEmin_u;         % sunlit for each layer
    piLd        =   vbEmin_h+vfEplu_h;              % shade leaf for each layer
    Fsmin       =   sfEs+sigfEmin_u+sigbEplu_u;     % Eq. 29a for sunlit leaf
    Fsplu       =   sbEs+sigbEmin_u+sigfEplu_u;     % Eq. 29b for sunlit leaf
    Fdmin       =   sigfEmin_h+sigbEplu_h;          % Eq. 29a for shade leaf
    Fdplu       =   sigbEmin_h+sigfEplu_h;          % Eq. 29b for shade leaf
    Femmin      =   iLAI*bsxfun(@times,Qs', Fsmin) +iLAI* bsxfun(@times,(1-Qs)',Fdmin);
    Femplu      =   iLAI*bsxfun(@times,Qs', Fsplu) +iLAI*bsxfun(@times,(1-Qs)',Fdplu);
    
    for j=nl:-1:1      % from bottom to top
        Y(j,:)  =(rho_dd(j,:).*U(j+1,:)+Femmin(:,j)')./(1-rho_dd(j,:).*R_dd(j+1,:));
        U(j,:) =tau_dd(j,:).*(R_dd(j+1,:).*Y(j,:)+U(j+1,:))+Femplu(:,j)';
    end
    for j=1:nl          % from top to bottom
        Fmin_(j+1,:,k)  = Xdd(j,:).*Fmin_(j,:,k)+Y(j,:);
        Fplu_(j,:,k)  = R_dd(j,:).*Fmin_(j,:,k)+U(j,:);
    end
    piLo1     = iLAI*piLs*Pso(1:nl);
    piLo2     = iLAI*piLd*(Po(1:nl)-Pso(1:nl));
    piLo3     = iLAI*(vb.*Fmin_(layers,:,k)  + vf.*Fplu_(layers,:,k))'*Po(1:nl);
    piLo4     = rs .* (Po(end)* Fmin_(end,:,k)');
    piLtot    = piLo1 + piLo2 + piLo3 + piLo4;
    LoF_(:,k)      = piLtot/pi;
end
Fhem_     = sum(Fplu_(1,:,:),3)';

%% output
rad.Lo_(iwlfi)          = rad.Lo_(iwlfi) + sum(LoF_,2);
rad.rso(iwlfi)          = rad.rso(iwlfi) + LoF_(:,1)./(rad.Esun_(iwlfi));
rad.rdo(iwlfi)          = rad.rdo(iwlfi) + LoF_(:,2)./(rad.Esky_(iwlfi));
rad.refl(iwlfi)         = pi*rad.Lo_(iwlfi)./(rad.Esky_(iwlfi)+rad.Esun_(iwlfi));     % [nwl] 
rad.Eout_(iwlfi)        = rad.Eout_(iwlfi) + Fhem_;

function Cx = Kn2Cx(Kn)
%Cx = 0.70*Kn;  % empirical fit by N Vilfan
Cx = 0.3187*Kn;  % empirical fit by N Vilfan (Vilfan et al, 2018, 2019)
return
