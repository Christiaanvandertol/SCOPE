function [rad,profiles] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles)  
 
% function 'RTMf' calculates the spectrum of fluorescent radiance in the
% observer's direction in addition to the total TOC spectral hemispherical upward Fs flux
%
% Authors:  Wout Verhoef and Christiaan van der Tol (c.vandertol@utwente.nl)
% Date:     12 Dec 2007
% Update:   26 Aug 2008 CvdT        Small correction to matrices
%           07 Nov 2008 CvdT        Changed layout
% Update:   19 Mar 2009 CvdT        Major corrections: lines 95-96,
%                                   101-107, and 119-120.
% Update:    7 Apr 2009 WV & CvdT   Major correction: lines 89-90, azimuth
%                                   dependence was not there in previous verions (implicit assumption of
%                                   azimuth(solar-viewing) = 0). This has been corrected
% Update:   May-June 2012 WV & CvdT Add calculation of hemispherical Fs
%                                   fluxes
% Update:   Jan-Feb 2013 WV         Inputs and outputs via structures for
%                                   SCOPE Version 1.40
% Update:   Jan 2015  CvdT          Added two contributions to SIF radiance cuased by rescattering of hemispherical SIF fluxes
% Update:   Jan 2015  JAK           (from SCOPE 1.53): Improved speed by factor of 9+! (by vectorizing the summation over the 60 layers)
% Update:   Jan 2015  WV            Rearranged some arrays to smoothen the vectorizations; adjusted some internal names
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
wlF          = spectral.wlF';       % Fluorescence wavelengths
wlE          = spectral.wlE';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF); %#ok<ASGLU>
nf           = length(iwlfo);
ne           = length(iwlfi);
nl           = canopy.nlayers;
LAI          = canopy.LAI;
litab        = canopy.litab;
lazitab      = canopy.lazitab;
lidf         = canopy.lidf;
nlinc        = length(litab);
nlazi        = length(lazitab);
nlori        = nlinc * nlazi;           % total number of leaf orientations
layers       = 1:nl;

vb           = rad.vb(iwlfo);		    % added for rescattering of SIF fluxes
vf           = rad.vf(iwlfo);
 
Ps           = gap.Ps;                  % [nl+1]
Po           = gap.Po;
Pso          = gap.Pso;

Qso         = (Pso(layers) + Pso(layers+1))/2;
Qs          = (Ps(layers)  + Ps(layers+1))/2;
Qo          = (Po(layers)  + Po(layers+1))/2;
Qsho        = Qo - Qso;
 
etah         = zeros(nl,1);
etau         = zeros(nl,nlori);   % modified dimensions to facilitate vectorization
 
[Mb,Mf]      = deal(zeros(nf,ne));
LoF_         = zeros(nf,2);
Fhem_        = zeros(nf,2);
Fiprofile    = zeros(nl+1,2);
[LoF_1,LoF_2,LoF_3,LoF_4] = deal(zeros(nf,2));
 
% for speed-up the calculation only uses wavelength i and wavelength o part of the spectrum

Esunf_       = rad.Esun_(iwlfi);
Eminf_       = rad.Emin_(:,iwlfi)';          % transpose into [nf,nl+1] matrix
Epluf_       = rad.Eplu_(:,iwlfi)';
iLAI         = LAI/nl;                       % LAI of a layer        
 
%% optical quantities

rho          = leafopt.refl(iwlfo);          % [nf]     leaf/needle reflectance
tau          = leafopt.tran(iwlfo);          % [nf]     leaf/needle transmittance

if isfield(leafopt,'MbI')
    MbI          = leafopt.MbI;
    MbII         = leafopt.MbII;
    MfI          = leafopt.MfI;
    MfII         = leafopt.MfII;
else
    MbII    = leafopt.Mb;
    MfII    = leafopt.Mf;
end

rs           = soil.refl(spectral.IwlF);                  % [nf]     soil reflectance

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
 
%% calculation of fluorescence in observation direction
 
% fluorescence efficiencies from ebal, after default fqe has been applied
 
etahi = profiles.etah;
etaur = permute(profiles.etau,[3 1 2]);    % make dimensions [nl,nlinc,nlazi]
etaui = reshape(etaur,nl,nlori);           % expand orientations in a vector >> [nl,nlori]
 
% fluorescence matrices and efficiencies for PSI and PSII

[Emin_,Eplu_] = deal(zeros(nl+1,nf));
[Fmin_,Fplu_] = deal(zeros(nf,nl+1,2));
Fem = zeros(nf,2);

for PS = 2:-1:2-isfield(leafopt,'MbI')   % Do for both photosystems II and I (or alternatively for only one PS, in this case PSII)

    switch PS
        case 1, Mb = MbI;  Mf = MfI;    etah(:) = 1;          etau(:) = 1;
        case 2, Mb = MbII; Mf = MfII;   etah(:) = etahi(:);   etau(:) = etaui(:);
    end
 
    Mplu = 0.5 * (Mb+Mf);     % [nf,ne]  
    Mmin = 0.5 * (Mb-Mf);     
 
    % in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.
    % resolution assumed is 1 nm
      
    MpluEmin  = Mplu * Eminf_;	    % [nf,nl+1]  = (nf,ne) * (ne,nl+1)
    MpluEplu  = Mplu * Epluf_;	
    MminEmin  = Mmin * Eminf_;	
    MminEplu  = Mmin * Epluf_;	
    
    MpluEsun  = Mplu * Esunf_;      % integration by inproduct
    MminEsun  = Mmin * Esunf_;    
 
    xdd2    = mean(ctl2' * lidf);                           % lidf-weighted cosine squared of leaf inclination
    mn_etau = mean(reshape(etau,nl,nlinc,nlazi),3) * lidf;  % lidf-weighted mean of etau per layer [nl]
 
    % we calculate the spectrum for all individual leaves, sunlit and
    % shaded
    
    [Fmin,Fplu]   = deal(zeros(nf,nl+1));
    [G1,G2]       = deal(zeros(nl+1,1));
     
    for i = 1 : nf
        
        Qso_wfEs = Qso * reshape(absfsfo * MpluEsun(i) + fsfo * MminEsun(i),1,nlori);  % [1,nlori]
	    Qs_sfEs  = Qs *reshape(absfs * MpluEsun(i) - fsctl * MminEsun(i),1,nlori);             
	    Qs_sbEs  = Qs * reshape(absfs * MpluEsun(i) + fsctl * MminEsun(i),1,nlori);             

       
        
        Mplu_i(layers) = MpluEmin(i,layers) + MpluEplu(i,layers+1);  % [1,nl]
        Mmin_i(layers) = MminEmin(i,layers) - MminEplu(i,layers+1);

	    sigfEmini_sigbEplui = Mplu_i' - xdd2 * Mmin_i';			     % [nl]
	    sigbEmini_sigfEplui = Mplu_i' + xdd2 * Mmin_i';
        
        Qso_Mplu  = Qso .* Mplu_i';                          % [nl]
        Qso_Mmin  = Qso .* Mmin_i';
        Qsho_Mplu = Qsho .* Mplu_i';
        Qsho_Mmin = Qsho .* Mmin_i';
        
        Qso_vEd  = Qso_Mplu * reshape(absfo,1,nlori) + Qso_Mmin * reshape(foctl,1,nlori);
        Qsh_vEd  = Qsho_Mplu * reshape(absfo,1,nlori) + Qsho_Mmin * reshape(foctl,1,nlori);
               
        % Directly observed fluorescence contributions from sunlit and
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
        
        
        % Transformed SIF fluxes calculated numerically
             
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
        
        Fplu_(i,:,PS) = (rinf*G1+G2)/(1-rinf^2);
        Fmin_(i,:,PS) = (rinf*G2+G1)/(1-rinf^2);
        
        Fhem_(i,PS) = Fplu_(i,1,PS);

	    % The following contributions are coming from:
        
        %   3) Rescattered SIF at observed leaves
 	    %   4) SIF flux reflected by observed soil
             
        piLo3       = iLAI * ((vb(i)*Fmin_(i,layers,PS) + vf(i)*Fplu_(i,layers+1,PS)) * Qo); 
        piLo4       = rs(i) * Fmin_(i,nl+1,PS) * Po(nl+1);
        
        piLtoti     = piLo1 + piLo2 + piLo3 + piLo4;
        LoF_(i,PS)  = piLtoti/pi;  
        
        LoF_1(i,PS) = piLo1/pi;
        LoF_2(i,PS) = piLo2/pi;
        LoF_3(i,PS) = piLo3/pi;
        LoF_4(i,PS) = piLo4/pi;
    end
   
    Fem(:,PS) = iLAI*sum(Fplu+Fmin,2);
    
    for ilayer = 1:nl+1
        Fiprofile(ilayer,PS) = 0.001 * Sint(Fplu(:,ilayer),spectral.wlF);
    end
end
rad.Fem_  = Fem(:,1) + Fem(:,2);
rad.Fhem_  = Fem(:,1) + Fem(:,2);
rad.LoF_  = LoF_(:,1)  + LoF_(:,2);
if isfield(leafopt,'MbI')
    rad.LoF1_ = LoF_(:,1);
    rad.LoF2_ = LoF_(:,2);
end
rad.Fhem_ = Fhem_(:,1) + Fhem_(:,2);
rad.Fmin  = sum(Fmin_,3);
rad.Fplu  = sum(Fplu_,3);
rad.LoF_sunlit      = LoF_1;
rad.LoF_shaded      = LoF_2;
rad.LoF_scattered   = LoF_3;
rad.LoF_soil        = LoF_4;

profiles.fluorescence   = Fiprofile(:,1) + Fiprofile(:,2);
 
rad.Eoutf = 0.001 * Sint(sum(Fhem_,2),spectral.wlF); 
rad.Eminf_ = Emin_;
rad.Epluf_ = Eplu_;
