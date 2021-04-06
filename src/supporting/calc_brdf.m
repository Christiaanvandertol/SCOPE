function directional = calc_brdf(constants,options,directional,spectral,angles,atmo,soil,leafopt,canopy,meteo,thermal,bcu,bch)

% simulates observations from a large number of viewing angles
% modified: 30 April 2020, CvdT, removed repeated angle combinations.

%% input
tts             = angles.tts;
psi_hot         = [0    ; 0     ;0     ;0     ;0    ;2  ;358];  % [noa_o]           angles for hotspot oversampling
tto_hot         = [tts  ; tts+02;tts+04;tts-02;tts-4;tts;tts];  % [noa_o]           angles for hotspot oversampling

psi_plane       = [000*ones(6,1);180*ones(6,1);090*ones(6,1);270*ones(6,1)];%       angles for plane oversampling
tto_plane       = [10:10:60     , 10:10:60    , 10:10:60    , 10:10:60]';   %       angles for plane oversampling

psi             = [directional.psi; psi_hot; psi_plane];
tto             = [directional.tto; tto_hot; tto_plane];

[~,u]           = unique([psi tto],'rows');
directional.psi = psi(u);
directional.tto = tto(u);
na              = length(u);

%% allocate memory
directional.brdf_       = zeros(length(spectral.wlS),na);      % [nwlS, no of angles]  
directional.Eoutte      = zeros(1,na);                         % [1, no of angles]  
directional.BrightnessT = zeros(1,na);                         % [1, no of angles] 
directional.LoF_        = zeros(length(spectral.wlF),na);      % [nwlF, no of angles] 
directional.Lot_        = zeros(length(spectral.wlT),na);      % [nwlF, no of angles] 

%% other preparations
directional_angles = angles;

%% loop over the angles
for j=1:na
    
    %optical BRDF
    directional_angles.tto = directional.tto(j);
    directional_angles.psi = directional.psi(j);
    [directional_rad,directional_gap] = RTMo(spectral,atmo,soil,leafopt,canopy,directional_angles,constants,meteo,options);
    directional.refl_(:,j)  = directional_rad.refl;         % [nwl]            reflectance (spectral) (nm-1)
    directional.rso_(:,j)  = directional_rad.rso;           % [nwl]            BRDF (spectral) (nm-1)
    
    % thermal directional brightness temperatures (Planck)
    
    if options.calc_planck
        directional_rad = RTMt_planck(spectral,directional_rad,...
            soil,leafopt,canopy,directional_gap,...
            thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh);
        directional.Lot_(:,j)      = directional_rad.Lot_(spectral.IwlT) + directional_rad.Lo_(spectral.IwlT);        % [nwlt]           emitted plus reflected diffuse     radiance at top
    end
    
    if options.calc_fluor      
        directional_rad = RTMf(constants,spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,bcu.eta,bch.eta);    
        directional.LoF_(:,j)              = directional_rad.LoF_;
    end % {if calcfluor}
    if options.calc_xanthophyllabs
        directional_rad = RTMz(constants,spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,bcu.Kn,bch.Kn);
        directional.Lo_ = directional_rad.Lo_;
    end
                
end % {for angles}
