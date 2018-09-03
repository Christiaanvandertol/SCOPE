function directional = calc_brdf(options,directional,spectral,angles,rad,atmo,soil,leafopt,canopy,meteo,profiles,thermal)

global constants

%% input
tts                 = angles.tts;
noa                 = directional.noa;

psi_hoversampling   = [0    ; 0     ;0     ;0     ;0    ;2  ;358];  % [noa_o]           angles for hotspot oversampling
tto_hoversampling   = [tts  ; tts+02;tts+04;tts-02;tts-4;tts;tts];  % [noa_o]           angles for hotspot oversampling

noah_o              = size(tto_hoversampling,1);                    % [1]               number of oversampling angles

psi_poversampling   = [000*ones(6,1);180*ones(6,1);090*ones(6,1);270*ones(6,1)];%       angles for plane oversampling
tto_poversampling   = [10:10:60     , 10:10:60    , 10:10:60    , 10:10:60]';   %       angles for plane oversampling

noap_o              = size(tto_poversampling,1);                    % [1]               number of oversampling angles

directional.psi     = [directional.psi;psi_hoversampling;psi_poversampling];   % [..]              observer azimuth angle
directional.tto     = [directional.tto;tto_hoversampling;tto_poversampling];   % [..]              observer zenith  angle

%% allocate memory
directional.brdf_       = zeros(length(spectral.wlS),noa + noah_o+noap_o);      % [nwlS, noa+noa_o+noap_o]  
directional.Eoutte      = zeros(1,noa + noah_o+noap_o);                         % [1, noa+noa_o+noap_o]  
directional.BrightnessT = zeros(1,noa + noah_o+noap_o);                         % [1, noa+noa_o+noap_o] 
directional.LoF_        = zeros(length(spectral.wlF),noa + noah_o+noap_o);      % [nwlF, noa+noa_o+noap_o] 
directional.Lot_        = zeros(length(spectral.wlT),noa + noah_o+noap_o);      % [nwlF, noa+noa_o+noap_o] 

%% other preparations
directional_angles = angles;

%% loop over the angles
for j=1:(noa+noah_o+noap_o)
    
    %optical BRDF
    directional_angles.tto = directional.tto(j);
    directional_angles.psi = directional.psi(j);   
    [directional_rad,directional_gap] = RTMo(spectral,atmo,soil,leafopt,canopy,directional_angles,meteo,rad,options);
    directional.brdf_(:,j)  = directional_rad.rso;%Lo_./E_tot;          % [nwl]            BRDF (spectral) (nm-1)
    
    % thermal directional brightness temperatures (Planck)
    if options.calc_planck   
        directional_rad = RTMt_planck(spectral,directional_rad,...
                                      soil,leafopt,canopy,directional_gap,directional_angles,...
                                      thermal.Tcu,thermal.Tch,thermal.Ts(1),thermal.Ts(1),1);
        directional.Lot_(:,j)      = directional_rad.Lot_(spectral.IwlT);        % [nwlt]           emitted diffuse     radiance at top      
        
    else            %thermal directional brightness temperatures (Stefan-Boltzmann)
        directional_rad                 = RTMt_sb(spectral,directional_rad,...
                                            soil,leafopt,canopy,directional_gap,directional_angles,thermal.Tcu,thermal.Tch,thermal.Ts(1),thermal.Ts(1),1);
        directional.Lot(j)              = directional_rad.Eoutte;
        directional.BrightnessT(j)      = (pi*rad.Lot/constants.sigmaSB)^0.25;
    end
    
    if options.calc_fluor
        directional_rad = RTMf(spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,profiles);    
        directional.LoF_(:,j)              = directional_rad.LoF_;
    end % {if calcfluor}
    if options.calc_xanthophyllabs
        [directional_rad] = RTMz(spectral,directional_rad,soil,leafopt,canopy,directional_gap,directional_angles,profiles);
    end
                
end % {for wavelengths}
