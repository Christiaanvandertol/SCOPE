function PA = tts_scaling2

%% this is the directory where the LUT is located
Dir = '../output/sens_tts_2015-05-11-2303/';


%% the files of the lookup table are loaded
Q = dlmread([Dir 'fluxes.dat'],'',2,0);                 % fluxes
p = dlmread([Dir 'pars_and_input_short.dat'],'',1,0);   % parameters
f = dlmread([Dir 'fluorescence.dat'],'',2,0);   % parameters

%%
PA = zeros(5,3);
s = zeros(3,1);
for target = 1:3
    
    switch target
        case 1
            Avec        = f(:,760-649);              % fluorescence
        case 2
            Avec        = f(:,685-649);              % fluorescence
        case 3
            Avec        = Q(:,11);              % photosynthesis
    end
    
    tts         = p(:,1);               % solar zenith angle
    
    Anorm       = Avec./Avec(tts == 30);
    Anorm       = Anorm(tts<75);
    tts         = tts(tts<75);
    tts         = tts/180*pi;
    
    F4 = figure(4); hold on
    set(F4,'Position',[360 216 265 706])
    s(target) = subplot(3,1,target);
    plot(sin(tts),Anorm,'kx');
    
    P = polyfit(sin(tts),Anorm,4);
    
    ttsi = sin((5:1:80)'/180*pi);
    Am = P(1)*ttsi.^4 + P(2)*ttsi.^3 + P(3)*ttsi.^2 + P(4)*ttsi + P(5);
    hold on
    plot(ttsi,Am,'k')
    
    switch target
        case 1, ylabel('f_\theta_{s,SIF760}')
        case 2, ylabel('f_\theta_{s,SIF685}')
        case 3, ylabel('f_\theta_{s,A}'), xlabel('sin(\theta_s)')
    end
    
    %ts = sin(ts/180*pi);
    %ftts = P(1)*ts.^4 + P(2)*ts.^3 + P(3)*ts.^2 + P(4)*ts + P(5);
    
    PA(:,target) = P;
    
end

resizefigure(s,1,3,.22,.07,.02,.03);
