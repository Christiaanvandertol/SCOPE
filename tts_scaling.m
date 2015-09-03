function [P,ftts,s3] = tts_scaling(target,wl,ts,s3)

%% this is the directory where the LUT is located
Dir = '../output/sens_tts_2015-05-11-2303/';

%% the files of the lookup table are loaded
Q = dlmread([Dir 'fluxes.dat'],'',2,0);                 % fluxes
p = dlmread([Dir 'pars_and_input_short.dat'],'',1,0);   % parameters
f = dlmread([Dir 'fluorescence.dat'],'',2,0);   % parameters
fhem = dlmread([Dir 'fluorescence_hemis.dat'],'',2,0);

%%
spi = 2*target+1 +(wl==685);
switch target
    case 0
        if wl>0
            Avec        = f(:,wl-639);              % fluorescence
        else
            Avec = Q(:,21);
%             Avec = zeros(size(fhem,1),1);
%             for k = 1:size(fhem,1)
%                 Avec(k)        = 0.001 * Sint(fhem(k,:),(640:850));
%             end
        end
    otherwise
        Avec        = Q(:,11);              % photosynthesis
end

tts         = p(:,1);               % solar zenith angle

Anorm       = Avec./Avec(tts == 30);
Anorm       = Anorm(tts<75);
tts         = tts(tts<75);
tts         = tts/180*pi;

P = polyfit(sin(tts),Anorm,4);

ttsi = sin((5:1:80)'/180*pi);
Am = P(1)*ttsi.^4 + P(2)*ttsi.^3 + P(3)*ttsi.^2 + P(4)*ttsi + P(5);
F3 = figure(3); hold on
%xlabel('sin(\theta_s')

ts = sin(ts/180*pi);
ftts = P(1)*ts.^4 + P(2)*ts.^3 + P(3)*ts.^2 + P(4)*ts + P(5);

if nargin>3
    
    set(F3,'Position',[360 216 265 706])
    s3(spi) = subplot(3,1,spi);
    plot(sin(tts),Anorm,'kx'); hold on
    plot(ttsi,Am,'k')
    if spi == 3
        resizefigure(s3,1,3,.22,.08,0,.05)
    end
else
    s3 = 0;
end

switch target
    case 0,
        switch wl
            case 760, ylabel('f_\theta_{s,SIF760}'),
            case 685, ylabel('f_\theta_{s,SIF685}')
        end
    case 1, ylabel('f_\theta_{s,A}'), xlabel('sin(\theta_s)')
end

