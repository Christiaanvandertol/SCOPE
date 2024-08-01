%% youtube video 6
% https://youtu.be/bIuc4kKLAnk?si=ZSJ_37z82vpG4NlG

%% input (to be adjusted by the user)

path_output =  fullfile('..', '..', Output_dir, 'directional');
% path_output =  "../../output/direction_2024-06-10-1150/directional";

lst = dir(path_output);
fnames = {lst.name};

Angles = load(fullfile(path_output, fnames{contains(fnames, 'Angles')}));

% optical
wl_target = 555; 
Values = load(fullfile(path_output, fnames{contains(fnames, 'refl')}));
% 
% % SIF
% wl_target = 760;  
% Values = load(fullfile(path_output, fnames{contains(fnames, 'Fluorescence')}));
% 
% % TIR
% wl_target = 4000;  
% Values = load(fullfile(path_output, fnames{contains(fnames, 'Thermal radiances')}));


%% grid creation
obs_zenith                  = Angles(1,:);
obs_azimuth               	= Angles(2,:);

obs_azimuth180              = obs_azimuth-360*(obs_azimuth>180);

obs_zenith_i                = 0:90;
obs_azimuth_i               = -180:1:180;

[Obs_Zenith_i,Obs_Azimuth_i]= meshgrid(obs_zenith_i,obs_azimuth_i);

% conversion to polar coordinates
% https://nl.mathworks.com/matlabcentral/answers/716343-creating-polar-mesh-in-matlab#answer_597478
% r = obs_zenith * pi/180
% phi = obs_azimuth180  * pi/180 + pi/2

x                           = obs_zenith  *pi/180 .* cos(obs_azimuth180  *pi/180+pi/2);
y                           = obs_zenith  *pi/180 .* sin(obs_azimuth180  *pi/180+pi/2);

X_i                         = Obs_Zenith_i*pi/180 .* cos(Obs_Azimuth_i*pi/180+pi/2);
Y_i                         = Obs_Zenith_i*pi/180 .* sin(Obs_Azimuth_i*pi/180+pi/2);


%% wavelength selection

wl = Values(:,1    );  % nm
i_wl = find(wl == wl_target);

assert(~isempty(i_wl), 'wavelength %d nm was not found. Adjust wl_target?', wl_target)

BRDF                        = Values(:,2:end);

% regridding 

BRDF_i                      = griddata(x,y,BRDF(i_wl,:),X_i,Y_i,'v4');


%% 2d plot

figure

xli = .5*pi*[0 -1.15  -.1 1];
yli = .5*pi*[1  0 -1.05 0];


subplot(2, 2, 1)
z = pcolor(X_i,Y_i,BRDF_i); 
hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color', 'black')
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',12.727272727272727,'Color', 'black');
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end

title(sprintf('BRDF @ %.2d nm', wl_target),'FontSize',14)
colorbar
colormap jet

axis off

%% 3d plot


% figure
subplot(2, 2, 2)
surf(X_i,Y_i,BRDF_i); 
shading interp
colormap jet
colorbar

title(sprintf('BRDF @ %.2d nm', wl_target),'FontSize',14)


%% principal plane: psi 0 -> 180


% figure
subplot(2, 2, 3)


% forward scattering
i = obs_azimuth == 180;
x_f = obs_zenith(i);
y_f = BRDF(i_wl,i);

% back scattering
i = obs_azimuth == 0;
x_b = obs_zenith(i) * -1;
y_b = BRDF(i_wl,i);

x = [x_f, x_b];
[~, i2] = sort(x);
y = [y_f, y_b];

plot(x(i2), y(i2))
xlabel('observation zenith angle (tts)')

title(sprintf('BRDF @ %.2d nm\nprincipal plane (phi 0 -> 180)\n"vertical" where the Sun is', wl_target),'FontSize',14)


%% cross-principal plane: psi 90 -> 270


% figure
subplot(2, 2, 4)

% forward scattering
i = obs_azimuth == 270;
x_f = obs_zenith(i);
y_f = BRDF(i_wl,i);

% back scattering
i = obs_azimuth == 90;
x_b = obs_zenith(i) * -1;
y_b = BRDF(i_wl,i);

x = [x_f, x_b];
[~, i2] = sort(x);
y = [y_f, y_b];

plot(x(i2), y(i2))
xlabel('observation zenith angle (tts)')

title(sprintf('BRDF @ %.2d nm\ncross-principal plane (phi 90 -> 270)\n perpendicular to principal', wl_target),'FontSize',14)

