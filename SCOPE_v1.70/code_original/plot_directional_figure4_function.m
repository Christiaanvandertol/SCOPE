function plot_directional_figure4_function(directory)
% Use: plot_directional_figure4(directory) makes BRDF, BFDF and
% bidirectional temperature polar plots from a SCOPE output directory
% (string 'directory') of directional data.

files = dir(directory);

spfig3 = zeros(4,1);
h = zeros(3,1);
for m = 1:3
    subplot(3,1,m)
    for k = 1:4
        plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',m~=6*[1 1 1])
        text(20*k/180*pi-.2,0,num2str(20*k),'FontSize',14,'Color',m~=6*[1 1 1]);
    end
    if m<4

    end
    axis off
end


Anglesfile          = files(3).name;
ValuesTfile         = files(6).name;
ValuesBRDFfile      = files(4).name;
ValuesFluorfile     = files(5).name;

Angles = load([directory Anglesfile]);
ValuesT = load([directory ValuesTfile]);
ValuesBRDF = load([directory ValuesBRDFfile]);
ValuesFluor = load([directory ValuesFluorfile]);

obs_zenith                  = Angles(1,:);
obs_azimuth               	= Angles(2,:);
sol_zenith                  = Angles(3,1);

wl                          = ValuesBRDF(:,1    )*1E-3;
wlF                         = ValuesFluor(:,1    )*1E-3;
%Tb                          = ValuesT(:,1:end);
BRDF                        = ValuesBRDF(:,2:end);
Fluor                       = ValuesFluor(:,2:end);

obs_azimuth                 = obs_azimuth-360*(obs_azimuth>180);
obs_zenith_i                = 0:90;
obs_azimuth_i               = -180:1:180;
[Obs_Zenith_i,Obs_Azimuth_i]= meshgrid(obs_zenith_i,obs_azimuth_i);

wl_i                        = 0.4;
[v,i_wl]                    = min(abs(wl-.8));


[v,j_wl]                    = min(abs(wlF-.685));
[v,j_wl2]                    = min(abs(wlF-.740));
[v,j_wl3]                    = min(abs(wlF-.755));


%obs_elevation               = 90-obs_zenith;
%obs_r                       = ones(size(obs_elevation));
%Obs_Elevation_i             = 90-Obs_Zenith_i;
%Obs_R_i                     = ones(size(Obs_Zenith_i));
% [x  ,y ,z ]                 = sph2cart(obs_azimuth*pi/180,obs_elevation*pi/180,obs_r);
% [X_i,Y_i]                   = sph2cart(Obs_Azimuth_i*pi/180,Obs_Elevation_i*pi/180,Obs_R_i);
x                           = obs_zenith  *pi/180 .* cos(obs_azimuth  *pi/180+pi/2);
y                           = obs_zenith  *pi/180 .* sin(obs_azimuth  *pi/180+pi/2);

X_i                         = Obs_Zenith_i*pi/180 .* cos(Obs_Azimuth_i*pi/180+pi/2);
Y_i                         = Obs_Zenith_i*pi/180 .* sin(Obs_Azimuth_i*pi/180+pi/2);

BRDF_i                      = griddata(x,y,BRDF(i_wl,:),X_i,Y_i,'v4');
%Tb_i                        = griddata(x,y,Tb,X_i,Y_i,'v4');
Fluor_i                     = griddata(x,y,Fluor(j_wl,:),X_i,Y_i,'v4');
Fluor_i2                     = griddata(x,y,Fluor(j_wl2,:),X_i,Y_i,'v4');
Fluor_i3                    = griddata(x,y,Fluor(j_wl3,:),X_i,Y_i,'v4');

%%
F3 = figure(3);



xli = .5*pi*[0 -1.15  -.1 1];
yli = .5*pi*[1  0 -1.05 0];

spfig3(1) = subplot(1,3,1);
z = pcolor(X_i,Y_i,BRDF_i); hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',j~=1*[1 1 1])
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',j~=1*[1 1 1]);
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end
%if j == 1
    text(-1.7,1.8,'BRDF','FontSize',14)
    h(3) = colorbar;
%end
axis off

spfig3(2) = subplot(1,3,2);
z = pcolor(X_i,Y_i,Fluor_i); hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',[1 1 1])
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',[1 1 1]);
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end
%if j == 1
    text(-1.7,1.8,'Fluor @ 685 nm (W m^{-2}\mum^{-1}sr^{-1})','FontSize',14)
    h(3) = colorbar;
%end
axis off

spfig3(3) = subplot(1,3,3);
z = pcolor(X_i,Y_i,Fluor_i2); hold on
set(z,'LineStyle','none')
%   set(gca,'clim',[830 900])
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',[1 1 1])
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',[1 1 1]);
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end
%if j == 1
    text(-1.7,1.8,'Fluor @ 740 nm (W m^{-2}\mum^{-1}sr^{-1})','FontSize',14)
    h(3) = colorbar;
%end
axis off
% 
% spfig3(4) = subplot(1,4,4);
% z = pcolor(X_i,Y_i,Fluor_i3); hold on
% set(z,'LineStyle','none')
% %   set(gca,'clim',[520 555])
% for k = 1:4
%     plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',[1 1 1])
%     text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',[1 1 1]);
%     text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
% end
% %if j == 1
%     text(-1.7,1.8,'Fluor @ 755 nm (W m^{-2}\mum^{-1}sr^{-1})','FontSize',14)
%     h(3) = colorbar;
% %end
axis off
%%
%set(h(:,:),'location','southoutside','FontSize',14)
resizefigure(spfig3,3,1,.07,.1,.1,.12, .9, .88)