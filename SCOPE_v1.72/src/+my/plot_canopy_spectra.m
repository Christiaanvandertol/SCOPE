%% canopy spectra

% optical

o_wl = [400:2400] / 1000;
o_i = 1:2001;
t_wl = [2500 : 100 : 15000, 16000 : 1000 : 50000] / 1000;
t_i = 2002:2162;

subplot(2,2,1)
plot(o_wl, rad.Eout_(o_i), 'r', 'LineWidth', 1.5)
xlabel('wavelength, um')
ylabel('W m-2 um-1')
xlim([0.4, 2.4])
ylim([0, 300])
title('Canopy hemispherical radiance (optical)')

subplot(2,2,2)
plot(o_wl, rad.Lo_(o_i), 'r', 'LineWidth', 1.5)
xlabel('wavelength, um')
ylabel('W m-2 um-1 sr-1')
xlim([0.4, 2.4])
ylim([0, 300])
title('Canopy directional radiance (optical)')

% thermal

subplot(2,2,3)
plot(t_wl, rad.Eoutte_(t_i), 'r', 'LineWidth', 1.5)
xlabel('wavelength, um')
ylabel('W m-2 um-1')
xlim([2.5, 50])
ylim([0, 30])
title('Canopy hemispherical radiance (thermal)')

subplot(2,2,4)
plot(t_wl, rad.Lot_(t_i), 'r', 'LineWidth', 1.5)
xlabel('wavelength, um')
ylabel('W m-2 um-1 sr-1')
xlim([2.5, 50])
ylim([0, 30])
title('Canopy directional radiance (thermal)')


set(findall(gcf,'-property','FontSize'),'FontSize',18)