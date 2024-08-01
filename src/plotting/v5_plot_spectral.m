%% youtube video 5
% https://youtu.be/6DvOiadMA_M?si=OWPUq3KQp2nl0eKi

%% reflectance
% components
% (1) soil
% (2) leaf (per layer)
% (3) canopy

figure

wl = spectral.wlS * 1e-3; 

plot(wl, leafopt.refl(1, :), 'o', 'DisplayName', 'leafopt.refl (leaf)')
hold on
plot(wl, soil.refl, 'o', 'DisplayName', 'soil.refl (soil)')

plot(wl, rad.refl, 'x', 'DisplayName', 'rad.refl (canopy)')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')
legend

xlabel('wavelength, \mum')
title('reflectance')


%% spectral regions (see MATLAB structure "spectral" in the workspace)
%
% (1) PROSPECT 400 - 2400 nm
% (2) TIR - 2400-50000 nm (spacing 100 nm)
% (3) SIF - 640-850 nm  (also SIF excitation 400-750 nm)
% (4) xanthophyll - 500 - 600 nm

%% spectral regions

figure

w = 10; % line width

plot(spectral.wlS*1e-3, repmat(1,size(spectral.wlS)), 'Color', "#0072BD", 'LineWidth', w)
hold on
plot(spectral.wlP*1e-3, repmat(2,size(spectral.wlP)), 'Color', "#77AC30", 'LineWidth', w)

% plot(spectral.wlE*1e-3, repmat(2.5,size(spectral.wlE)), 'o', 'Color', "#77AC30", 'DisplayName', 'spectral.wlF (SIF excitation)')

plot(spectral.wlF*1e-3, repmat(3,size(spectral.wlF)), 'Color', "#EDB120", 'LineWidth', w)

plot(spectral.wlT*1e-3, repmat(4,size(spectral.wlT)), 'Color', "#D95319", 'LineWidth', w)

plot(spectral.wlZ*1e-3, repmat(5,size(spectral.wlZ)), 'Color', "#7E2F8E", 'LineWidth', w)

ylim([0, 6])
yticks([1,2,3,4,5])
yticklabels({'wlS (SCOPE)','wlP (FLUSPECT)','wlF (SIF)','wlT (TIR)', 'wlZ (PRI)'})

% legend
title('spectral regions of the SCOPE model', 'spectral.wlS')
xlabel('wavelength, \mum')
set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% radiance background
% (1) hemispherically integrated flux, E
% (2) flux in the observation direction, L, sr-1
%% L

figure

subplot(1,2,1)

wl = spectral.wlS * 1e-3; 

plot(wl, rad.Lo_, 'DisplayName', ['rad.Lo_' newline 'reflected'])
hold on

plot(wl, rad.Lot_, 'DisplayName', ['rad.Lot_' newline 'emitted in TIR wl'])

wl = spectral.wlF * 1e-3;  % um
plot(wl, rad.LoF_, 'o', 'DisplayName', ['rad.LoF_' newline 'emitted in SIF wl'])


xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')
title('Radiance (L)', 'in the observation direction')

legend('Interpreter', 'none')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% E
% figure
subplot(1,2,2)

wl = spectral.wlS * 1e-3; 

plot(wl, rad.Eout_, 'DisplayName', ['rad.Eout_' newline 'reflected in FLUSPECT wl'])
hold on

plot(wl, rad.Eoutte_, 'DisplayName', ['rad.Eoutte_' newline 'emitted in TIR wl'])

wl = spectral.wlF * 1e-3;  % um
plot(wl, rad.EoutF_, 'o', 'DisplayName', ['rad.EoutF_' newline 'emitted in SIF wl'])


xlabel('wavelength, \mum')
ylabel('E, W m-2 \mum-1')
title('Radiance (E)', 'hemispherically integrated')

legend('Interpreter', 'none')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

%% apparent reflectance

figure

subplot(2,2,1)
wl = spectral.wlS * 1e-3; 

plot(wl, rad.refl, 'DisplayName', 'rad.refl (reflectance)')
hold on
plot(wl, rad.reflapp, 'DisplayName', 'rad.reflapp (apparent reflectance)')

set(gca, 'XScale', 'log')

title('canopy reflectances')
xlabel('wavelength, \mum')
legend('Interpreter', 'none')


subplot(2,2,2)
plot(wl, rad.reflapp - rad.refl, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true reflectance', 'rad.reflapp - rad.refl', 'Interpreter','none')
xlabel('wavelength, \mum')


subplot(2,2,3)

plot(wl, rad.Lotot_, 'DisplayName', 'rad.Lotot_ (radiance)')
hold on
plot(wl, rad.Lototf_, 'DisplayName', 'rad.Lototf_ (apparent radiance)')

set(gca, 'XScale', 'log')

title('canopy radiance in observation direction')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')
legend('Interpreter', 'none')

subplot(2,2,4)

plot(wl, rad.Lototf_ - rad.Lotot_, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true radiance', 'rad.Lototf_ - rad.Lotot_', 'Interpreter','none')
xlabel('wavelength, \mum')
ylabel('L, W m-2 \mum-1 sr-1')


set(findall(gcf,'-property','FontSize'),'FontSize', 14)


for i=1:4
    subplot(2,2,i)
    xlim([640, 850]*1e-3)
end



%% xanthophyll cycle (leaf reflectance only)
% this contribution is summed with Lo_, so individual plotting is not
% easlity possible

figure

wl = spectral.wlS * 1e-3; 

plot(wl, leafopt.refl(1, :), 'DisplayName', ['leafopt.refl' newline 'leaf reflectance 100% viola'])
hold on
plot(wl, leafopt.reflZ(1, :), 'DisplayName', ['leafopt.reflZ' newline 'leaf reflectance 100% zea [light stress]'])


legend
set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')

title('Leaf reflectance', 'violaxanthin cycle')
xlabel('wavelength, \mum')

% xlim([500, 600]*1e-3)


% p = patch([wl; flip(wl)], [leafopt.refl(1, :), fliplr(leafopt.reflZ(1, :))],  'blue', 'DisplayName', 'effect');
% p.FaceColor = '#7E2F8E';
% p.EdgeColor = 'none';