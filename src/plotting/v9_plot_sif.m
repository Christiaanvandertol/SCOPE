%% youtube video 9
% https://youtu.be/CpE07VEjs4Q?si=r46ltIpa94ljWff_

%% fluorescence background
% emitted 
% (1) at photosystem level (inside chloroplasts)
% (2) in all directions (hemispherical)
% reabsorbed 
% (1) within the leaf it was emitted
% (2) within the canopy (by other leaves)
% (3) by soil


%% hemispherical (theoretical + real)

figure
wl = spectral.wlF;  % nm

% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')
area(wl, rad.EoutFrc_, 'DisplayName', ['rad.EoutFrc_', newline ...
    'photosystem level (without within-leaf reabsorption)' newline ...
    'fluorescence_ReabsCorr.csv'])
hold on

area(wl, rad.Femleaves_, 'DisplayName', ['rad.Femleaves_', newline ...
    'all leaves forward+backward (without canopy/soil reabsorption)' newline ...
    'fluorescence_AllLeaves.csv'])

area(wl, rad.EoutF_, 'DisplayName', ['rad.EoutF_', newline ...
    'SIF hemispherical (at the top of canopy)' newline...
    'fluorescence_hemis.csv'])


legend('Interpreter','none')


xlabel('wavelength, nm')
ylabel('E, mW m-2 nm-1')  % == W m-2 um-1

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

title({'SIF', 'hemispherically integrated'})

%% stacked: in observation direction

figure
wl = spectral.wlF;  % nm

plot(wl, rad.LoF_, 'o', 'DisplayName', 'in observation direction')
hold on

stacked = [rad.LoF_sunlit, rad.LoF_shaded, rad.LoF_scattered, rad.LoF_soil];
area(wl, stacked)

legend({['rad.LoF_' newline 'fluorescence.csv'],...
    'rad.LoF_sunlit', 'rad.LoF_shaded', 'rad.LoF_scattered', 'rad.LoF_soil'})

legend('Interpreter','none')

xlabel('wavelength, nm')
ylabel('L, mW m-2 nm-1 sr-1')  % == W m-2 um-1

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

title({'SIF', 'in the observation direction'})


%% escape probability sigmaF

figure

subplot(2, 2, 1)

wl = spectral.wlF;  % nm

plot(wl, rad.LoF_ * pi, 'o')


xlabel('wavelength, nm')
ylabel('L * pi, mW m-2 nm-1')

title({'SIF in the observation direction (measured)', 'rad.LoF_ multiplied by pi'}, 'Interpreter', 'none')

%
subplot(2, 2, 3)

plot(wl, rad.sigmaF)

xlabel('wavelength, nm')
ylabel('sigmaF, -')
title({'SIF escape probability (rad.sigmaF)', 'escape from photosystem to sensor'}, 'Interpreter', 'none')

% plot(wl, rad.EoutFrc_ .* rad.sigmaF' / pi)
% hold on
% plot(wl,  rad.LoF_)


%%%%% reabsorption corrected

subplot(2, 2, [2, 4])


% from ETR and fqe of PSII
plot(wl, rad.EoutFrc_, 'x')


xlabel('wavelength, nm')
ylabel('E, mW m-2 nm-1')
title({'SIF hemispherically integrated (theoretical)', ...
    'rad.EoutFrc_ = rad.LoF_ * pi / rad.sigmaF'}, 'Interpreter', 'none')


set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%% apparent reflectance and radiance

figure

subplot(2,2,1)
wl = spectral.wlS * 1e-3; 

plot(wl, rad.refl, 'DisplayName', 'rad.refl (reflectance)')
hold on
plot(wl, rad.reflapp, 'DisplayName', 'rad.reflapp (apparent reflectance)')

set(gca, 'XScale', 'log')

title('canopy reflectances')
xlabel('wavelength, \mum')
ylabel('reflectance, -')
legend('Interpreter', 'none')


subplot(2,2,2)
plot(wl, rad.reflapp - rad.refl, 'Color', "#EDB120")

set(gca, 'XScale', 'log')

title('apparent minus true reflectance', 'rad.reflapp - rad.refl', 'Interpreter','none')
xlabel('wavelength, \mum')
ylabel('reflectance, -')


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