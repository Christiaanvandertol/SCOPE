path_output =  fullfile('..', '..', Output_dir);

spectra = dlmread(fullfile(path_output, 'apparent_reflectance.csv'), ',', 2, 0);
pars = readtable(fullfile(path_output, 'pars_and_input_short.csv'), "VariableNamesLine", 1, 'CommentStyle', '#');

wl = spectral.wlS * 1e-3; 

figure

plot(wl, spectra)

leg = legend(num2str(pars{:, 2}));
title(leg, pars.Properties.VariableNames{2})

xlabel('wavelength, \mum')
set(findall(gcf,'-property','FontSize'),'FontSize', 14)
set(gca, 'XScale', 'log')


%% Cs or Cbrown?


% prospect_pro = dataSpec_PRO();   %% in input\fluspect_parameters

figure

plot(optipar.wl, optipar.Ks, 'o-')
hold on
plot(prospect_pro(:, 1), prospect_pro(:, 6), 'x-')

legend({'SAC Ks', 'SAC Kbrown'})

xlabel('wavelength, nm')
set(findall(gcf,'-property','FontSize'),'FontSize', 14)


