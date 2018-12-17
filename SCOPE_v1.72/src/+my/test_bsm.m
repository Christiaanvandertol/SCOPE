load bsm_input.mat

wl = 400:2400;


%% from BSM

smc = [5, 10, 20, 30, 55];
smc = [0.5, 0.10, 0.20, 0.30, 0.55];
for el=smc
    soilpar.SMC = el;
    res = BSM(soilpar,spec,emp);
    subplot(1,2,1)
    plot(wl, res, 'LineWidth', 1.5)
    hold on
end
lgd=legend(arrayfun(@string, smc));
title(lgd,'Soil moisture')
xlabel('wavelength, nm')
xlim([400, 2400])
ylim([0, 1])
ylabel('soil reflectance')
title({'Soil reflectance simulated with BSM model', ...
    'B=0.5, lat=25, lon=45, SMC=25, film=0.015'})


%% from file 

subplot(1,2,2)
plot(wl, rsfile(:,2), 'LineWidth', 1.5)
hold on
plot(wl, rsfile(:,3), 'LineWidth', 1.5)
plot(wl, rsfile(:,4), 'LineWidth', 1.5)
lgd=legend('1', '2', '3');
title(lgd,'Column number')
xlabel('wavelength, nm')
ylim([0, 1])
xlim([400, 2400])
ylabel('soil reflectance')
title('Soil reflectance from \it../data/input/soil\_spectrum/soilnew.txt')
set(findall(gcf,'-property','FontSize'),'FontSize',18)

