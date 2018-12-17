%% plot fluorescence
% wl = spectral.wlF;
% subplot(1,2,1)
% p = plot(wl, rad.Fhem_, 'o-', 'MarkerFaceColor', 'r');
% hold on
% plot(wl, rad.Femtot, 'o-', 'MarkerFaceColor', 'b')
% plot(wl, rad.Fem_, 'o-', 'MarkerFaceColor', 'g')
% lgd = legend('canopy', 'photosystems', 'leaves');
% title(lgd,'Level')
% xlabel('wavelength, nm')
% xlim([640, 850])
% ylim([0, 90])
% ylabel('hemispherical fluorescence, W m-2 um-1')
% title('Simulated hemispherical fluorescence quantiles')
% 
% % c = get(p,'Color')  % to get corresponding line color for canopy
% 
% subplot(1,2,2)
% plot(wl, rad.LoF_, 'o-', 'MarkerFaceColor', 'r') % , 'Color', [0.9290    0.6940    0.1250])
% xlabel('wavelength, nm')
% xlim([640, 850])
% ylim([0, 90])
% ylabel('directional fluorescence, W m-2 um-1 sr-1')
% title({'Top of canopy (TOC) directional fluorescence'})
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',18)

%% fluorescence and its contributors (with PSI)
% wl = spectral.wlF;
% 
% subplot(1,2,1)
% plot(wl, rad.LoF_, 'o-', 'MarkerFaceColor', 'r')
% hold on
% plot(wl, sum(rad.LoF_sunlit,2), 'o-', 'MarkerFaceColor', 'y')
% plot(wl, sum(rad.LoF_shaded,2), 'o-', 'MarkerFaceColor', 'b')
% LoF_scat = sum(rad.LoF_scattered,2)+sum(rad.LoF_soil,2);
% plot(wl, LoF_scat, 'o-', 'MarkerFaceColor', 'g')
% lgd = legend('measured TOC', 'sunlit leaves', 'shaded leaves', 'scattered');
% title(lgd,'contributors')
% xlabel('wavelength, nm')
% xlim([640, 850])
% ylim([0, Inf])
% ylabel('directional fluorescence, W m-2 um-1 sr-1')
% title({'TOC directional fluorescence', 'and its contributing parts'})
% 
% subplot(1,2,2)
% plot(wl, rad.LoF_, 'o-', 'MarkerFaceColor', 'r')
% hold on
% plot(wl, rad.LoF1_, 'o-', 'MarkerFaceColor', 'y')
% plot(wl, rad.LoF2_, 'o-', 'MarkerFaceColor', 'b')
% lgd = legend('measured TOC', 'PSI', 'PSII');
% title(lgd,'contributors')
% xlabel('wavelength, nm')
% xlim([640, 850])
% ylim([0, Inf])
% ylabel('directional fluorescence, W m-2 um-1 sr-1')
% title({'TOC directional fluorescence', 'and its contributing parts'})
% 
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% 

