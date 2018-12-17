%% plot leaf quantiles

% wl = 400:2400;
% 
% % subplot(1,2,1)
% plot(wl, leafopt.reflZ(1:2001), 'r', 'LineWidth', 1.5)
% hold on
% plot(wl, leafopt.refl(1:2001), 'b', 'LineWidth', 1.5)
% 
% plot(wl, 1 - leafopt.tranZ(1:2001), 'r', 'LineWidth', 1.5)
% plot(wl, 1 - leafopt.tran(1:2001), 'b', 'LineWidth', 1.5)
% 
% lgd = legend('zeaxanthin', 'violaxanthin');
% title(lgd,'xanthopyll state')
% xlabel('wavelength, nm')
% ylabel('R, 1-T')
% xlim([400, 2400])
% title('Leaf optical properties simulated with Fluspect')
% 
% % subplot(1,2,2)
% % text(0.5, 0.5, {'Here can be your model', 'of leaf thermal emittance'} , 'HorizontalAlignment', 'center', 'FontSize',14)
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',18)