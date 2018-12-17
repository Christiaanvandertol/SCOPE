%% plot vertical profiles

% Note layer 1 is top, layer 60 is bottom: to check plot Rn, it is higher on top
% plot(profiles.Rn1d, 1:60, 'o-')

% plot(profiles.Pn1d, 1:60, 'o-')  % strange shape
plot(flip(profiles.Tc1d), 1:60, 'o-', 'MarkerFaceColor', 'r')
hold on
plot(flip(profiles.Tcu1d), 1:60, 'o-', 'MarkerFaceColor', 'y')
plot(flip(profiles.Tch), 1:60, 'o-', 'MarkerFaceColor', 'b')

lgd = legend('canopy', 'sunlit leaves', 'shaded leaves', 'Location', 'eastoutside');
title(lgd,'component')
xlabel('Temperature ^{o}C')
% xlim([640, 850])
% ylim([1, Inf])
ylabel('flipped layer number')
title({'Vertical profile of canopy temperature'})

% annotation('textbox', [0, 0.5, 0.5, 0.5], 'string', 'My Text')

set(findall(gcf,'-property','FontSize'),'FontSize',18)
