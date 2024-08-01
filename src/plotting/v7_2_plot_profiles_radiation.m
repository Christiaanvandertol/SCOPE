%% youtube video 7
% https://youtu.be/lhmP78f5tYg?si=O2i3m-WipOSlY8xP

%% plot profiles
% this script works with the structures available in the workspace after
% the SCOPE run

%% radiation
% we just sample top (1) and bottom (end) 
% Emin_ - downward diffuse
% Eplu_ - upwar diffuse

figure

wl = spectral.wlS * 1e-3; % um

subplot(2, 2, 1)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling (E-) radiation', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 2)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) radiation', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')



wl = spectral.wlT * 1e-3;

subplot(2, 2, 3)

per_wl = rad.Emint_;

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling thermal (E-) radiation', 'rad.Emint_ (from RTMt_sb or RTMt_planck)', 'emitted by canopy'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 4)

per_wl = rad.Eplut_;

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) thermal radiation', 'rad.Eplut_ (from RTMt_sb or RTMt_planck)', 'emitted by canopy'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%% per layer

figure

wl = spectral.wlS * 1e-3; % um

for p=1:2
    per_wl = rad.Emin_;
    per_wl(:, spectral.IwlT) = per_wl(:, spectral.IwlT) + rad.Emint_;
    what = {'downwelling (E-) radiation', 'rad.Emin_ + rad.Emint_'};
    
    if p == 2
        per_wl = rad.Eplu_;
        per_wl(:, spectral.IwlT) = per_wl(:, spectral.IwlT) + rad.Eplut_;
        what = {'upwelling (E+) radiation', 'rad.Eplu_ + rad.Eplut_'};
    end
    

    subplot(1, 2, p)
    plot(wl, per_wl(1, :), 'DisplayName', 'top')
    hold on
    for i=1:(nl - 1)
        if rem(i, 5) == 0  % each 5th layer
            plot(wl, per_wl(i, :), 'DisplayName', sprintf('layer %d', i))
        end
    
    end
    plot(wl, per_wl(end, :), 'DisplayName', 'bottom')
    
    legend
    set(gca, 'XScale', 'log')
    title(what, 'Interpreter', 'none')
    
    ylabel('W m-2 \mum-1')
    xlabel('wavelength, \mum')
end


set(findall(gcf,'-property','FontSize'),'FontSize', 14)


%% components

figure

wl = spectral.wlS * 1e-3; % um

subplot(2, 2, 1)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'downwelling (E-) radiation', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')


subplot(2, 2, 2)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

plot(wl, top, 'DisplayName', 'top')
hold on
plot(wl, bottom, 'DisplayName', 'bottom')

legend

ylabel('W m-2 \mum-1')
xlabel('wavelength, \mum')
title({'upwelling (E+) radiation', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')

set(gca, 'XScale', 'log')

wl = spectral.wlS * 1e-3; % um


subplot(2, 2, 3)

per_wl = rad.Emin_;
% per_wl = rad.Emins_;  % diffuse that used to be direct
% per_wl = rad.Emind_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);



% plot(wl, top)
hold on
plot(wl, rad.Esun_, 'Color', "#EDB120")
plot(wl, rad.Esky_, 'Color', "#0072BD")

stacked = [rad.Emins_(end,:); rad.Emind_(end,:)]; 
a = area(wl, stacked');
a(1).FaceColor = "#EDB120";
a(2).FaceColor = "#0072BD";

set(gca, 'XScale', 'log')

legend({'rad.Esun_ (direct solar light)', 'rad.Esky_ (diffuse solar light)',...
    ['rad.Eplus_' newline 'diffuse that used to be direct'], ...
    ['rad.Emind_' newline 'diffuse that came from the Sun']}, 'Interpreter', 'none')

title({'components of downwelling (E-) radiation [bottom]', 'rad.Emin_ (from RTMo)', 'diffuse solar and atmosphere'}, ...
    'Interpreter', 'none')
ylabel('W m-2 \mum-1')

subplot(2, 2, 4)

per_wl = rad.Eplu_;
% per_wl = rad.Eplus_;  % diffuse that used to be direct
% per_wl = rad.Eplud_;  % diffuse that came from the Sun

top = per_wl(1, :);
bottom = per_wl(end, :);

stacked = [rad.Eplus_(1,:); rad.Eplud_(1,:)]; 

% plot(wl, top)
% hold on
a = area(wl, stacked');
a(1).FaceColor = "#EDB120";
a(2).FaceColor = "#0072BD";

set(gca, 'XScale', 'log')

legend({['rad.Eplus_' newline 'diffuse that used to be direct'], ...
    ['rad.Eplud_' newline 'diffuse that came from the Sun']}, 'Interpreter', 'none')

title({'components of upwelling (E+) radiation [top]', 'rad.Eplu_ (from RTMo)', 'reflected solar and atmosphere'}, ...
    'Interpreter', 'none')
ylabel('W m-2 \mum-1')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)