%% youtube video 8
% https://youtu.be/f543CkCSIUg?si=etRXbAaf7dUtxtX1

%% mSCOPE compute different reflectance for leaf layers

figure

wl = spectral.wlS * 1e-3; 

for i=1:canopy.nlayers
    plot(wl, leafopt.refl(i, :), 'DisplayName', 'leaf reflectance')
    hold on
end


set(gca, 'XScale', 'log')


xlabel('wavelength, \mum')
ylabel('reflectance, -')
title('leaf reflectance', 'per layer')

set(findall(gcf,'-property','FontSize'),'FontSize', 14)

%%

nl = canopy.nlayers;
Ps = gap.Ps;

y = canopy.xl(2:end); % vertical depth in canopy, 0 - top, -1 - bottom 

%% radiation

figure

subplot(2, 2, 1)

plot(Ps, canopy.xl, 'DisplayName', 'gap.Ps (sunlit)')
hold on
plot(1-Ps, canopy.xl, 'DisplayName', '1 - gap.Ps (shaded)')

legend
ylabel('depth in canopy')
xlabel('probability, -')
title('probability gap fraction', 'per leaf layer')



for p = 2:4
    
    subplot(2, 2, p)

    % aPAR in umol photons m-2 s-1
    name = 'aPAR';
    full_name = {"absorbed PAR", "rad.Pnu, rad.Pnh", "per leaf layer"};
    units = '\mumol photons m-2 s-1';
    sunlit = rad.Pnu;
    shaded = rad.Pnh;
    % fapar canopy.Rntot_PAR / rad.EPAR


    if p == 3
        name = 'Ag';
        full_name =  {"gross photosynthesis", "bcu.Ag, bch.Ag", "per leaf layer"};
        units = '\mumol CO_2 m-2 s-1';
        sunlit = bcu.Ag;
        shaded = bch.Ag;
    elseif p == 4
        % eta
        name = 'eta';
        full_name = {"eta = Fs / Fo", "bcu.eta, bch.eta", "per leaf layer"};
        units = '-';
        sunlit = bcu.eta;
        shaded = bch.eta;
    end
    
    if options.lite == 0
        sunlit = meanleaf(canopy, sunlit, 'angles');
    end

    plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
    hold on
    plot(shaded, y, 'o-', 'DisplayName', 'shaded')
    
    sunlit_w  = Ps(1:nl) .* sunlit;
    shaded_w = (1-Ps(1:nl)) .* shaded;
    

    % plot(sunlit_w, y, '--', 'DisplayName', 'sunlit1d weighted', 'Color', "#0072BD")
    % hold on
    % plot(shaded_w, y, '--', 'DisplayName', 'shaded weighted', 'Color', 	"#D95319")
    
    plot(sunlit_w + shaded_w, y, 'x-', 'DisplayName', 'sum proflie')
    
    legend
    
    ylabel('depth in canopy')
    xlabel(sprintf('%s, %s', name, units))
    title(full_name)
end

set(findall(gcf,'-property','FontSize'),'FontSize', 14)