%% youtube video 7
% https://youtu.be/lhmP78f5tYg?si=O2i3m-WipOSlY8xP

%% plot profiles
% this script works with the structures available in the workspace after
% the SCOPE run

%% contribution of sunlit and shaded leaves depends on their fraction
% computed as the weighted [by the gap fraction gap.Ps] sum L386 RTMo 
% profiles.Pn1d  = ((1-Ps(1:nl)).*Pnhc  + Ps(1:nl).*(Pnu1d));  

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
        % net radiation
        name = 'Rn';
        full_name = {"shortwave net radiation", "rad.Rnuc, rad.Rnhc" "per leaf layer"};
        units = 'W m-2';
        sunlit = rad.Rnuc;
        shaded = rad.Rnhc;
    elseif p == 4
        % thermal net radiation
        name = 'Rnt';
        full_name = {"thermal net radiation", "rad.Rnuct, rad.Rnhct", "per leaf layer"};
        units = 'W m-2';
        sunlit = rad.Rnuct;
        shaded = rad.Rnhct;
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

%% biochemical output

figure


for p = 1:5
    
    subplot(2, 3, p)

    % Ag (gross leaf photosynthesis) in umol CO2 m-2 s-1
    name = 'Ag';
    full_name = {"gross photosynthesis", "bcu.Ag, bch.Ag", "per leaf layer"};
    units = '\mumol CO_2 m-2 s-1';
    sunlit = bcu.Ag;
    shaded = bch.Ag;

    if p == 2
        name = 'Vcmax';
        full_name = {"maximum carboxylation capactiy", "bcu.Vcmax, bch.Vcmax", "per leaf layer"};
        units = '\mumol CO_2 m-2 s-1';
        sunlit = bcu.Vcmax;
        shaded = bch.Vcmax;
    elseif p == 3
        % eta
        name = 'eta';
        full_name = {"eta = Fs / Fo", "bcu.eta, bch.eta", "per leaf layer"};
        units = '-';
        sunlit = bcu.eta;
        shaded = bch.eta;

    elseif p == 4
        % RH
        name = 'RH_{leaf}';
        full_name = {"relative humidity (RH=0.64)", "bcu.RH, bch.RH", "per leaf layer"};
        units = '-';
        sunlit = bcu.RH;
        shaded = bch.RH;
    elseif p == 5
        % gsw
        name = 'gsw';
        full_name = {"stomatal conductance", "bcu.gs, bch.gs", "per leaf layer"};
        units = 'mol H_2O m-2 s-1';
        sunlit = bcu.gs;
        shaded = bch.gs;
    
    % elseif p == 6
    %     % Ja
    %     name = 'Ja';
    %     full_name = {"ETR actual", "per leaf layer"};
    %     units = '\mumol photons m-2 s-1';
    %     sunlit = bcu.Ja;
    %     shaded = bch.Ja;
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
    
    % legend
    
    ylabel('depth in canopy')
    xlabel(sprintf('%s, %s', name, units))
    title(full_name)
end


%% thermal

subplot(2,3,6)

name = 'T_{leaf}';
full_name = sprintf("temperatures (Ta = %.0f^oC)", meteo.Ta);
units = 'C';
sunlit = thermal.Tcu;
shaded = thermal.Tch;

if options.lite == 0
    sunlit = meanleaf(canopy, sunlit, 'angles');
end


plot(sunlit, y, 'o-', 'DisplayName', 'sunlit1d')
hold on
plot(shaded, y, 'o-', 'DisplayName', 'shaded')


sunlit_w  = Ps(1:nl) .* (sunlit .^ 4);
shaded_w = (1-Ps(1:nl)) .* (shaded .^ 4);

plot((sunlit_w + shaded_w).^ (1/4), y, 'x-', 'DisplayName', 'weighted proflie')

legend

ylabel('depth in canopy')
xlabel(sprintf('%s, %s', name, units))
title({full_name, 'thermal.Tcu, thermal.Tch', 'per leaf layer'})


set(findall(gcf,'-property','FontSize'),'FontSize', 14)






