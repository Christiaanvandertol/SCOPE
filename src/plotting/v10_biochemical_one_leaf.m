%% youtube video 10
% https://youtu.be/COM89KOGKMo?si=1xHLTl1sPibPoyCP

%% canopy
fV = 1;  % the exponent of the vertical decline of Vcmax in the canopy

%% options_adj

options_adj = struct();

options_adj.apply_T_corr = 0;

%% constant

constants_adj = struct();

constants_adj.R    = 8.314;      % [J mol-1K-1]  Molar gas constant
constants_adj.rhoa = 1.2047;     % [kg m-3]      Specific mass of air
constants_adj.Mair = 28.96;      % [g mol-1]     Molecular mass of dry air


%% meteo

meteo_leaf = struct();

% 1 W m-2 ~ 4.6 umol photons m-2 s-1
meteo_leaf.Q = 500;  % umol photons m-2 s-1, aPAR by chlorophyll
meteo_leaf.T = 20;     % deg C,     leaf temperature
meteo_leaf.Cs = 410;   % ppm,       leaf boundary layer CO2 concentration
meteo_leaf.eb = 15;    % hPa,       actual atmospheric vapour pressure  meteo.ea
meteo_leaf.Oa = 209;  % per mille, atmospheric O2 concentraion  %meteo.Oa;
meteo_leaf.p = 970;    % hPa,       atmospheric pressure, meteo.p;

%% leaf biochemical parameters

leafbio_adj = struct();

% in Magniani C4 CO2 response of ps, Fs, eta is not present
leafbio_adj.Type = 'C3'; 

leafbio_adj.Vcmax25 = 60;
leafbio_adj.RdPerVcmax25 = 0.015;

leafbio_adj.BallBerrySlope = 8;
leafbio_adj.BallBerry0 = 0.01;


% CvdT model
% Ф(Kp) + Ф(Kf) + Ф(Kd) + Ф(Kn) = 1
% Kn is modelled based on the degree of light saturation: 'x'
leafbio_adj.Kn0 = 2.48;
leafbio_adj.Knalpha = 2.83;
leafbio_adj.Knbeta = 0.1140;

leafbio_adj.stressfactor = 1;  % for Vcmax

leafbio_adj.TDP = define_temp_response_biochem; % temperature response C3 and C4 according to CLM4 model


% Magniani model
% Ф(Kp) + Ф(Kf) + Ф(Kd) + Ф(Kn) = 1
% Kp and Kd are adjusted witht the respective stress factors
leafbio_adj.qLs = 1; % fraction of functional reaction centres (RC)
leafbio_adj.kNPQs = 0; % rate constant of sustained thermal dissipation, normalized to (Kf+Kd)

leafbio_adj.beta = 0.51; % fraction of photons partitioned to PSII

leafbio_adj.Tyear = 15; % mean annual temperature



%% light response curves

meteo_leaf.Q = linspace(1, 2000, 100);
meteo_leaf.Cs = 410;
meteo_leaf.T = 20; 
x_values = meteo_leaf.Q;
x_label = 'aPAR\_Cab, \mumol photons m-2 s-1';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% CO2 response curves

meteo_leaf.Q = 500;
meteo_leaf.Cs = linspace(0, 1000, 100); 
meteo_leaf.T = 20; 
x_values = meteo_leaf.Cs;
x_label = 'Cs, ppm CO_2';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% Temperature response

meteo_leaf.Q = 500;
meteo_leaf.Cs = 410; 
meteo_leaf.T = linspace(-10, 50, 10);
x_values = meteo_leaf.T;
x_label = 'T_{leaf}, ^oC';

leaf_out_vdt = biochemical(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);
leaf_out_md12 = biochemical_MD12(leafbio_adj,meteo_leaf,options_adj,constants_adj,fV);

plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)

%% 

function plot_leaf(x_values, x_label, leaf_out_vdt, leaf_out_md12, constants_adj)
    
    figure

    subplot(2,3,1)
    plot(x_values, leaf_out_vdt.A, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.A, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('A, \mumol CO_2 m-2 s-1')
    title('Net photosynthesis')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,2)
    plot(x_values, leaf_out_vdt.Ja, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.Ja, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('Ja, \mumol electrons m-2 s-1')
    title('Electron transport rate')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    
    subplot(2,3,3)
    
    num = constants_adj.rhoa./(constants_adj.Mair*1E-3);
    % gs = (rhoa./(Mair*1E-3)) / rcw  % mol H_2O m-2 s-1
    plot(x_values, num ./ leaf_out_vdt.rcw, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, num./leaf_out_md12.rcw, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    % ylabel('rcw, s m-1')
    ylabel('gsw, mol H_2O m-2 s-1')
    title('Stomatal conductance')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,4)
    plot(x_values, leaf_out_vdt.ps, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.ps, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('ps, \Phi_{PSII}')
    title('PSII photochemcal quantum yield')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,5)
    plot(x_values, leaf_out_vdt.fs, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.fs, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('fs, \Phi_F_s')
    title('Steady-state fluorescence')
    
    % legend('location', 'best', 'Interpreter','none')
    
    
    subplot(2,3,6)
    plot(x_values, leaf_out_vdt.eta, 'DisplayName', ['Fluorescence_model == 0', ...
        newline, 'Van der Tol et al., 2014'])
    hold on
    plot(x_values, leaf_out_md12.eta, 'DisplayName', ['Fluorescence_model == 1', ...
        newline, 'Magniani'])
    
    xlabel(x_label)
    ylabel('eta, F_s / F_o')
    title('eta, F_s / F_o')
    
    legend('location', 'best', 'Interpreter','none')
    
    set(findall(gcf,'-property','FontSize'),'FontSize', 14)
end

