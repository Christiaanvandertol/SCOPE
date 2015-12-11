function [const]=define_constants()

    const.A         = 6.02214E23; % [mol-1]       Constant of Avogadro
    const.h         = 6.6262E-34; % [J s]         Planck's constant
    const.c         = 299792458;  % [m s-1]       Speed of light
    const.cp        = 1004;       % [J kg-1 K-1]  Specific heat of dry air
    const.R         = 8.314;      % [J mol-1K-1]  Molar gas constant
    const.rhoa      = 1.2047;     % [kg m-3]      Specific mass of air
    const.g         = 9.81;       % [m s-2]       Gravity acceleration
    const.kappa     = 0.4;        % []            Von Karman constant
    const.MH2O      = 18;         % [g mol-1]     Molecular mass of water
    const.Mair      = 28.96;      % [g mol-1]     Molecular mass of dry air
    const.MCO2      = 44;         % [g mol-1]     Molecular mass of carbon dioxide
    const.sigmaSB   = 5.67E-8;    % [W m-2 K-4]   Stefan Boltzman constant  
    const.deg2rad   = pi/180;     % [rad]         Conversion from deg to rad
    const.C2K       = 273.15;     % [K]           Melting point of water
    
end