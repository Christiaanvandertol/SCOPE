%% SCOPE.m (script)

%     SCOPE is a coupled radiative transfer and energy balance model
%     Copyright (C) 2015  Christiaan van der Tol
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%
clc
%clear all

%% 0. globals
global constants

%% 1. define constants
[constants] = io.define_constants();

%% 2. simulation options
path_of_code                = cd;
run ../set_parameter_filenames; 
% parameter_file = {'input_data.xlsx'};  % for .exe

if length(parameter_file)>1, useXLSX = 0; else useXLSX = 1; end

if ~useXLSX
    run(['../' parameter_file{1}])
    
    options.calc_ebal           = N(1);    % calculate the energy balance (default). If 0, then only SAIL is executed!
    options.calc_vert_profiles  = N(2);    % calculate vertical profiles of fluxes
    options.calc_fluor          = N(3);    % calculate chlorophyll fluorescence in observation direction
    options.calc_planck         = N(4);    % calculate spectrum of thermal radiation
    options.calc_directional    = N(5);    % calculate BRDF and directional temperature
    options.calc_xanthophyllabs = N(6);    % include simulation of reflectance dependence on de-epoxydation state
    options.calc_PSI            = N(7);    % 0: optipar 2017 file with only one fluorescence spectrum vs 1: Franck et al spectra for PSI and PSII
    options.rt_thermal          = N(8);    % 1: use given values under 10 (default). 2: use values from fluspect and soil at 2400 nm for the TIR range
    options.calc_zo             = N(9);
    options.soilspectrum        = N(10);    %0: use soil spectrum from a file, 1: simulate soil spectrum with the BSM model
    options.soil_heat_method    = N(11);    % 0: calculated from specific heat and conductivity (default), 1: empiricaly calibrated, 2: G as constant fraction of soil net radiation
    options.Fluorescence_model  = N(12);     %0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
    options.calc_rss_rbs        = N(13);    % 0: calculated from specific heat and conductivity (default), 1: empiricaly calibrated, 2: G as constant fraction of soil net radiation
    options.apply_T_corr        = N(14);    % correct Vcmax and rate constants for temperature in biochemical.m
    options.verify              = N(15);
    options.save_headers        = N(16);    % write headers in output files
    options.makeplots           = N(17);
    options.simulation          = N(18);    % 0: individual runs (specify all input in a this file)
    % 1: time series (uses text files with meteo input as time series)
    % 2: Lookup-Table (specify the values to be included)
    % 3: Lookup-Table with random input (specify the ranges of values)
else
    options = io.readStructFromExcel(['../' char(parameter_file)], 'options', 3, 1);
end

if options.simulation>2 || options.simulation<0, fprintf('\n simulation option should be between 0 and 2 \r'); return, end

%% 3. file names
if ~useXLSX
    run(['../' parameter_file{2}])
else
    [dummy,X]                       = xlsread(['../' char(parameter_file)],'filenames');
    j = find(~strcmp(X(:,2),{''}));
    X = X(j,(1:end));
end

F = struct('FileID',{'Simulation_Name','soil_file','leaf_file','atmos_file'...
    'Dataset_dir','t_file','year_file','Rin_file','Rli_file'...
    ,'p_file','Ta_file','ea_file','u_file','CO2_file','z_file','tts_file'...
    ,'LAI_file','hc_file','SMC_file','Vcmax_file','Cab_file','LIDF_file'});
for i = 1:length(F)
    k = find(strcmp(F(i).FileID,strtok(X(:,1))));
    if ~isempty(k)
        F(i).FileName = strtok(X(k,2));
        %if i==4, F(i).FileName = strtok(X(k,2:end)); end
    end
end

%% 4. input data

if ~useXLSX
    X                           = textread(['../' parameter_file{3}],'%s'); %#ok<DTXTRD>
    N                           = str2double(X);
else
    [N,X]                       = xlsread(['../' char(parameter_file)],'inputdata', '');
    X                           = X(9:end,1);
end
V                           = io.assignvarnames();
options.Cca_function_of_Cab = 0;

for i = 1:length(V)
    j = find(strcmp(strtok(X(:,1)),V(i).Name));
    if ~useXLSX, cond = isnan(N(j+1)); else cond = sum(~isnan(N(j,:)))<1; end
    if isempty(j) || cond
        if i==2
            fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
            fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
            options.Cca_function_of_Cab = 1;
        else
            
            if ~(options.simulation==1) && (i==30 || i==32)
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
                fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum as it is');
            else
                if (options.simulation == 1 || (options.simulation~=1 && (i<46 || i>50)))
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet');
                    if (options.simulation ==1 && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                        fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',char(F(5).FileName),'"');
                    else
                        if (i== 24 || i==25)
                            fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                            options.calc_zo = 1;
                        else
                            if (i>38 && i<44)
                                fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                options.calc_zo = 0;
                            else
                                if ~(options.simulation ==1 && (i==30 ||i==32))
                                    fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                    return
                                else
                                    fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',char(F(5).FileName), '"');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if ~useXLSX
        j2 = []; j1 = j+1;
        while 1
            if isnan(N(j1)), break, end
            j2 = [j2; j1]; %#ok<AGROW>
            j1 = j1+1;
        end
        if isempty(j2)
            V(i).Val            = -999;
        else
            V(i).Val            = N(j2);
        end
        
        
    else
        if sum(~isnan(N(j,:)))<1
            V(i).Val            = -999;
        else
            V(i).Val            = N(j,~isnan(N(j,:)));
        end
    end
end

%% 5. Declare paths
path_input      = '../../data/input/';          % path of all inputs

%% 6. Numerical parameters (iteration stops etc)
iter.maxit           = 100;                          %                   maximum number of iterations
iter.maxEBer         = 1;                            %[W m-2]            maximum accepted error in energy bal.
iter.Wc              = 1;                         %                   Weight coefficient for iterative calculation of Tc

%% 7. Load spectral data for leaf and soil
%opticoef    = xlsread([path_input,'fluspect_parameters/',char(F(3).FileName)]);  % file with leaf spectral parameters
%xlsread([path_input,'fluspect_parameters/',char(F(3).FileName)]);  % file with leaf spectral parameters
load([path_input,'fluspect_parameters/',char(F(3).FileName)]);
rsfile      = load([path_input,'soil_spectrum/',char(F(2).FileName)]);        % file with soil reflectance spectra
% Optical coefficient data used by fluspect
% optipar.nr    = opticoef(:,2);
% optipar.Kab   = opticoef(:,3);
% optipar.Kca   = opticoef(:,4);
% optipar.Ks    = opticoef(:,5);
% optipar.Kw    = opticoef(:,6);
% optipar.Kdm   = opticoef(:,7);
% optipar.nw   = opticoef(:,8);
% optipar.phiI  = opticoef(:,9);
% optipar.phiII = opticoef(:,10);
% optipar.GSV1  = opticoef(:,11); 
% optipar.GSV2  = opticoef(:,12);
% optipar.GSV3  = opticoef(:,13);
% optipar.KcaV   = opticoef(:,14);
% optipar.KcaZ   = opticoef(:,15);

%% 8. Load directional data from a file
if options.calc_directional
    anglesfile          = load([path_input,'directional/brdf_angles2.dat']); %     Multiple observation angles in case of BRDF calculation
    directional.tto     = anglesfile(:,1);              % [deg]             Observation zenith Angles for calcbrdf
    directional.psi     = anglesfile(:,2);              % [deg]             Observation zenith Angles for calcbrdf
    directional.noa     = length(directional.tto);      %                   Number of Observation Angles
end

%% 9. Define canopy structure
canopy.nlayers  = 60;
nl              = canopy.nlayers;
canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
canopy.xl       = [0; canopy.x];                 % add top level
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row

%% 10. Define spectral regions
[spectral] = io.define_bands();

wlS  = spectral.wlS;    % SCOPE 1.40 definition
wlP  = spectral.wlP;    % PROSPECT (fluspect) range
wlT  = spectral.wlT;    % Thermal range
wlF  = spectral.wlF;    % Fluorescence range

I01  = find(wlS<min(wlF));   % zero-fill ranges for fluorescence
I02  = find(wlS>max(wlF));
N01  = length(I01);
N02  = length(I02);

nwlP = length(wlP);
nwlT = length(wlT);

nwlS = length(wlS);

spectral.IwlP = 1 : nwlP;
spectral.IwlT = nwlP+1 : nwlP+nwlT;
spectral.IwlF = (640:850)-399;

[rho,tau,rs] = deal(zeros(nwlP + nwlT,1));

%% 11. load time series data
if options.simulation == 1
    vi = ones(length(V),1);
    [soil,leafbio,canopy,meteo,angles,xyt]  = io.select_input(V,vi,canopy,options);
    [V,xyt,canopy]  = io.load_timeseries(V,leafbio,soil,canopy,meteo,constants,F,xyt,path_input,options);
else
    soil = struct;
end

%% 12. preparations
if options.simulation==1
    diff_tmin           =   abs(xyt.t-xyt.startDOY);
    diff_tmax           =   abs(xyt.t-xyt.endDOY);
    I_tmin              =   find(min(diff_tmin)==diff_tmin);
    I_tmax              =   find(min(diff_tmax)==diff_tmax);
    if options.soil_heat_method<2
        if (isempty(meteo.Ta) || meteo.Ta<-273), meteo.Ta = 20; end
        soil.Tsold = meteo.Ta*ones(12,2);
    end
end

nvars = length(V);
vmax = ones(nvars,1);
for i = 1:nvars
    vmax(i) = length(V(i).Val);
end
vmax([14,27],1) = 1; % these are Tparam and LIDFb
vi      = ones(nvars,1);
switch options.simulation
    case 0, telmax = max(vmax);  [xyt.t,xyt.year]= deal(zeros(telmax,1));
    case 1, telmax = size(xyt.t,1);
    case 2, telmax  = prod(double(vmax)); [xyt.t,xyt.year]= deal(zeros(telmax,1));
end
[rad,thermal,fluxes] = io.initialize_output_structures(spectral);
atmfile     = [path_input 'radiationdata/' char(F(4).FileName(1))];
atmo.M      = helpers.aggreg(atmfile,spectral.SCOPEspec);

%% 13. create output files
Output_dir = io.create_output_files(parameter_file, F, path_of_code, options, V, vmax, spectral);

%% 14. Run the model
fprintf('\n The calculations start now \r')
calculate = 1;

for k = 1:telmax
    
    if options.simulation == 1, vi(vmax>1) = k; end
    if options.simulation == 0, vi(vmax==telmax) = k; end
    [soil,leafbio,canopy,meteo,angles,xyt] = io.select_input(V,vi,canopy,options,xyt,soil);
    if options.simulation ~=1
        fprintf('simulation %i ', k );
        fprintf('of %i \n', telmax);
    else
        calculate = 0;
        if k>=I_tmin && k<=I_tmax
            quality_is_ok   = ~isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli);
            fprintf('time = %4.2f \n', xyt.t(k));
            if quality_is_ok
                calculate = 1;
            end
        end
    end
    
    if calculate
        
        iter.counter = 0;
        
        LIDF_file            = char(F(22).FileName);
        if  ~isempty(LIDF_file)
            canopy.lidf     = dlmread([path_input,'leafangles/',LIDF_file],'',3,0);
        else
            canopy.lidf     = equations.leafangles(canopy.LIDFa,canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
        end
        
        if options.calc_PSI
            fversion = @fluspect_B_CX;
        else
            fversion = @fluspect_B_CX_PSI_PSII_combined;
        end
        leafbio.V2Z = 0;
        leafopt  = fversion(spectral,leafbio,optipar);
        leafbio.V2Z = 1;
        leafoptZ = fversion(spectral,leafbio,optipar);
        
        IwlP     = spectral.IwlP;
        IwlT     = spectral.IwlT;
        
        rho(IwlP)  = leafopt.refl;
        tau(IwlP)  = leafopt.tran;
        rlast    = rho(nwlP);
        tlast    = tau(nwlP);
        
        if options.soilspectrum == 0
            rs(IwlP) = rsfile(:,soil.spectrum+1);
        else
            soilemp.SMC   = 25;        % empirical parameter (fixed)
            soilemp.film  = 0.015;     % empirical parameter (fixed)
            rs(IwlP) = BSM(soil,optipar,soilemp);
        end
        rslast   = rs(nwlP);
        
        switch options.rt_thermal
            case 0
                rho(IwlT) = ones(nwlT,1) * leafbio.rho_thermal;
                tau(IwlT) = ones(nwlT,1) * leafbio.tau_thermal;
                rs(IwlT)  = ones(nwlT,1) * soil.rs_thermal;
            case 1
                rho(IwlT) = ones(nwlT,1) * rlast;
                tau(IwlT) = ones(nwlT,1) * tlast;
                rs(IwlT)  = ones(nwlT,1) * rslast;
        end
        leafopt.refl = rho;     % extended wavelength ranges are stored in structures
        leafopt.tran = tau;
        
        reflZ = leafopt.refl;
        tranZ = leafopt.tran;
        reflZ(1:300) = leafoptZ.refl(1:300);
        tranZ(1:300) = leafoptZ.tran(1:300);
        leafopt.reflZ = reflZ;
        leafopt.tranZ = tranZ;
        
        soil.refl    = rs;
        
        soil.Ts     = meteo.Ta * ones(2,1);       % initial soil surface temperature
        
        if length(F(4).FileName)>1 && options.simulation==0
            atmfile     = [path_input 'radiationdata/' char(F(4).FileName(k))];
            atmo.M      = helpers.aggreg(atmfile,spectral.SCOPEspec);
        end
        atmo.Ta     = meteo.Ta;
        
        [rad,gap,profiles]   = RTMo(spectral,atmo,soil,leafopt,canopy,angles,meteo,rad,options);
        
        switch options.calc_ebal
            case 1
                [iter,fluxes,rad,thermal,profiles,soil]                          ...
                    = ebal(iter,options,spectral,rad,gap,                       ...
                    leafopt,angles,meteo,soil,canopy,leafbio,xyt,k,profiles);
                
                if options.calc_fluor
                    if options.calc_vert_profiles
                        [rad,profiles] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                    else
                        [rad] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                    end
                end
                if options.calc_xanthophyllabs
                    [rad] = RTMz(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                end
                
                if options.calc_planck
                    rad         = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,angles,thermal.Tcu,thermal.Tch,thermal.Ts(2),thermal.Ts(1),1);
                end
                
                if options.calc_directional
                    directional = calc_brdf(options,directional,spectral,angles,rad,atmo,soil,leafopt,canopy,meteo,profiles,thermal);
                end
                
            otherwise
                Fc              = (1-gap.Ps(1:end-1))'/nl;      %           Matrix containing values for Ps of canopy
                fluxes.aPAR     = canopy.LAI*(Fc*rad.Pnh        + equations.meanleaf(canopy,rad.Pnu    , 'angles_and_layers',gap.Ps));% net PAR leaves
                fluxes.aPAR_Cab = canopy.LAI*(Fc*rad.Pnh_Cab    + equations.meanleaf(canopy,rad.Pnu_Cab, 'angles_and_layers',gap.Ps));% net PAR leaves
                [fluxes.aPAR_Wm2,fluxes.aPAR_Cab_eta] = deal(canopy.LAI*(Fc*rad.Rnh_PAR    + equations.meanleaf(canopy,rad.Rnu_PAR, 'angles_and_layers',gap.Ps)));% net PAR leaves
                if options.calc_fluor
                    profiles.etah = ones(60,1);
                    profiles.etau = ones(13,36,60);
                    if options.calc_vert_profiles
                        [rad,profiles] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                    else
                        [rad] = RTMf(spectral,rad,soil,leafopt,canopy,gap,angles,profiles);
                    end
                end
        end
        if options.calc_fluor % total emitted fluorescence irradiance (excluding leaf and canopy re-absorption and scattering)
            if options.calc_PSI
                rad.Femtot = 1E3*(leafbio.fqe(2)* optipar.phiII(spectral.IwlF) * fluxes.aPAR_Cab_eta +leafbio.fqe(1)* optipar.phiI(spectral.IwlF)  * fluxes.aPAR_Cab);
            else
                rad.Femtot = 1E3*leafbio.fqe* optipar.phi(spectral.IwlF) * fluxes.aPAR_Cab_eta;
            end
        end    
        io.output_data(Output_dir, options, k, iter, xyt, fluxes, rad, thermal, gap, meteo, spectral, V, vi, vmax, profiles)
    end
    if options.simulation==2 && telmax>1, vi  = helpers.count(nvars,vi,vmax,1); end
end

if options.verify
    io.output_verification(Output_dir)
end

if options.makeplots
    plot.plots(Output_dir)
end

%% for Compiler
% catch ME
%     disp(['ERROR: ' ME.message])
% end
% fprintf('\nThe run is finished. Press any key to close the window')
% fprintf('\nIf no error message was produced navigate to ./SCOPE_v1.70/output to see the results')
% pause
