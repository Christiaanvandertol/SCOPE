%% SCOPE.m (script)

%     SCOPE is a coupled radiative transfer and energy balance model.
%     Option 'lite' runs a computationally lighter variation of the model,
%     with the net radiation and leaf temperatures of leaf
%     inclination classes are averaged. SCOPE_lite is developed by C. van
%     der Tol of the Ts_sunlitUniversity of Twente, under subcontract of Magellium,
%     funded by the Europan Space Agency under contract FLEXL2-PFT-CCN2
%
%     Copyright (C) 2021  Christiaan van der Tol
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

clear all %#ok<CLALL>
restoredefaultpath
addpath src/RTMs
addpath src/supporting
addpath src/fluxes
addpath src/IO 

%% 1. define constants
constants = define_constants;

%% 2. paths
path_input      = 'input/';          % path of all inputs
path_of_code    = cd;

%% 3. simulation options
fid = fopen('set_parameter_filenames.csv','r');
parameter_file = textscan(fid,'%s','Delimiter', ',');
fclose(fid);

fid = fopen([path_input parameter_file{1}{1}],'r');
Ni = textscan(fid,'%d%s','Delimiter',',');%,'Whitespace','');
fclose(fid);
N = double(Ni{1});

options.lite                = N(1);    % lite version
options.calc_fluor          = N(2);    % calculate chlorophyll fluorescence in observation direction
options.calc_planck         = N(3);    % calculate spectrum of thermal radiation
options.calc_xanthophyllabs = N(4);    % include simulation of reflectance dependence on de-epoxydation state
options.soilspectrum        = N(5);    % 0: use soil reflectance from file; 1: calculate soil reflectance with BSM
options.Fluorescence_model  = N(6);     %0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
options.apply_T_corr        = N(7);     % correct Vcmax and rate constants for temperature in biochemical.m
options.verify              = N(8);
options.saveCSV             = N(9);
options.mSCOPE              = N(10);
options.simulation          = N(11);    % 0: individual runs (specify all input in a this file)
% 1: time series (uses text files with meteo input as time series)
% 2: Lookup-Table (specify the values to be included)
% 3: Lookup-Table with random input (specify the ranges of values)
options.calc_directional     = N(12);    % 0: calculate full BRDF (many angles)
options.calc_vert_profiles   = N(13);
options.soil_heat_method     = N(14);  % 0 - GAM=Soil_Inertia0(lambdas), 1 - GAM=Soil_Inertia1(SMC), 2 - G=0.35*Rn (always in no TS)
options.calc_rss_rbs          = N(15);  % 0 - fixed, 1 calc
options.MoninObukhov        = N(16);

if options.simulation>2 || options.simulation<0, fprintf('\n simulation option should be between 0 and 2 \r'); return, end
options.Cca_function_of_Cab = 0; % this will change to 1 if Cca is not provided in the input.

switch options.lite
    case 0, integr = 'angles_and_layers';
    otherwise, integr = 'layers';
end

%% 3. file names
f_names = {'Simulation_Name','soil_file','optipar_file','atmos_file', 'Dataset_dir',...
    'meteo_ec_csv', 'vegetation_retrieved_csv', 'LIDF_file', 'verification_dir', ...
    'mSCOPE_csv', 'nly'};  % must be in this order
cols = {'t', 'year', 'Rin','Rli', 'p','Ta','ea','u','RH', 'VPD', 'tts','tto', 'psi' ...  % expected from EC file as well as ('Ca','SMC')
    'Cab','Cca','Cdm','Cw','Cs','Cant','N'...  % leaf
    'SMC','BSMBrightness', 'BSMlat', 'BSMlon',...  % soil
    'LAI', 'hc', 'LIDFa', 'LIDFb',...  % canopy
    'z','Ca', ...  % meteo
    'Vcmo', 'm',...  % biochemistry;
    'atmos_names' 
    };

fnc = [f_names, cols];
F = struct('FileID', fnc);

fid = fopen([path_input parameter_file{1}{2}],'r');
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(line)
        charline = char(line);
        if ~(charline(1) == '%')
            X = textscan(line,'%s%s','Delimiter', ',', 'Whitespace','\t');
            x = X{1}; y = X{2};
            k = find(strcmp(fnc,x{:}));
            if ~isempty(k) && ~isempty(y)
                F(k).FileName = y{:};
            end
        end
    end
end
fclose(fid);

%% 4. input data
k = 1;
fid = fopen([path_input parameter_file{1}{3}], 'r');
clear('X')
while ~feof(fid)
    line    = fgetl(fid);
    y       = textscan(line,'%s', 'Delimiter', ',', 'TreatAsEmpty',  ' ');
    varnames(k)  = y{1}(1); %#ok<SAGROW>
    X(k).Val   = str2double(y{:});
    k       = k+1;
end
fclose(fid);
V                           = assignvarnames();

for i = 1:length(V)
    j = find(strcmp(varnames,V(i).Name));
    if isempty(j)
        if i==2
            fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
            fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
            options.Cca_function_of_Cab = 1;
        else
            if ~(options.simulation==1) && (i==30 || i==32)
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data...');
                fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum as it is');
            else
                if (options.simulation == 1 || (~options.simulation && (i<46 || i>50 )))
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input data');
                    if (options.simulation ==1)% && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                        fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',F(5).FileName,'"');
                    else
                        if (i== 24 || i==25)
                            fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                            options.calc_zo = 1;
                        else
                            if (i>38 && i<44)
                                fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                options.calc_zo = 0;
                            else
                                if ~((options.simulation ==1 && (i==30 ||i==32)))
                                    fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                    return
                                elseif (options.simulation ==1 && (i==30 ||i==32))
                                    fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',F(5).FileName, '"');
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        k = find(~isnan(X(j).Val));
        if ~isempty(k)
            V(i).Val = X(j).Val(k);
        else
            V(i).Val            = -999;
        end
    end
end

%% 6. Load spectral data for leaf and soil
load([path_input,'fluspect_parameters/', F(3).FileName]);
if options.soilspectrum ==0
    rsfile = load([path_input,'soil_spectra/', F(2).FileName]);        % file with soil reflectance spectra
end

%% 8. Define canopy structure and other 'fixed' parameters
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row
soilemp.SMC     = 25;        % empirical parameter (fixed) for BSM
soilemp.film    = 0.015;     % empirical parameter (fixed) for BMS
LIDF_file = F(8).FileName;
if  ~isempty(LIDF_file)
    canopy.lidf     = dlmread([path_input,'leafangles/',LIDF_file],'',3,0);
end

%% 10. Define spectral regions
[spectral] = define_bands;

%% 11. load time series data
if options.simulation == 1
    vi = ones(length(V),1);
    [soil,leafbio,canopy,meteo,angles,xyt]  = select_input(V,vi,canopy,options,constants);
    [V, xyt, mly_ts, atmo_paths]  = load_timeseries(V, F, xyt, path_input);
else
    soil = struct;
end

%% 12. preparations
%% soil heat
if options.simulation==1
    if options.soil_heat_method<2
        if (isempty(meteo.Ta) || meteo.Ta<-273), meteo.Ta = 20; end
        soil.Tsold = meteo.Ta*ones(12,2);
    end
end
%% temperature sensitivity of photosynthesis parameters
leafbio.TDP = define_temp_response_biochem; % temperature response C3 and C4 according to CLM4 model

%% variables
nvars = length(V);
vmax = cellfun(@length, {V.Val})';
vmax(27,1) = 1; % these are Tparam and LIDFb
vi      = ones(nvars,1);
switch options.simulation
    case 0, telmax = max(vmax);  [xyt.t,xyt.year]= deal(zeros(telmax,1));
    case 1, telmax = size(xyt.t,1);
    case 2, telmax  = prod(double(vmax)); [xyt.t,xyt.year]= deal(zeros(telmax,1));
end
% [rad,thermal,fluxes] = initialize_output_structures(spectral);

if options.calc_directional
    anglesfile          = load([path_input,'directional/brdf_angles2.dat']); %     Multiple observation angles in case of BRDF calculation
    directional.tto     = anglesfile(:,1);              % [deg]             Observation zenith Angles for calcbrdf
    directional.psi     = anglesfile(:,2);              % [deg]             Observation zenith Angles for calcbrdf
    directional.noa     = length(directional.tto);      %                   Number of Observation Angles
else
    directional = NaN;
end

%% irradiance
atmfile = fullfile(path_input, 'radiationdata', F(4).FileName);
if options.simulation == 1 && ~isempty(atmo_paths)
    atmfile = atmo_paths{1};
end
atmo = load_atmo(atmfile, spectral.SCOPEspec);

%% 13. create output files
[Output_dir, f, fnames] = create_output_files_binary(parameter_file, F, path_of_code, path_input, spectral,options);

%% 14. Run the models
fprintf('\n The calculations start now \r')
calculate = 1;
tic
        
for k = 1:telmax
    if options.simulation == 1, vi(vmax>1) = k; end
    if options.simulation == 0, vi(vmax==telmax) = k; end
    [soil,leafbio,canopy,meteo,angles,xyt] = select_input(V,vi,canopy,options,constants,xyt,soil,leafbio);
    canopy.nlayers  = ceil(10*canopy.LAI);
    canopy.nlayers = max(2, canopy.nlayers);  % patch for LAI < 0.1
    nl              = canopy.nlayers;
	x        = (-1/nl : -1/nl : -1)';         % a column vector
    canopy.xl       = [0; x];                 % add top level
  % canopy.xl(1:end-1) = canopy.xl(1:end-1)+canopy.xl(1:end-1)-1/(2*nl); % middle of the thin layer
    
    if options.simulation ~=1
        fprintf('simulation %i ', k );
        fprintf('of %i \n', telmax);
    else
        calculate = ~isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli);
        fprintf('time = %s: %i / %i\n', datestr(xyt.t(k)), k, telmax)
        if isnan(meteo.p*meteo.Ta*meteo.ea*meteo.u.*meteo.Rin.*meteo.Rli)
            warning('run is invalid: there is NaN somewhere in meteo input [p, Ta, ea, u, Rin, Rli]')
        end
    end
    
    if calculate
        if isempty(LIDF_file)
            canopy.lidf     = leafangles(canopy.LIDFa,canopy.LIDFb);    % This is 'ladgen' in the original SAIL model,
        end
        
        %% leaf radiative transfer model FLUSPECT
        leafbio.emis        = 1-leafbio.rho_thermal-leafbio.tau_thermal;
        leafbio.V2Z         = 0;
        
        if options.simulation == 1 && ~isempty(fieldnames(mly_ts))  % means that options.simulation == 1 
           mly.nly    = mly_ts.nly;
           mly.pLAI   = mly_ts.pLAI(k, :);
           mly.totLAI = sum(mly.pLAI);
           mly.pCab   = mly_ts.pCab(k, :);
           mly.pCca   = mly_ts.pCca(k, :);
           mly.pCdm   = mly_ts.pCw(k, :);
           mly.pCw    = mly_ts.pCw(k, :);
           mly.pCs    = mly_ts.pCs(k, :);
           mly.pN     = mly_ts.pN(k, :);
        elseif k == 1 && options.mSCOPE
           mly = input_mSCOPE(fullfile('input', 'mSCOPE.csv'));
        else
           if options.mSCOPE
                warning('I do not know how to use mSCOPE layers with multiple but non time series runs, so I will not use it')
           end
           mly.nly      = 1;
           mly.pLAI     = canopy.LAI;
           mly.totLAI   = canopy.LAI;
           mly.pCab     = leafbio.Cab;
           mly.pCca     = leafbio.Cca;
           mly.pCdm     = leafbio.Cdm;
           mly.pCw      = leafbio.Cw;
           mly.pCs      = leafbio.Cs;
           mly.pN       = leafbio.N;
        end
        
        if options.simulation == 1 && ~isempty(atmo_paths) && k > 1
            atmfile_k = atmo_paths{k};
            if ~strcmp(atmfile_k, atmo_paths{k-1})
                atmo = load_atmo(atmfile_k, spectral.SCOPEspec);
            end
        end
        
        leafopt = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl); 
        leafopt.refl(:, spectral.IwlT) = leafbio.rho_thermal;               
        leafopt.tran(:, spectral.IwlT) = leafbio.tau_thermal;
        
        if options.calc_xanthophyllabs
            leafbio.V2Z     = 1;
            leafoptZ        = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl); 
            leafopt.reflZ   = leafopt.refl;
            leafopt.tranZ   = leafopt.tran;
            leafopt.reflZ(:, spectral.IwlP) = leafoptZ.refl(:, spectral.IwlP);
            leafopt.tranZ(:, spectral.IwlP) = leafoptZ.tran(:, spectral.IwlP);
        end
        
        %% soil reflectance model BSM
        if options.soilspectrum == 0
            soil.refl       = rsfile(:,soil.spectrum+1);
        else
            soil.refl       = BSM(soil,optipar,soilemp);
        end
        soil.refl(spectral.IwlT) = soil.rs_thermal;
        
        %% four stream canopy radiative transfer model for incident radiation
        [rad,gap,profiles]       = RTMo(spectral,atmo,soil,leafopt,canopy,angles,constants,meteo,options);
        
        %% energy balance
        [iter,rad,thermal,soil,bcu,bch,fluxes,resistance,meteo]             ...
            = ebal(constants,options,rad,gap,  ...
            meteo,soil,canopy,leafbio, k, xyt,integr);
        
        %% fluorescence radiative transfer model
        if options.calc_fluor
            [rad]           = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles,bcu.eta,bch.eta);
        end
        
        %% radiative transfer model for PRI effects
        if options.calc_xanthophyllabs
            [rad] = RTMz(constants,spectral,rad,soil,leafopt,canopy,gap,angles,bcu.Kn,bch.Kn);
        end
                       
        rad  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh,1,spectral);
        if options.calc_planck
            rad = RTMt_planck(spectral,rad,soil,leafopt,canopy,gap,thermal.Tcu,thermal.Tch,thermal.Tsu,thermal.Tsh);
        end
        
        %% computation of data products
        % aPAR, LST, NPQ, ETR, photosynthesis, SIF-reabsorption correction
        % aPAR [umol m-2 s-1, total canopy and total chlorphyll]
 
        Ps = gap.Ps(1:nl);
        Ph = (1-Ps);
        
        canopy.LAIsunlit = canopy.LAI*mean(Ps);
        canopy.LAIshaded = canopy.LAI-canopy.LAIsunlit;
        
        canopy.Pnsun_Cab = canopy.LAI*meanleaf(canopy,rad.Pnu_Cab,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha_Cab = canopy.LAI*meanleaf(canopy,rad.Pnh_Cab,'layers',Ph); % net PAR Cab shaded leaves (photons)
        canopy.Pntot_Cab = canopy.Pnsun_Cab+canopy.Pnsha_Cab; % net PAR Cab leaves (photons)

        canopy.Pnsun_Car = canopy.LAI*meanleaf(canopy,rad.Pnu_Car,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha_Car = canopy.LAI*meanleaf(canopy,rad.Pnh_Car,'layers',Ph); % net PAR Cab shaded leaves (photons)
        canopy.Pntot_Car = canopy.Pnsun_Car+canopy.Pnsha_Car; % net PAR Cab leaves (photons)
        
        canopy.Pnsun = canopy.LAI*meanleaf(canopy,rad.Pnu,integr,Ps); % net PAR sunlit leaves (photons)
        canopy.Pnsha = canopy.LAI*meanleaf(canopy,rad.Pnh,'layers',Ph); % net PAR shaded leaves (photons)
        canopy.Pntot = canopy.Pnsun+canopy.Pnsha; % net PAR leaves (photons)
        
        canopy.Rnsun_Cab = canopy.LAI*meanleaf(canopy,rad.Rnu_Cab,integr,Ps); % net PAR Cab sunlit leaves (radiance)
        canopy.Rnsha_Cab = canopy.LAI*meanleaf(canopy,rad.Rnh_Cab,'layers',Ph); % net PAR Cab sunlit leaves (radiance)
        canopy.Rntot_Cab = canopy.Rnsun_Cab+canopy.Rnsha_Cab; % net PAR Cab leaves (radiance)
        
        canopy.Rnsun_Car = canopy.LAI*meanleaf(canopy,rad.Rnu_Car,integr,Ps); % net PAR Cab sunlit leaves (radiance)
        canopy.Rnsha_Car = canopy.LAI*meanleaf(canopy,rad.Rnh_Car,'layers',Ph); % net PAR Cab sunlit leaves (radiance)
        canopy.Rntot_Car = canopy.Rnsun_Car+canopy.Rnsha_Car; % net PAR Cab leaves (radiance)
        
        canopy.Rnsun_PAR = canopy.LAI*meanleaf(canopy,rad.Rnu_PAR,integr,Ps); % net PAR sunlit leaves (radiance)
        canopy.Rnsha_PAR = canopy.LAI*meanleaf(canopy,rad.Rnh_PAR,'layers',Ph); % net PAR sunlit leaves (radiance)
        canopy.Rntot_PAR = canopy.Rnsun_PAR+canopy.Rnsha_PAR; % net PAR leaves (radiance)
        
        % LST [K] (directional, but assuming black-body surface!)
        canopy.LST      = (pi*(rad.Lot+rad.Lote)./(constants.sigmaSB*rad.canopyemis)).^0.25;
        canopy.emis     = rad.canopyemis;
        
        % photosynthesis [mumol m-2 s-1]
        canopy.A        = canopy.LAI*(meanleaf(canopy,bch.A,'layers',Ph)+meanleaf(canopy,bcu.A,integr,Ps)); % photosynthesis
        
        % electron transport rate [mumol m-2 s-1]
        canopy.Ja       = canopy.LAI*(meanleaf(canopy,bch.Ja,'layers',Ph)+meanleaf(canopy,bcu.Ja,integr,Ps)); % electron transport
        
        % non-photochemical quenching (energy) [W m-2]
        canopy.ENPQ     = canopy.LAI*(meanleaf(canopy,rad.Rnh_Cab.*bch.Phi_N,'layers',Ph)+meanleaf(canopy,rad.Rnu_Cab.*bcu.Phi_N,integr,Ps)); % NPQ energy;
        canopy.PNPQ     = canopy.LAI*(meanleaf(canopy,rad.Pnh_Cab.*bch.Phi_N,'layers',Ph)+meanleaf(canopy,rad.Pnu_Cab.*bcu.Phi_N,integr,Ps)); % NPQ energy;
        
        % computation of re-absorption corrected fluorescence
        % Yang and Van der Tol (2019); Van der Tol et al. (2019)
        %aPAR_Cab_eta    = canopy.LAI*(meanleaf(canopy,bch.eta .* rad.Rnh_Cab,'layers',Ph)+meanleaf(canopy,bcu.eta .* rad.Rnu_Cab,integr,Ps)); %
        aPAR_Cab_eta    = canopy.LAI*(meanleaf(canopy,bch.eta .* rad.Pnh_Cab,'layers',Ph)+meanleaf(canopy,bcu.eta .* rad.Pnu_Cab,integr,Ps)); %
        if options.calc_fluor
            ep              = constants.A*ephoton(spectral.wlF'*1E-9,constants);
            rad.PoutFrc     = leafbio.fqe*aPAR_Cab_eta;
            rad.EoutFrc_    = 1E-3*ep.*(rad.PoutFrc*optipar.phi(spectral.IwlF)); %1E-6: umol2mol, 1E3: nm-1 to um-1
            rad.EoutFrc     = 1E-3*Sint(rad.EoutFrc_,spectral.wlF);
            sigmaF          = pi*rad.LoF_./rad.EoutFrc_;
            rad.sigmaF      = interp1(spectral.wlF(1:4:end),sigmaF(1:4:end),spectral.wlF);
            canopy.fqe      = rad.PoutFrc./canopy.Pntot_Cab;
        end
        
        rad.Lotot_      = rad.Lo_+rad.Lot_;
        rad.Eout_       = rad.Eout_+rad.Eoutte_;
        if options.calc_fluor
            rad.Lototf_     = rad.Lotot_;
            rad.Lototf_(spectral.IwlF') = rad.Lototf_(spectral.IwlF)+rad.LoF_;            
            rad.reflapp = rad.refl;
            rad.reflapp(spectral.IwlF) =pi*rad.Lototf_(spectral.IwlF)./(rad.Esun_(spectral.IwlF)+rad.Esky_(spectral.IwlF));
        end
        
        if options.calc_directional
            directional = calc_brdf(constants,options,directional,spectral,angles,atmo,soil,leafopt,canopy,meteo,thermal,bcu,bch);
            savebrdfoutput(options,directional,angles,spectral,Output_dir)
        end
        
        rad.Lo = 0.001 * Sint(rad.Lo_(spectral.IwlP),spectral.wlP);
        %% write output
        n_col = output_data_binary(f, k, xyt, rad, canopy, V, vi, vmax, options, fluxes, meteo, iter,resistance);
        
        %% update input
        if options.simulation==2 && telmax>1, vi  = count_k(nvars,vi,vmax,1); end
    end
end
toc
if options.saveCSV
    bin_to_csv(fnames, V, vmax, n_col, telmax)
end
fclose('all');

if options.verify
    output_verification_csv(Output_dir, F(9).FileName)
end

