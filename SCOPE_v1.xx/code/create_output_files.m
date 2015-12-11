%% Create DATA files
% author J.timmermans
% last modified      4 Aug 2008: Added the creation of log file (file with input parameters)
%                    4 Aug 2008: j.timmermans: included variable output directories
%                   31 Jul 2008: (CvdT) added layer_pn.dat
%                   19 Sep 2008: (CvdT) added spectrum.dat
%                   16 Apr 2009: (CvdT) added layer_rn.dat
%                   18 Nov 2013: (CvdT) several updates.

%% Create Output dir

% a=mfilename('fullpath');
% c=strfind(a,filesep);
% pathsalida=a(1:c(end-1)-1);
% string          = datestr(now,30);
% 
% Outdir_Name     = char(F(1).FileName);
% Output_dir      = fullfile(pathsalida,'output',string);

string          = clock;

Outdir_Name     = char(F(1).FileName);
Output_dir      = sprintf(['../output/',Outdir_Name,'_%4.0f-%02.0f-%02.0f-%02.0f%02.0f/'],[string(1) string(2) string(3) string(4) string(5)]);

warning('off','MATLAB:DELETE:FileNotFound')
if any(~exist(Output_dir,'dir'))
    mkdir(Output_dir)
    mkdir([Output_dir,'Parameters/'])
    mkdir([Output_dir,'Directional/'])
    mkdir([Output_dir,'figures/'])
end


%% Log File
for i = 1:length(parameter_file)  
    copyfile(['../' parameter_file{i}],[Output_dir,'Parameters/', parameter_file{i}],'f')
end
fidpath          = fopen([Output_dir,'Parameters/SCOPEversion.txt'],'w');      % complete path of the SCOPE code
fprintf(fidpath,'%s', path_of_code);
%copyfile(['../' parameter_file],[Output_dir,'Parameters/', parameter_file ],'f')

%% Normal Output
fidf            = fopen([Output_dir,'fluxes.dat'],'w');       %     fluxes
fidt            = fopen([Output_dir,'surftemp.dat'],'w');         % surftemp
fidra           = fopen([Output_dir,'aerodyn.dat'],'w');          % aerodyn
fidr            = fopen([Output_dir,'radiation.dat'],'w');        % radiation
fidwl           = fopen([Output_dir,'wl.dat'],'w');               % wavelength
fidsi           = fopen([Output_dir,'irradiance_spectra.dat'],'w');     % Fluorescence
fidfho          = fopen([Output_dir,'spectrum_hemis_optical.dat'],'w');   % spectrum hemispherical
fidfoo          = fopen([Output_dir,'spectrum_obsdir_optical.dat'],'w');  % spectrum observation direction
fidref          = fopen([Output_dir,'reflectance.dat'],'w');      % reflectance spectrum

if options.calc_ebal
    fidto           = fopen([Output_dir,'spectrum_obsdir_BlackBody.dat'],'w');  % spectrum observation direction
end

if ~(options.simulation==1)
    fidv           = fopen([Output_dir,'pars_and_input.dat'],'w');               % wavelength
    for j = 1:length(V)
        fprintf(fidv,'%s\t',V(j).Name);
    end
    fprintf(fidv,'\r');
end

if ~(options.simulation==1)
    fidvs         = fopen([Output_dir,'pars_and_input_short.dat'],'a');
    for j = find(vmax>1)
        fprintf(fidvs,'%s\t',V(vmax>1).Name);
    end
    fprintf(fidvs,' \r');
end

%% Optional Output
if options.calc_vert_profiles
    fidgp       = fopen([Output_dir,'gap.dat'],'w');              % gap
    fidtc       = fopen([Output_dir,'leaftemp.dat'],'w');         % leaftemp
    fidhl       = fopen([Output_dir,'layer_H.dat'],'w');          % vertical profile
    fidlel      = fopen([Output_dir,'layer_lE.dat'],'w');         % latent heat
    fidal       = fopen([Output_dir,'layer_A.dat'],'w');          %
    fidpl       = fopen([Output_dir,'layer_aPAR.dat'],'w');       %
    fidplC      = fopen([Output_dir,'layer_aPAR_Cab.dat'],'w');       %
    fidrn       = fopen([Output_dir,'layer_Rn.dat'],'w');         %
    if options.calc_fluor
        fidfll = fopen([Output_dir,'layer_fluorescence.dat'],'w');   
        fidfllem = fopen([Output_dir,'layer_fluorescenceEm.dat'],'w'); 
        fidNPQ = fopen([Output_dir,'layer_NPQ.dat'],'w');
    end
    
else
    delete([Output_dir,'../output/leaftemp.dat'])
    delete([Output_dir,'../output/layer_H.dat'])
    delete([Output_dir,'../output/layer_lE.dat'])
    delete([Output_dir,'../output/layer_A.dat'])
    delete([Output_dir,'../output/layer_aPAR.dat'])
    delete([Output_dir,'../output/layer_Rn.dat'])
end

if options.calc_fluor
    fidfl       = fopen([Output_dir,'fluorescence.dat'],'w');     % Fluorescence
    fidfl1       = fopen([Output_dir,'fluorescencePSI.dat'],'w');     % Fluorescence
    fidfl2       = fopen([Output_dir,'fluorescencePSII.dat'],'w');     % Fluorescence
    fidflh      = fopen([Output_dir,'fluorescence_hemis.dat'],'w');     % Fluorescence
else
    delete([Output_dir,'fluorescence.dat'])
end

if options.calc_directional
    delete([Output_dir,'BRDF/*.dat'])
end

if options.calc_planck && options.calc_ebal
    fidplancko      = fopen([Output_dir,'spectrum_obsdir_thermal.dat'],'w');  % spectrum observation direction
    fidplanckh      = fopen([Output_dir,'spectrum_hemis_thermal.dat'],'w');  % spectrum hemispherically integrated
end

%% write headers
if options.save_headers
    fprintf(fidf,'timestep counter year t Rntot lEtot Htot Rnctot lEctot Hctot Actot Rnstot lEstot Hstot Gtot Resp aPAR aPAR_Cab faPAR aPAR_energyunits');
    if options.calc_fluor
        fprintf(fidf,' fluortot fluor_yield');
    end
    fprintf(fidf,'\r');
    fprintf(fidf,'""        ""      ""  JulianDay  Wm-2   Wm-2 Wm-2 Wm-2 Wm-2   Wm-2 umolm-2s-1 Wm-2 Wm-2 Wm-2 Wm-2 umolm-2s-1 umolm-2s-1 umolumol-1 Wm-2 ');
    if options.calc_fluor
        fprintf(fidf,' W m-2 WW^{-1}');
    end
    fprintf(fidf,'\r');
    
    fprintf(fidt,'timestep year t Ta Tss(1) Tss(2) Tcave Tsave \r');
    fprintf(fidt,'""        ""  JulianDay  ^oC ^oC ^oC ^oC ^oC \r');

    fprintf(fidra, 'raa rawc raws ustar \r');
    fprintf(fidra, 'sm-1 sm-1 sm-1 ms-1 \r');

    fprintf(fidr, 'timestep year t HemisOutShort HemisOutLong HemisOutTot  \r');
    fprintf(fidr, '""  ""  JulianDay  Wm-2 Wm-2 Wm-2\r');

    fprintf(fidfho, 'hemispherically integrated radiation spectrum \r');
    fprintf(fidfho, 'W m-2 um-1 \r');

    fprintf(fidfoo, 'reflectance spectrum in observation direction \r');
    fprintf(fidfoo, 'W m-2 sr-1 um-1 \r');

    
    if options.calc_ebal
        fprintf(fidto, 'thermal BlackBody emission spectrum in observation direction \r');
        fprintf(fidto, 'W m-2 sr-1 um-1 \r');
        if options.calc_planck
            fprintf(fidplancko, 'thermal emission spectrum in observation direction \r');
            fprintf(fidplancko, 'W m-2 sr-1 um-1 \r');
            
            fprintf(fidplanckh, 'thermal emission spectrum in hemispherical direction \r');
            fprintf(fidplanckh, 'W m-2 sr-1 um-1 \r');
        end
    end
    
    fprintf(fidwl, 'wavelengths of the spectral output files \r');
    fprintf(fidwl, 'um \r');
    
    fprintf(fidsi, 'irradiance \r');
    fprintf(fidsi, 'W m-2 um-1\r');
    
    fprintf(fidref, 'reflectance \r');
    fprintf(fidref, 'fraction of radiation in observation direction *pi / irradiance \r');
    
    if options.calc_fluor
        fprintf(fidfl, 'fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution \r');
        fprintf(fidfl, 'W m-2 um-1 sr-1\r');
        fprintf(fidfl1, 'fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution, for PSI only \r');
        fprintf(fidfl1, 'W m-2 um-1 sr-1\r');
        fprintf(fidfl2, 'fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution, for PSII only \r');
        fprintf(fidfl2, 'W m-2 um-1 sr-1\r');
        fprintf(fidflh, 'hemispherically integrated fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution \r');
        fprintf(fidflh, 'W m-2 um-1 \r');
        
    end
    if options.calc_vert_profiles
        fprintf(fidgp, 'Fraction leaves in the sun, fraction of observed, fraction of observed&visible per layer \r');
        fprintf(fidgp, ' rows: simulations or time steps, columns: layer numbers  \r');
        fprintf(fidtc, 'leaf temperature of sunlit leaves, shaded leaves, and weighted average leaf temperature per layer \r');
        fprintf(fidtc, '^oC ^oC ^oC \r');
        fprintf(fidhl, 'sensible heat flux per layer \r');
        fprintf(fidhl, 'Wm-2\r');
        fprintf(fidlel, 'latent heat flux per layer \r');
        fprintf(fidlel, 'Wm-2\r');
        fprintf(fidal, 'photosynthesis per layer\r');
        fprintf(fidal, 'umol-2s-1\r');
        fprintf(fidpl, 'aPAR per leaf layer \r');
        fprintf(fidpl, 'umol-2s-1 \r');
        fprintf(fidplC, 'aPAR by Cab per leaf layer \r');
        fprintf(fidplC, 'umol-2s-1 \r');
        fprintf(fidrn,'net radiation per leaf layer \r');
        fprintf(fidrn,'Wm-2\r');
        if options.calc_fluor
            fprintf(fidfll, 'upward fluorescence per layer\r');
            fprintf(fidfll, 'W m^{-2}\r');

            fprintf(fidNPQ, 'average NPQ = 1-(fm-fo)/(fm0-fo0), per layer \r');
            fprintf(fidNPQ, '\r');
        end
    end
end
%%
fprintf(fidwl,'%9.5f ',spectral.wlS);
warning('on','MATLAB:DELETE:FileNotFound')
fclose all;