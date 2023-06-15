function [Output_dir, f, fnames] = create_output_files_binary(parameter_file, F, path_of_code,input_path,spectral,options)
%% Create Output dir
string          = clock;
simulation_name = char(F(1).FileName);
outdir_name     = sprintf('%s_%4.0f-%02.0f-%02.0f-%02.0f%02.0f', simulation_name, string(1:5));
Output_dir = [fullfile('output', outdir_name) filesep];

warning('off','MATLAB:DELETE:FileNotFound')
if any(~exist(Output_dir,'dir'))
    mkdir(Output_dir)
    mkdir([Output_dir,'Parameters' filesep])
end

%% Log File
for i = 1:length(parameter_file{1})
    copy_name = [strrep(parameter_file{1}{i}, '.csv', '') '_' outdir_name '.csv'];
    copyfile([input_path parameter_file{1}{i}],[Output_dir,'Parameters/', copy_name],'f')
end
fidpath          = fopen([Output_dir,'Parameters/SCOPEversion.txt'],'w');      % complete path of the SCOPE code
fprintf(fidpath,'%s', path_of_code);

%% Filenames, will become .csv if options is on

fnames.pars_file               = fullfile(Output_dir,'pars_and_input_short.bin');
fnames.apar_file               = fullfile(Output_dir,'aPAR.bin');
fnames.veg_file                = fullfile(Output_dir,'vegetation.bin');
fnames.flu_file                = fullfile(Output_dir,'fluxes.bin');
fnames.rad_file                = fullfile(Output_dir,'radiation.bin');
if options.calc_fluor
    fnames.fluor_file              = fullfile(Output_dir,'fluorescence_scalars.bin');
    fnames.fluor_spectrum_file     = fullfile(Output_dir,'fluorescence.bin');
    fnames.sigmaF_file             = fullfile(Output_dir,'sigmaF.bin');
    fnames.fhemis_file             = fullfile(Output_dir,'fluorescence_hemis.bin');
    fnames.fRC_file                = fullfile(Output_dir,'fluorescence_ReabsCorr.bin');
    fnames.fRCL_file               = fullfile(Output_dir,'fluorescence_AllLeaves.bin'); 
    fnames.Lo2_file                = fullfile(Output_dir,'Lo_spectrum_inclF.bin');
    fnames.rapp_file               = fullfile(Output_dir,'apparent_reflectance.bin');
end
if options.save_spectral
    fnames.r_file                  = fullfile(Output_dir,'reflectance.bin');
    fnames.rsd_file                = fullfile(Output_dir,'rsd.bin');
    fnames.rdd_file                = fullfile(Output_dir,'rdd.bin');
    fnames.rso_file                = fullfile(Output_dir,'rso.bin');
    fnames.rdo_file                = fullfile(Output_dir,'rdo.bin');
    fnames.Eout_file               = fullfile(Output_dir,'Eout_spectrum.bin');
    fnames.Lo_file                 = fullfile(Output_dir,'Lo_spectrum.bin');
    fnames.Esun_file               = fullfile(Output_dir,'Esun.bin');
    fnames.Esky_file               = fullfile(Output_dir,'Esky.bin');
end
fnames.resist_file             = fullfile(Output_dir,'resistances.bin');

%% Open files for writing
f = structfun(@(x) fopen(x, 'w'), fnames, 'UniformOutput',false);

%% write wl
wlS = spectral.wlS; %#ok<*NASGU>
wlF = spectral.wlF;

save([Output_dir 'wlS.txt'], 'wlS', '-ascii');
save([Output_dir 'wlF.txt'], 'wlF', '-ascii');
end
