function Output_dir = create_output_files_csv(parameter_file, F, path_of_code, options, V, vmax, spectral)
    %% Create Output dir
    string          = clock;
    Outdir_Name     = char(F(1).FileName);
    Output_dir      = sprintf(['../output/',Outdir_Name,'_%4.0f-%02.0f-%02.0f-%02.0f%02.0f/'],[string(1) string(2) string(3) string(4) string(5)]);
    warning('off','MATLAB:DELETE:FileNotFound')
    if any(~exist(Output_dir,'dir'))
        mkdir(Output_dir)
        mkdir([Output_dir,'Parameters/'])
%         mkdir([Output_dir,'Directional/'])
%         mkdir([Output_dir,'figures/'])
    end


    %% Log File
    for i = 1:length(parameter_file)  
        copyfile(['../' parameter_file{i}],[Output_dir,'Parameters/', parameter_file{i}],'f')
    end
    fidpath          = fopen([Output_dir,'Parameters/SCOPEversion.txt'],'w');      % complete path of the SCOPE code
    fprintf(fidpath,'%s', path_of_code);
    %copyfile(['../' parameter_file],[Output_dir,'Parameters/', parameter_file ],'f')

    %% Normal Output
    fidf            = fopen([Output_dir,'fluxes.csv'],'w');       %     fluxes
    fidt            = fopen([Output_dir,'surftemp.csv'],'w');         % surftemp
    fidra           = fopen([Output_dir,'aerodyn.csv'],'w');          % aerodyn
    fidr            = fopen([Output_dir,'radiation.csv'],'w');        % radiation
    fidwl           = fopen([Output_dir,'wl.csv'],'w');               % wavelength
    fidsi           = fopen([Output_dir,'irradiance_spectra.csv'],'w');     % Fluorescence
    fidfho          = fopen([Output_dir,'spectrum_hemis_optical.csv'],'w');   % spectrum hemispherical
    fidfoo          = fopen([Output_dir,'spectrum_obsdir_optical.csv'],'w');  % spectrum observation direction
    fidref          = fopen([Output_dir,'reflectance.csv'],'w');      % reflectance spectrum
    fidp            = fopen([Output_dir,'BOC_irradiance.csv'],'w');

    if options.calc_ebal
        fidto           = fopen([Output_dir,'spectrum_obsdir_BlackBody.csv'],'w');  % spectrum observation direction
    end

    fidv           = fopen([Output_dir,'pars_and_input.csv'],'w');
    for j = 1:length(V)
        fprintf(fidv,'%s,',V(j).Name);
    end
    fprintf(fidv,'\n');

    fidvs         = fopen([Output_dir,'pars_and_input_short.csv'],'a');
    for j = find(vmax>1)
        fprintf(fidvs,'%s,',V(vmax>1).Name);
    end
    fprintf(fidvs,' \n');

    %% Optional Output
    if options.calc_vert_profiles
        fidgp       = fopen([Output_dir,'gap.csv'],'w');              % gap
        fidtc       = fopen([Output_dir,'leaftemp.csv'],'w');         % leaftemp
        fidhl       = fopen([Output_dir,'layer_H.csv'],'w');          % vertical profile
        fidlel      = fopen([Output_dir,'layer_lE.csv'],'w');         % latent heat
        fidal       = fopen([Output_dir,'layer_A.csv'],'w');          %
        fidpl       = fopen([Output_dir,'layer_aPAR.csv'],'w');       %
        fidplC      = fopen([Output_dir,'layer_aPAR_Cab.csv'],'w');       %
        fidrn       = fopen([Output_dir,'layer_Rn.csv'],'w');         %
        if options.calc_fluor
            fidfll = fopen([Output_dir,'layer_fluorescence.csv'],'w');   
            fidfllem = fopen([Output_dir,'layer_fluorescenceEm.csv'],'w'); 
            fidNPQ = fopen([Output_dir,'layer_NPQ.csv'],'w');
        end
%     else
%         delete([Output_dir,'../output/leaftemp.dat'])
%         delete([Output_dir,'../output/layer_H.dat'])
%         delete([Output_dir,'../output/layer_lE.dat'])
%         delete([Output_dir,'../output/layer_A.dat'])
%         delete([Output_dir,'../output/layer_aPAR.dat'])
%         delete([Output_dir,'../output/layer_Rn.dat'])
    end

    if options.calc_fluor
        fidfl       = fopen([Output_dir,'fluorescence.csv'],'w');     % Fluorescence
        if options.calc_PSI
            fidfl1       = fopen([Output_dir,'fluorescencePSI.csv'],'w');     % Fluorescence
            fidfl2       = fopen([Output_dir,'fluorescencePSII.csv'],'w');     % Fluorescence
        end
        fidflh      = fopen([Output_dir,'fluorescence_hemis.csv'],'w');     % Fluorescence
        fidfle         = fopen([Output_dir,'fluorescence_emitted_by_all_leaves.csv'],'w');
        fidfrc         = fopen([Output_dir,'fluorescence_emitted_by_all_photosystems.csv'],'w');
        fidflsu      = fopen([Output_dir,'fluorescence_sunlit.csv'],'w');     % Fluorescence
        fidflsh      = fopen([Output_dir,'fluorescence_shaded.csv'],'w');     % Fluorescence
        fidflsc      = fopen([Output_dir,'fluorescence_scattered.csv'],'w');     % Fluorescence
%     else
%         delete([Output_dir,'fluorescence.dat'])
    end

    if options.calc_directional
        delete([Output_dir,'BRDF/*.dat'])
    end

    if options.calc_planck && options.calc_ebal
        fidplancko      = fopen([Output_dir,'spectrum_obsdir_thermal.csv'],'w');  % spectrum observation direction
        fidplanckh      = fopen([Output_dir,'spectrum_hemis_thermal.csv'],'w');  % spectrum hemispherically integrated
    end

    %% write headers
    if options.save_headers
        fprintf(fidf,'timestep,counter,year,t,Rntot,lEtot,Htot,Rnctot,lEctot,Hctot,Actot,Rnstot,lEstot,Hstot,Gtot,Resp,aPAR,aPAR_Cab,faPAR,aPAR_energyunits,iPAR');
        if options.calc_fluor
            fprintf(fidf,',fluortot fluor_yield');
        end
        fprintf(fidf,'\n');
        
        fprintf(fidf,'#-,-,-,JulianDay,W m-2,W m-2,W m-2,W m-2,W m-2,W m-2,umol m-2 s-1,W m-2,W m-2,W m-2,W m-2,umol m-2 s-1,umol m-2 s-1,umol umol-1,-,W m-2,umol m-2 s-1');
        if options.calc_fluor
            fprintf(fidf,',W m-2 WW^{-1}');
        end
        fprintf(fidf,'\n');

        fprintf(fidt,'timestep,year,t,Ta,Tss(1),Tss(2),Tcave,Tsave\n');
        fprintf(fidt,'#-,-,JulianDay,^oC,^oC,^oC,^oC,^oC\n');

        fprintf(fidra, 'raa,rawc,raws,ustar\n');
        fprintf(fidra, '#sm-1,sm-1,sm-1,ms-1\n');

        fprintf(fidr, 'timestep,year,t,ShortIn,LongIn,HemisOutShort,HemisOutLong,HemisOutTot,Net\n');
        fprintf(fidr, '#-,-,JulianDay,Wm-2,Wm-2,Wm-2,Wm-2,Wm-2,Wm-2\n');

        fprintf(fidfho, '# hemispherically integrated radiation spectrum \n');
        fprintf(fidfho, '# W m-2 um-1 \n');

        fprintf(fidfoo, '# radiance spectrum in observation direction \n');
        fprintf(fidfoo, '# W m-2 sr-1 um-1 \n');

        if options.calc_ebal
            fprintf(fidto, '# thermal BlackBody emission spectrum in observation direction \n');
            fprintf(fidto, '# W m-2 sr-1 um-1 \n');
            if options.calc_planck
                fprintf(fidplancko, '# thermal emission spectrum in observation direction \n');
                fprintf(fidplancko, '# W m-2 sr-1 um-1 \n');

                fprintf(fidplanckh, '# thermal emission spectrum in hemispherical direction \n');
                fprintf(fidplanckh, '# W m-2 sr-1 um-1 \n');
            end
        end

        fprintf(fidwl, '# wavelengths of the spectral output files \n');
        fprintf(fidwl, '# um \n');

        fprintf(fidsi, '# irradiance \n');
        fprintf(fidsi, '# W m-2 um-1\n');

        fprintf(fidref, '# reflectance \n');
        fprintf(fidref, '# fraction of radiation in observation direction *pi / irradiance \n');

        fprintf(fidp, '# Bottom of canopy irradiance in the shaded fraction, and average BOC irradiance \n');
        fprintf(fidp, '# First 2162 columns: shaded fraction. Last 2162 columns: average BOC irradiance. Unit: Wm-2 um-1 \n');

       % fprintf(fidref2, 'reflectance including dynamic Xanthophyll effects \n');
       % fprintf(fidref2, 'fraction of radiation in observation direction *pi / irradiance \n');

        if options.calc_fluor
            fprintf(fidfl, '# fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidfl, '# W m-2 um-1 sr-1\n');
            if options.calc_PSI
                fprintf(fidfl1, '# fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution, for PSI only \n');
                fprintf(fidfl1, '# W m-2 um-1 sr-1\n');
                fprintf(fidfl2, '# fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution, for PSII only \n');
                fprintf(fidfl2, '# W m-2 um-1 sr-1\n');
            end
            fprintf(fidflh, '# hemispherically integrated fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidflh, '# W m-2 um-1 \n');
            fprintf(fidfle, '# total emitted fluorescence by all leaves for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidfle, '# W m-2 um-1 \n');
            fprintf(fidfrc, '# total emitted fluorescence by all photosystems for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidfrc, '# W m-2 um-1 \n');
            fprintf(fidflsu, '# TOC fluorescence contribution from sunlit leaves for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidflsu, '# W m-2 um-1 sr^{-1} \n');
            fprintf(fidflsh, '# TOC fluorescence contribution from shaded leaves for wavelengths of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidflsh, '# W m-2 um-1 sr^{-1} \n');
            fprintf(fidflsc, '# TOC fluorescence contribution from from leaves and soil after scattering for wavelenghts of 640 to 850 nm, with 1 nm resolution \n');
            fprintf(fidflsc, '# W m-2 um-1 sr^{-1} \n');

        end
        if options.calc_vert_profiles
            fprintf(fidgp, '# Fraction leaves in the sun, fraction of observed, fraction of observed&visible per layer \n');
            fprintf(fidgp, '# rows: simulations or time steps, columns: layer numbers  \n');
            fprintf(fidtc, '# leaf temperature of sunlit leaves, shaded leaves, and weighted average leaf temperature per layer \n');
            fprintf(fidtc, '# ^oC ^oC ^oC \n');
            fprintf(fidhl, '# sensible heat flux per layer \n');
            fprintf(fidhl, '# Wm-2\n');
            fprintf(fidlel, '# latent heat flux per layer \n');
            fprintf(fidlel, '# Wm-2\n');
            fprintf(fidal, '# photosynthesis per layer\n');
            fprintf(fidal, '# umol-2s-1\n');
            fprintf(fidpl, '# aPAR per leaf layer \n');
            fprintf(fidpl, '# umol-2s-1 \n');
            fprintf(fidplC, '# aPAR by Cab per leaf layer \n');
            fprintf(fidplC, '# umol-2s-1 \n');
            fprintf(fidrn,'# net radiation per leaf layer \n');
            fprintf(fidrn,'# Wm-2\n');
            if options.calc_fluor
                fprintf(fidfll, '# upward fluorescence per layer\n');
                fprintf(fidfll, '# W m^{-2}\n');

                fprintf(fidNPQ, '# average NPQ = 1-(fm-fo)/(fm0-fo0), per layer \n');
                fprintf(fidNPQ, '# \n');
            end
        end
    end
    %%
%     fprintf(fidwl,'%9.5f ',spectral.wlS);
    dlmwrite([Output_dir,'wl.csv'], spectral.wlS, '-append', 'precision', '%.5f')
    fclose all;
%     dlmwrite([Output_dir,'wl.csv'], spectral.wlS, '-append', 'precision', '%.5f')
end