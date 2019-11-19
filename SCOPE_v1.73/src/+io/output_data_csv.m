function output_data_csv(Output_dir, options, k, iter, xyt, fluxes, rad, thermal, gap, meteo, spectral, V, vi, vmax, profiles, directional, angles)
%% OUTPUT DATA
% author C. Van der Tol
% modified:      31 Jun 2008: (CvdT) included Pntot in output fluxes.dat
% last modified: 04 Aug 2008: (JT)   included variable output directories
%                31 Jul 2008: (CvdT) added layer_pn.dat
%                19 Sep 2008: (CvdT) spectrum of outgoing radiation
%                19 Sep 2008: (CvdT) Pntot added to fluxes.dat
%                15 Apr 2009: (CvdT) Rn added to vertical profiles
%                03 Oct 2012: (CvdT) included boolean variabel calcebal
%                04 Oct 2012: (CvdT) included reflectance and fPAR
%                10 Mar 2013: (CvdT) major revision: introduced structures
%                22 Nov 2013: (CvdT) added additional outputs
%%
if isdatetime(xyt.t)
    get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
    V(46).Val = get_doy(io.timestamp2datetime(xyt.startDOY));
    V(47).Val = get_doy(io.timestamp2datetime(xyt.endDOY));
    xyt.t = get_doy(xyt.t);
end

write_csv = @(filename, precision, matrix) dlmwrite(filename, matrix, '-append', 'precision', precision); % , 'precision', 9);
% write_csv = @(filename, matrix) writematrix(matrix, filename);  % can't append

%% Standard output

% fluxes
flux_out = [k iter.counter xyt.year(k) xyt.t(k) fluxes.Rntot fluxes.lEtot fluxes.Htot fluxes.Rnctot fluxes.lEctot, ...
    fluxes.Hctot fluxes.Actot fluxes.Rnstot fluxes.lEstot fluxes.Hstot fluxes.Gtot fluxes.Resp, ...
    1E6*fluxes.aPAR  1E6*fluxes.aPAR_Cab fluxes.aPAR/rad.PAR fluxes.aPAR_Wm2 1E6*rad.PAR];
if options.calc_fluor
    flux_out = [flux_out rad.Eoutf,  rad.Eoutf./fluxes.aPAR_Wm2];
end
write_csv([Output_dir,'fluxes.csv'], '%.4f' , flux_out)

% surftemp
write_csv([Output_dir,'surftemp.csv'], '%.4f', ...
    [k xyt.year(k) xyt.t(k) thermal.Ta thermal.Ts(1) thermal.Ts(2) thermal.Tcave thermal.Tsave]);

% aerodyn
write_csv([Output_dir,'aerodyn.csv'], '%.4f', [thermal.raa, thermal.rawc, thermal.raws, thermal.ustar]);

% radiation
write_csv([Output_dir,'radiation.csv'], '%.4f', ...
    [k xyt.year(k) xyt.t(k) meteo.Rin meteo.Rli rad.Eouto ... 
    rad.Eoutt + rad.Eoutte  rad.Eouto+rad.Eoutt + rad.Eoutte fluxes.Rntot]);

% spectrum (added on 19 September 2008)
write_csv([Output_dir,'spectrum_hemis_optical.csv'], '%.5f', rad.Eout_'); 
write_csv([Output_dir,'spectrum_obsdir_optical.csv'], '%.5f', rad.Lo_');

if options.calc_ebal
    write_csv([Output_dir,'spectrum_obsdir_BlackBody.csv'], '%.2f', rad.LotBB_');
    
    if options.calc_planck
        write_csv([Output_dir,'spectrum_hemis_thermal.csv'], '%.2f', rad.Eoutte_');
        write_csv([Output_dir,'spectrum_obsdir_thermal.csv'], '%.2f', rad.Lot_');
    end
end

write_csv([Output_dir,'irradiance_spectra.csv'], '%.4f', meteo.Rin*(rad.fEsuno+rad.fEskyo)');
write_csv([Output_dir,'BOC_irradiance.csv'], '%.0f', [rad.Emin_(61,:),rad.Emin_(61,:)+(rad.Esun_*gap.Ps(61)')']);

reflectance = pi*rad.Lo_./(rad.Esun_+rad.Esky_);
reflectance(spectral.wlS>3000) = NaN;
write_csv([Output_dir,'reflectance.csv'], '%.5f', reflectance');

% input and parameter values (added June 2012)
V_this = nan(size(V));
for i = 1:length(V)
    V_this(i) = V(i).Val(vi(i));
end
write_csv([Output_dir,'pars_and_input.csv'], '%.3f', V_this);


k2 = find(vmax>1);  % really needed for the first one, later vi > 1
V_short = nan(1,length(k2));
for i = 1:length(k2)
    V_short(i) = V(k2(i)).Val(vi(k2(i)));
end
write_csv([Output_dir,'pars_and_input_short.csv'], '%.5f', V_short);

%% Optional Output

if options.calc_vert_profiles
    % gap
    write_csv([Output_dir,'gap.csv'], '%.2f', [gap.Ps gap.Po gap.Pso]);
    write_csv([Output_dir,'layer_aPAR.csv'], '%.2f', [1E6*profiles.Pn1d' 0]);
    write_csv([Output_dir,'layer_aPAR_Cab.csv'], '%.2f', [1E6*profiles.Pn1d_Cab' 0]);
    
    if options.calc_ebal
        % leaftemp
        write_csv([Output_dir,'leaftemp.csv'], '%.2f', [profiles.Tcu1d' profiles.Tch' profiles.Tc1d']);
        write_csv([Output_dir,'layer_h.csv'], '%.2f', [profiles.Hc1d' fluxes.Hstot]);
        write_csv([Output_dir,'layer_le.csv'], '%.2f', [profiles.lEc1d' fluxes.lEstot]);
        write_csv([Output_dir,'layer_a.csv'], '%.2f', [profiles.A1d' fluxes.Resp]);
        write_csv([Output_dir,'layer_NPQ.csv'], '%.2f', [profiles.qE' 0]); 
        write_csv([Output_dir,'layer_rn.csv'], '%.2f', [profiles.Rn1d' fluxes.Rnstot]);
    end
    
    if options.calc_fluor
        write_csv([Output_dir,'layer_fluorescence.csv'], '%.2f', profiles.fluorescence');
    end
end

if options.calc_fluor% && options.calc_ebal
%     for j=1:size(spectral.wlF,1)  % always == 1
    % added here ' everywhere
    write_csv([Output_dir,'fluorescence.csv'], '%.4f', rad.LoF_');
    if options.calc_PSI
        write_csv([Output_dir,'fluorescencePSI.csv'], '%.4f', rad.LoF1_');
        write_csv([Output_dir,'fluorescencePSII.csv'], '%.4f', rad.LoF2_');
    end
    write_csv([Output_dir,'fluorescence_hemis.csv'], '%.4f', rad.Fhem_');
    write_csv([Output_dir,'fluorescence_emitted_by_all_leaves.csv'], '%.4f', rad.Fem_');
    write_csv([Output_dir,'fluorescence_emitted_by_all_photosystems.csv'], '%.4f', rad.Femtot');
    write_csv([Output_dir,'fluorescence_sunlit.csv'], '%.4f', sum(rad.LoF_sunlit,2)');
    write_csv([Output_dir,'fluorescence_shaded.csv'], '%.4f', sum(rad.LoF_shaded,2)');
    write_csv([Output_dir,'fluorescence_scattered.csv'], '%.4f', [sum(rad.LoF_scattered,2)+sum(rad.LoF_soil,2)]');
end

%%
if options.calc_directional && options.calc_ebal
    Output_angle    =   [directional.tto_ov';  directional.psi_ov'; angles.tts*ones(size(directional.psi_ov'))];
    Output_brdf     =   [spectral.wlS'     directional.brdf_];
    if options.calc_planck
        Output_temp =   [spectral.wlT'    directional.Lot_];
    else
        Output_temp =   [directional.BrightnessT];
    end
    if options.calc_fluor
        Output_fluor = [spectral.wlF'     directional.LoF_];
    end
    
    save([Output_dir,'Directional/',sprintf('BRDF (SunAngle %2.2f degrees).dat',angles.tts)],'Output_brdf' ,'-ASCII','-TABS')
    save([Output_dir,'Directional/',sprintf('Angles (SunAngle %2.2f degrees).dat',angles.tts)],'Output_angle','-ASCII','-TABS')
    save([Output_dir,'Directional/',sprintf('Temperatures (SunAngle %2.2f degrees).dat',angles.tts)],'Output_temp','-ASCII','-TABS')
    
    if options.calc_fluor
        save([Output_dir,'Directional/',sprintf('Fluorescence (SunAngle %2.2f degrees).dat',angles.tts)],'Output_fluor','-ASCII','-TABS')
    end
    
    fiddirtir       =   fopen([Output_dir,'Directional/','read me.txt'],'w');
    fprintf(fiddirtir,'The Directional data is written in three files: \r\n');
    fprintf(fiddirtir,'\r\n- Angles: contains the directions. \r\n');
    fprintf(fiddirtir,'   * The 1st row gives the observation zenith  angles\r\n');
    fprintf(fiddirtir,'   * The 2nd row gives the observation azimuth angles\r\n');
    fprintf(fiddirtir,'   * The 3rd row gives the solar       zenith  angles\r\n');
    fprintf(fiddirtir,'\r\n- Temperatures: contains the directional brightness temperature. \r\n');
    fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the brightness temperature values (except for broadband)\r\n');
    fprintf(fiddirtir,'   * The 2nd column gives the Tb values corresponding to the directions given by first column in the Angles file\r\n');
    fprintf(fiddirtir,'\r\n- BRDF: contains the bidirectional distribution functions values. \r\n');
    fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the BRDF values\r\n');
    fprintf(fiddirtir,'   * The 2nd column gives the BRDF values corresponding to the directions given by first column in the Angles file\r\n');
    fclose(fiddirtir);
end

%%
fclose all;
end