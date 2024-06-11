function bin_to_csv(fnames, V, vmax, n_col, ns)

%% pars
if sum(vmax>1)
  write_output(['n_pars', {V(vmax>1).Name}], {''}, fnames.pars_file, n_col.pars, ns)
end

%% aPAR
apar_names = {'simulation_number', 'year', 'DoY', 'iPAR', 'iPARE', 'LAIsunlit', 'LAIshaded'...
    'aPARtot', 'aPARsun', 'aPARsha',...
    'aPARCabtot', 'aPARCabsun', 'aPARCabsha',...
    'aPARCartot', 'aPARCarsun', 'aPARCarsha',...
    'aPARtotE', 'aPARsunE', 'aPARshaE',...
    'aPARCabtotE', 'aPARCabsunE', 'aPARCabshaE',...
    'aPARCartotE', 'aPARCarsunE', 'aPARCarshaE',};

apar_units = {'', '', '', 'umol m-2 s-1','W m-2', 'm2m-2','m2m-2',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'umol m-2 s-1', 'umol m-2 s-1', 'umol m-2 s-1',...
    'W m-2','W m-2','W m-2',...
    'W m-2','W m-2','W m-2',...
    'W m-2','W m-2','W m-2'};

write_output(apar_names, apar_units, fnames.apar_file, n_col.apar, ns)

%% veg
veg_names = {'simulation_number', 'year', 'DoY', 'Photosynthesis', 'Electron_transport', 'NPQ_energy',  'NPQ_photon', 'canopy_level_FQE','LST','emis', 'GPP'};
veg_units = {'', '', '', 'umol CO2 m-2 s-1', 'umol m-2 s-1', 'W m-2', 'umol m-2 s-1', 'umol photons (umol photons)-1', 'K','', 'umol CO2 m-2 s-1'};
write_output(veg_names, veg_units, fnames.veg_file, n_col.veg, ns)

%% flu
flu_names = {'simulation_number','nu_iterations', 'year','DoY',...
    'Rnctot','lEctot','Hctot','Actot','Tcave', ...
    'Rnstot','lEstot','Hstot','Gtot','Tsave',...
    'Rntot','lEtot','Htot'};
flu_units = {'', '', '', '',  ...
    'W m-2','W m-2','W m-2','umol m-2 s-1','C',...
    'W m-2','W m-2','W m-2','W m-2','C',...
    'W m-2',' W m-2','W m-2'};
write_output(flu_names, flu_units, fnames.flu_file, n_col.flu, ns)

%% rad
rad_names = {'simulation_number','year','DoY','ShortIn','LongIn','HemisOutShort','HemisOutLong','Lo','Lot','Lote'};
rad_units = {'','','','W m-2','W m-2','W m-2','W m-2','W m-2 sr-1','W m-2 sr-1','W m-2 sr-1'};
write_output(rad_names, rad_units, fnames.rad_file, n_col.rad, ns)

%% fluor
if isfield(fnames, 'fluor_file')
    fluor_names = {'F_1stpeak', 'wl_1stpeak', 'F_2ndpeak', 'wl_2ndpeak', 'F684', 'F761', 'LFtot', 'EFtot', 'EFtot_RC'};
    fluor_units = {'W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','nm','W m-2 um-1 sr-1','W m-2 um-1 sr-1','W m-2 sr-1','W m-2','W m-2'};
    write_output(fluor_names, fluor_units, fnames.fluor_file, n_col.fluor, ns)

    write_output({'fluorescence_spectrum 640:1:850 nm'}, {'W m-2 um-1 sr-1'}, ...
        fnames.fluor_spectrum_file, n_col.fluor_spectrum, ns, true)

    write_output({'escape probability 640:1:850 nm'}, {''}, ...
        fnames.sigmaF_file, n_col.sigmaF, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm hemispherically integrated'}, {'W m-2 um-1'}, ...
        fnames.fhemis_file, n_col.fhemis, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm reabsorption corrected'}, {'W m-2 um-1'}, ...
        fnames.fRC_file, n_col.fhemis, ns, true)

    write_output({'fluorescence_spectrum 640:1:850 nm emission by all leaves'}, {'W m-2 um-1'}, ...
        fnames.fRCL_file, n_col.fhemis, ns, true)

    write_output({'upwelling radiance including fluorescence'}, {'W m-2 um-1 sr-1'}, ...
        fnames.Lo2_file, n_col.Lo2, ns, true)

    write_output({'apparent reflectance'}, {''}, ...
        fnames.rapp_file, n_col.rapp, ns, true)

end

%% reflectance
if isfield(fnames, 'r_file')
    write_output({'reflectance'}, {'pi*upwelling radiance/irradiance'}, ...
        fnames.r_file, n_col.r, ns, true)

    write_output({'rsd'}, {'directional-hemispherical reflectance factor'}, ...
        fnames.rsd_file, n_col.rsd, ns, true)

    write_output({'rdd'}, {'bi-hemispherical reflectance factor'}, ...
        fnames.rdd_file, n_col.rdd, ns, true)

    write_output({'rso'}, {'bi-directional reflectance factor'}, ...
        fnames.rso_file, n_col.rso, ns, true)

    write_output({'rdo'}, {'hemispherical-directional reflectance factor'}, ...
        fnames.rdo_file, n_col.rdo, ns, true)

    %% radiance
    write_output({'hemispherically integrated upwelling radiance'}, {'W m-2 um-1'}, ...
        fnames.Eout_file, n_col.Eout, ns, true)

    write_output({'upwelling radiance excluding fluorescence'}, {'W m-2 um-1 sr-1'}, ...
        fnames.Lo_file, n_col.Lo, ns, true)

    write_output({'direct solar irradiance'}, {'W m-2 um-1'}, ...
        fnames.Esun_file, n_col.Esun, ns, true)

    write_output({'diffuse solar irradiance'}, {'W m-2 um-1'}, ...
        fnames.Esky_file, n_col.Esky, ns, true)
end

%% resistances
write_output({'aerodynamicresistance(ra)', 'raforsoil', 'rss','ustar'}, ...
    {'s m-1','s m-1','s m-1','m s-1'}, ...
    fnames.resist_file, n_col.resist, ns, true)

fclose('all');

%% deleting .bin
structfun(@delete, fnames)
end

function write_output(header, units, bin_path, f_n_col, ns, not_header)
    if nargin == 5
        not_header = false;
    end
    n_csv = strrep(bin_path, '.bin', '.csv');

    f_csv = fopen(n_csv, 'w');
    header_str = [strjoin(header, ','), '\n'];
    if not_header
        header_str = ['#' header_str];
    else
        % it is a header => each column must have one
        assert(length(header) == f_n_col, 'Less headers than lines `%s` or n_col is wrong', bin_path)
    end
    fprintf(f_csv, header_str);
    fprintf(f_csv, ['#' strjoin(units, ','), '\n']);

    f_bin = fopen(bin_path, 'r');
    out = fread(f_bin, 'double');
%     fclose(f_bin);  % + some useconds to execution
    out_2d = reshape(out, f_n_col, ns)';
%     dlmwrite(n_csv, out_2d, '-append', 'precision', '%d'); % SLOW!
    for k=1:ns
        fprintf(f_csv, '%d,', out_2d(k, 1:end-1));
        fprintf(f_csv, '%d\n', out_2d(k, end));  % saves from extra comma
    end
%     fclose(f_csv);
end
