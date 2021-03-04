function n_col = output_data_binary(f, k, xyt, rad,  canopy, V, vi, vmax, options, fluxes, meteo, iter)
%% OUTPUT DATA
% author C. Van der Tol
% date:      30 Nov 2019 

%%
if isdatetime(xyt.t)
    get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
    V(46).Val = get_doy(timestamp2datetime(xyt.startDOY));
    V(47).Val = get_doy(timestamp2datetime(xyt.endDOY));
    xyt.t = get_doy(xyt.t);
end

%% Vegetation products
veg_out = [k xyt.year(k) xyt.t(k) canopy.Pntot canopy.Pntot_Cab canopy.Rntot_Cab canopy.A canopy.Ja canopy.ENPQ  canopy.LST];
n_col.veg = length(veg_out);
fwrite(f.veg_file,veg_out,'double');

%% Fluxes product
flu_out = [k iter.counter xyt.year(k) xyt.t(k) cell2mat(struct2cell(fluxes))'];
n_col.flu = length(flu_out);
fwrite(f.flu_file,flu_out,'double');

%% Radiation
rad_out = [k xyt.year(k) xyt.t(k) meteo.Rin, meteo.Rli, rad.Eouto, rad.Eoutt + rad.Eoutte, ...
    rad.Lo, rad.Lot, rad.Lote];
n_col.rad = length(rad_out);
fwrite(f.rad_file,rad_out,'double');

%% Fluorescence scalar outputs
if options.calc_fluor
    fluor_out = [rad.F685  rad.wl685 rad.F740 rad.wl740 rad.F687 rad.F760 ...
        rad.LoutF rad.EoutF rad.EoutFrc];
    n_col.fluor = length(fluor_out);
    fwrite(f.fluor_file,fluor_out,'double');
    
    %% Fluorescence spectral outputs
    % fluorescence radiance (L) in observation direction [mW m-2 nm-1 sr-1]
    n_col.fluor_spectrum = length(rad.LoF_);
    fwrite(f.fluor_spectrum_file, rad.LoF_, 'double');
    
    n_col.sigmaF = length(rad.sigmaF);
    fwrite(f.sigmaF_file, rad.sigmaF, 'double');
    
    n_col.fhemis = length(rad.EoutF_);
    fwrite(f.fhemis_file,rad.EoutF_, 'double');
    
    n_col.Lo2 = length(rad.Lototf_);
    fwrite(f.Lo2_file, rad.Lototf_,'double');
end

if options.save_spectra 
    %% reflectance
    n_col.r = length(canopy.reflectance);
    fwrite(f.r_file,canopy.reflectance,'double');

    n_col.rsd = length(rad.rsd);
    fwrite(f.rsd_file,rad.rsd,'double');

    n_col.rdd = length(rad.rdd);
    fwrite(f.rdd_file,rad.rdd,'double');

    n_col.rso = length(rad.rso);
    fwrite(f.rso_file,rad.rso,'double');

    n_col.rdo = length(rad.rdo);
    fwrite(f.rdo_file,rad.rdo,'double');

    %% Radiance
    n_col.Eout = length(rad.Eout_);
    fwrite(f.Eout_file,rad.Eout_,'double');

    n_col.Lo = length(rad.Lotot_);
    fwrite(f.Lo_file, rad.Lotot_, 'double');

    n_col.Esun = length(rad.Esun_);
    fwrite(f.Esun_file, rad.Esun_, 'double');

    n_col.Esky = length(rad.Esky_);
    fwrite(f.Esky_file, rad.Esky_, 'double');
end

%% pars
k2 = find(vmax>1);  % really needed for the first one, later vi > 1
V_short = nan(1,length(k2)+1);
V_short(1) = length(k2);
for i = 1:length(k2)
    V_short(i+1) = V(k2(i)).Val(vi(k2(i)));
end
n_col.pars = length(V_short);
fwrite(f.pars_file, V_short,'double');

end
