function n_col = output_data_binary(f, k, xyt, rad,  canopy, V, vi, vmax, options, fluxes, meteo, iter,resistance)
%% OUTPUT DATA
% author C. Van der Tol
% date:      30 Nov 2019
% update:    22 Jan 2021 (additional output variables)

%%
if isdatetime(xyt.t)
    get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
    V(46).Val = get_doy(timestamp2datetime(xyt.startDOY));
    V(47).Val = get_doy(timestamp2datetime(xyt.endDOY));
    xyt.t = get_doy(xyt.t);
end

%% Vegetation products
apar_out= [k xyt.year(k) xyt.t(k)  rad.PAR rad.EPAR canopy.LAIsunlit  canopy.LAIshaded...
    canopy.Pntot     canopy.Pnsun     canopy.Pnsha ...
    canopy.Pntot_Cab canopy.Pnsun_Cab canopy.Pnsha_Cab ...
    canopy.Pntot_Car canopy.Pnsun_Car canopy.Pnsha_Car ...
    canopy.Rntot_PAR canopy.Rnsun_PAR canopy.Rnsha_PAR ...
    canopy.Rntot_Cab canopy.Rnsun_Cab canopy.Rnsha_Cab ...
    canopy.Rntot_Car canopy.Rnsun_Car canopy.Rnsha_Car];
n_col.apar = length(apar_out);
fwrite(f.apar_file,apar_out,'double');

veg_out = [k xyt.year(k) xyt.t(k) canopy.A canopy.Ja canopy.ENPQ  canopy.PNPQ canopy.fqe canopy.LST canopy.emis canopy.GPP];
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
    fluor_out = [rad.F685  rad.wl685 rad.F740 rad.wl740 rad.F684 rad.F761 ...
        rad.LoutF rad.EoutF rad.EoutFrc];
    n_col.fluor = length(fluor_out);
    fwrite(f.fluor_file,fluor_out,'double');

    %% Fluorescence spectral outputs
    % fluorescence radiance (L) in observation direction [mW m-2 nm-1 sr-1]
    n_col.fluor_spectrum = length(rad.LoF_);
    fwrite(f.fluor_spectrum_file, rad.LoF_, 'double');

    n_col.sigmaF = length(rad.sigmaF);
    fwrite(f.sigmaF_file, rad.sigmaF, 'double');

    n_col.fRC = length(rad.EoutFrc_);
    fwrite(f.fRC_file, rad.EoutFrc_, 'double');

    n_col.fRCL = length(rad.Femleaves_);
    fwrite(f.fRCL_file, rad.Femleaves_, 'double');

    n_col.fhemis = length(rad.EoutF_);
    fwrite(f.fhemis_file,rad.EoutF_, 'double');

    n_col.Lo2 = length(rad.Lototf_);
    fwrite(f.Lo2_file, rad.Lototf_,'double');

    n_col.rapp = length(rad.reflapp);
    fwrite(f.rapp_file, rad.reflapp,'double');


end

%% reflectance
if isfield(f, 'r_file')
    n_col.r = length(rad.refl);
    fwrite(f.r_file,rad.refl,'double');

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

%% Resistances
resist_out = [resistance.raa, resistance.raws, resistance.rss, resistance.ustar];
n_col.resist = length(resist_out);
fwrite(f.resist_file, resist_out, 'double');

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
