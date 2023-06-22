function [canopy, fluxes, rad, resistance, iter] = fill_output_with_nans(canopy, spectral)

    %% canopy
    
    canopy_names = {"LAIsunlit", "LAIshaded", ...
        "Pntot", "Pnsun", "Pnsha", ...
        "Pntot_Cab", "Pnsun_Cab", "Pnsha_Cab", ...
        "Pntot_Car", "Pnsun_Car", "Pnsha_Car", ...
        "Rntot_PAR", "Rnsun_PAR", "Rnsha_PAR", ...
        "Rntot_Cab", "Rnsun_Cab", "Rnsha_Cab", ...
        "Rntot_Car", "Rnsun_Car", "Rnsha_Car", ...
        "A", "Ja", "ENPQ", "PNPQ", "fqe", "LST", "emis"};
    
    for i=1:length(canopy_names)
        canopy.(canopy_names{i}) = nan;
    end
    
    %% fluxes
    
    fluxes_names = {'Rnctot','lEctot','Hctot','Actot','Tcave', ...
        'Rnstot','lEstot','Hstot','Gtot','Tsave',...
        'Rntot','lEtot','Htot'};
    
    fluxes = struct();
    for i=1:length(fluxes_names)
        fluxes.(fluxes_names{i}) = nan;
    end
    
    %% rad scalars
    
    rad_names = {"PAR", "EPAR",...
        "Eouto","Eoutt","Eoutte", "Lo","Lot","Lote",...
        "F685","wl685","F740","wl740","F684","F761", "LoutF","EoutF","EoutFrc"};
    
    rad = struct();
    for i=1:length(rad_names)
        rad.(rad_names{i}) = nan;
    end
    
    %% rad spectral
    
    rad_names_spectral = {"reflapp", "refl", "rsd", "rdd", "rso", "rdo", ...
        "Eout_", "Lotot_", "Lototf_", "Esun_", "Esky_"};
    
    for i=1:length(rad_names_spectral)
        rad.(rad_names_spectral{i}) = nan(size(spectral.wlS));
    end
    
    %% sif spectral
    
    sif_names = {"LoF_", "sigmaF", "EoutFrc_", "Femleaves_", "EoutF_"};
    
    for i=1:length(sif_names)
        rad.(sif_names{i}) = nan(size(spectral.wlF));
    end
    
    %% resistances
    
    resistances_names = {'raa','raws','rss','ustar'};
    
    resistance = struct();
    for i=1:length(resistances_names)
        resistance.(resistances_names{i}) = nan;
    end

    %% iter
    % in case NaN is in the first row

    iter.counter = nan;

end
