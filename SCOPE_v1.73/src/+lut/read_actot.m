function actot = read_actot(fluxes_path)

    if verLessThan('matlab', '9.1')  % < 2016b
        opt = detectImportOptions(fluxes_path);
        flu = readtable(fluxes_path, opt);
        actot = flu.Actot;
    else  % less error-stable
        c_num_actot = 11;
        flu = dlmread(fluxes_path, ',', 2, 0);
        actot = flu(:, c_num_actot);
    end

end