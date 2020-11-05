function use_rmse(lut_in_path, lut_out_path, full_set_path, var_name)

    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        full_set_path = fullfile(in_dir, 'full_set.csv');
        var_name = 'Actot';
    end
    csv_out = fullfile(fileparts(full_set_path), 'results_rmse.csv');
    map_out = fullfile(fileparts(full_set_path), 'results_rmse.tif');
    fig_out = fullfile(fileparts(full_set_path), 'results_rmse.png');
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    assert(exist(full_set_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you flattened the images with lut.image2csv()?'], full_set_path)

    lut_in = readtable(lut_in_path);
    opt = detectImportOptions(lut_out_path);
    flu = readtable(lut_out_path, opt);
    lut_out = flu.(var_name);
    val_in = readtable(full_set_path);
    
    [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    
    tab = array2table([res, res_std], 'VariableNames', {'Actot', 'Actot_sd'});
    writetable(tab, csv_out)
    fprintf('saved `%s`\n', csv_out)
    
    im = lut.csv2image_plane(val_in, res);
    lut.write_tiff(im, map_out)
    
    lut.plot_image(im, 'RMSE', fig_out)
end