function validate_rmse(lut_in_path, lut_out_path, val_in_path, val_out_path)

    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        val_in_path = fullfile(in_dir, 'validation_in.csv');
        val_out_path = fullfile(in_dir, 'validation_out.csv');
    end
    fig_path = fullfile(fileparts(val_in_path), 'validation_rmse.png');
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    assert(exist(val_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the validation input with lut.pick100()?'], val_in_path)
    assert(exist(val_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?'], val_out_path)

    lut_in = readtable(lut_in_path);
    lut_out = lut.read_actot(lut_out_path);
    val_in = readtable(val_in_path);
    val_out = lut.read_actot(val_out_path);
    
    [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    
    lut.plot_1to1(val_out, res, 'RMSE', fig_path)

end