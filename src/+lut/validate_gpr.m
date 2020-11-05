function validate_gpr(gpr_path, val_in_path, val_out_path)

    if nargin == 0
        in_dir = '../exercise';
        gpr_path = fullfile(in_dir, 'gpr.mat');
        val_in_path = fullfile(in_dir, 'validation_in.csv');
        val_out_path = fullfile(in_dir, 'validation_out.csv');
    end
    fig_path = fullfile(fileparts(val_in_path), 'validation_gpr.png');
    
    assert(exist(gpr_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you trained the gaussian process regression (GPR) with lut.train_gpr()?'], gpr_path)
    assert(exist(val_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the validation input with lut.pick100()?'], val_in_path)
    assert(exist(val_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on validation_in.csv?'], val_out_path)

    gpr = load(gpr_path);
    gprMdl = gpr.gprMdl;
    val_in = readtable(val_in_path);
    val_out = lut.read_actot(val_out_path);
    
    tic
    res = predict(gprMdl, val_in);
    toc
    
    lut.plot_1to1(val_out, res, 'GPR', fig_path)

end