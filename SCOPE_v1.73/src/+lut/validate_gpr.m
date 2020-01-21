function validate_gpr(gpr_path, val_in_path, val_out_path)

    if nargin == 0
        in_dir = '../exercise';
        gpr_path = fullfile(in_dir, 'gpr.mat');
        val_in_path = fullfile(in_dir, 'validation_in.csv');
        val_out_path = fullfile(in_dir, 'validation_out.csv');
    end
    fig_path = fullfile(fileparts(val_in_path), 'validation_gpr.png');

    gpr = load(gpr_path);
    gprMdl = gpr.gprMdl;
    val_in = readtable(val_in_path);
    val_out = lut.read_actot(val_out_path);
    
    res = predict(gprMdl, val_in);
    
    lut.plot_1to1(val_out, res, 'GPR', fig_path)

end