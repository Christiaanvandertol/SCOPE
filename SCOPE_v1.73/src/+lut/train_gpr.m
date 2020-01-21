function train_gpr(lut_in_path, lut_out_path)
    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
    end
    mat_out = fullfile(fileparts(lut_in_path), 'gpr.mat');
    
    lut_in = readtable(lut_in_path);
    lut_in.Actot = lut.read_actot(lut_out_path);
    
    gprMdl = fitrgp(lut_in, 'Actot', 'Standardize', 1);
    save(mat_out, 'gprMdl')
    fprintf('GPR is saved in `%s`\n', mat_out)
end