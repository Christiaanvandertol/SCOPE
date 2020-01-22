function train_gpr(lut_in_path, lut_out_path)
    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
    end
    mat_out = fullfile(fileparts(lut_in_path), 'gpr.mat');
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    
    lut_in = readtable(lut_in_path);
    lut_in.Actot = lut.read_actot(lut_out_path);
    
    fprintf('Started training gaussian process regression.\nUsually takes 1 minute.\n')
    gprMdl = fitrgp(lut_in, 'Actot', 'Standardize', 1);
    save(mat_out, 'gprMdl')
    fprintf('GPR is saved in `%s`\n', mat_out)
end