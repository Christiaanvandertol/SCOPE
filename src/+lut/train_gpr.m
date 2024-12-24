function train_gpr(lut_in_path, lut_out_path, var_name)
    if nargin == 0
        in_dir = fullfile('..', 'exercise');
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        var_name = 'Actot';
    end
    mat_out = fullfile(fileparts(lut_in_path), ['gpr_' var_name '.mat']);
    
    assert(exist(lut_in_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you generated the LUT input with lut.generate_lut_input()?'], lut_in_path)
    assert(exist(lut_out_path, 'file') ~= 0, ['Did not find `%s` file.\n'...
        'Have you copied the `fluxes.csv` from `../output` after SCOPE run on lut_in.csv?'], lut_out_path)
    
    lut_in = readtable(lut_in_path);
    
    opt = detectImportOptions(lut_out_path);
    flu = readtable(lut_out_path, opt);
    lut_in.(var_name) = flu.(var_name);
    
    % Cross varidation (train: 70%, test: 30%)
    cv = cvpartition(size(lut_in,1),'HoldOut',0.3);
    idx = cv.test;
    % Separate to training and test data
    lut_train = lut_in(~idx,:);
    lut_test  = lut_in(idx,:);
    
    fprintf('Started training gaussian process regression.\nUsually takes 1 minute.\n')
    gprMdl = fitrgp(lut_train, var_name, 'Standardize', 1);
    % we can overfit but not predict on new data with cvGPR (kfoldGPR)
%     gprMdl = fitrgp(lut_train, var_name, 'Standardize', 1, 'CrossVal', 'on', 'Verbose', 1);
%     res = kfoldPredict(gprMdl);
    save(mat_out, 'gprMdl')
    fprintf('GPR is saved in `%s`\n', mat_out)

    res = predict(gprMdl, lut_test);
    fig_path = fullfile(fileparts(lut_in_path), [var_name '_gpr.png']);
    lut.plot_1to1(lut_test.(var_name), res, 'GPR', fig_path, var_name)
    
end