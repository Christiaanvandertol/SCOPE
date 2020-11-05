function change_detection(rmse_im_path, gpr_im_path)

    if nargin == 0
        in_dir = '../exercise';
        rmse_im_path = fullfile(in_dir, 'results_rmse.tif');
        gpr_im_path = fullfile(in_dir, 'results_gpr.tif');
    end
    out_path = fullfile(fileparts(rmse_im_path), 'results_rmse_gpr.png');
    
    assert(exist(rmse_im_path, 'file') ~= 0, ['Did not find `%s` image.\n'... 
        'Have you retrieved on full_set.csv with lut.use_rmse()?'], rmse_im_path)
    assert(exist(gpr_im_path, 'file') ~= 0, ['Did not find `%s` image.\n'... 
        'Have you retrieved on full_set.csv with lut.use_gpr()?'], gpr_im_path)
    
    lut_im = imread(rmse_im_path);
    gpr_im = imread(gpr_im_path);
    
    im_dif = lut_im - gpr_im;
    lut.plot_image(im_dif, '(RMSE - GPR)', out_path)
end