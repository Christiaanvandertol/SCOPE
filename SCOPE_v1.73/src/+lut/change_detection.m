function change_detection(lut_im_path, gpr_im_path)

    if nargin == 0
        in_dir = '../exercise';
        lut_im_path = fullfile(in_dir, 'results_rmse.tif');
        gpr_im_path = fullfile(in_dir, 'results_gpr.tif');
    end
    out_path = fullfile(fileparts(lut_im_path), 'results_rmse_gpr.png');
    
    lut_im = imread(lut_im_path);
    gpr_im = imread(gpr_im_path);
    
    im_dif = lut_im - gpr_im;
    lut.plot_image(im_dif, '(RMSE - GPR)', out_path)
end