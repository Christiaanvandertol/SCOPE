function change_detection(lut_im_path, gpr_im_path)

    if nargin == 0
        in_dir = '../exercise';
        lut_im_path = fullfile(in_dir, 'results_lut.tif');
        gpr_im_path = fullfile(in_dir, 'results_gpr.tif');
    end
    out_path = fullfile(fileparts(lut_im_path), 'results_lut_gpr.png');
    
    lut_im = imread(lut_im_path);
    gpr_im = imread(gpr_im_path);
    
    im_dif = lut_im - gpr_im;
%     im_dif(round(im_dif) == 0) = nan;
    lut.plot_image(im_dif, '(LUT - GPR)', out_path)

end