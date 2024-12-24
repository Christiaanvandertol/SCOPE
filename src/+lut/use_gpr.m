function use_gpr(gpr_path, csv_im_path, im_crs_path)

    if nargin == 0
        in_dir = fullfile('..', 'exercise');
        gpr_path = fullfile(in_dir, 'gpr_Actot.mat');
        csv_im_path = fullfile(in_dir, 'full_set.csv');
        im_crs_path = fullfile(in_dir, 'images', 'LAI.tif');
    end
    
    assert(exist(gpr_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you trained the gaussian process regression (GPR) with lut.train_gpr()?'], gpr_path)
    assert(exist(csv_im_path, 'file') ~= 0, ['Did not find `%s` file.\n'... 
        'Have you flattened the images with lut.image2csv()?'], csv_im_path)

    gpr = load(gpr_path);
    gprMdl = gpr.gprMdl;
    var_out_name = gprMdl.ResponseName;
    val_in = readtable(csv_im_path);

    fprintf('Working, usually takes around 1 minute\n')
    res = predict(gprMdl, val_in);

    csv_out = fullfile(fileparts(csv_im_path), sprintf('results_gpr_%s.csv', var_out_name));
    map_out = fullfile(fileparts(csv_im_path), sprintf('results_gpr_%s.tif', var_out_name));
    fig_out = fullfile(fileparts(csv_im_path), sprintf('results_gpr_%s.png', var_out_name));

    tab = array2table(res, 'VariableNames', {var_out_name}); 
    writetable(tab, csv_out)
    fprintf('saved `%s`\n', csv_out)
    
    im = lut.csv2image_plane(val_in, res);
    lut.write_tiff(im, map_out, im_crs_path)
    
    lut.plot_image(im, 'GPR', fig_out)
end