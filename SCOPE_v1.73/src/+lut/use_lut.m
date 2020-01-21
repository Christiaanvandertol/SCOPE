function use_lut(lut_in_path, lut_out_path, full_set_path)

    if nargin == 0
        in_dir = '../exercise';
        lut_in_path = fullfile(in_dir, 'lut_in.csv');
        lut_out_path = fullfile(in_dir, 'lut_out.csv');
        full_set_path = fullfile(in_dir, 'full_set.csv');
    end
    csv_out = fullfile(fileparts(full_set_path), 'results_lut.csv');
    map_out = fullfile(fileparts(full_set_path), 'results_lut.tif');
    fig_out = fullfile(fileparts(full_set_path), 'results_lut.png');

    lut_in = readtable(lut_in_path);
    lut_out = lut.read_actot(lut_out_path);
    val_in = readtable(full_set_path);
    
    [res, res_std] = lut.lut_search(val_in, lut_in, lut_out);
    
    tab = array2table([res, res_std], 'VariableNames', {'Actot', 'Actot_sd'});
    writetable(tab, csv_out)
    fprintf('saved `%s`', csv_out)
    
    im = lut.csv2image_plane(val_in, res);
    imwrite(im, map_out)
    
    figure
    imagesc(im)
    cb = colorbar;
    title(cb, '[\mumol CO_2 cm^{-2} s^{-1}]')
    title('LUT GPP')
    saveas(gcf, fig_out)
end