function im = csv2image_plane(val_in, res)
% to read r, c from column name ind_r_c
    vars = val_in.Properties.VariableNames;
    i_ind = ~cellfun(@isempty, strfind(vars, 'ind'));  % ignore readability, or check > 2016b
    splt =  strsplit(vars{i_ind}, '_');
    r = str2num(splt{2});
    c = str2num(splt{3});
    im = nan(r, c);
    pix_ind = table2array(val_in(:, i_ind));
    im(pix_ind) = res;
end