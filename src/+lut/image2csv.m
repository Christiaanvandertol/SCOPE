function image2csv(im_path, out_dir, gpr_name, csv_name)
    
    if nargin == 0
        out_dir = fullfile('..', 'exercise');
        im_path = fullfile(out_dir, 'images');
        gpr_name = 'gpr_Actot.mat';  % to read variable names
        csv_name = 'full_set.csv';
    end
    assert(exist(im_path, 'dir') == 7, '`%s` does not exist. Please, put images (LAI.tif ...) into that folder', im_path)

    gpr_path = fullfile(out_dir, gpr_name);
    gpr = load(gpr_path);
    gpr = gpr.gprMdl;
    var_names = gpr.PredictorNames;
    
    [~, ~, ext] = fileparts(im_path);
    if strcmp(ext, '.nc')
        get_val = @(x) ncread(im_path, x);
    else
        get_val = @(x) imread(fullfile(im_path, sprintf('%s.tif', x)));
    end

    [r, c] = size(get_val(var_names{1}));
    
    %% flattening
    fprintf('Flattening %d rows, %d columns\n', r, c)
    vals = nan(r*c, length(var_names));
    for i=1:length(var_names)
        var = var_names{i};
        fprintf('%s ', var)
        v = get_val(var);
        [r_i, c_i] = size(v);
        if (r_i ~= 1) & (c_i ~= 1)
            assert((r_i == r) & (c_i == c), ...
                'The number of rows (r_i=%i) or columns (c_i=%i) in %s.tif does not match the expected (r=%d, c=%d)', ...
                r_i, c_i, var, r, c)
        end
        vals(:, i) = v(:);  % flattening
    end
    
    i_nans = any(isnan(vals), 2);
    df = array2table(vals, 'VariableNames', var_names);
    df_clean = df(~i_nans, :);
    df_clean.(['ind_' num2str(r), '_', num2str(c)]) = find(~i_nans);  % or fprintf()
    
    fprintf('\nfound %d not-nan pixels\n', sum(~i_nans))
    
    %% saving
%     out = struct();
%     out.i_nans = i_nans;
%     out.df_clean = df_clean;
%     out.r = r;
%     out.c = c;
% %     reshape(i_nans, r, c);
%     save(fullfile(outdir, sprintf('%s.mat', name)), 'out')
    
    csv_out_path = fullfile(out_dir, csv_name);
    writetable(df_clean, csv_out_path)
    fprintf('saved to `%s`\n', csv_out_path)
end
