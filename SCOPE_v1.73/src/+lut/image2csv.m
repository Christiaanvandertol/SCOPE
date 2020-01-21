function image2csv(im_path, outdir)

%     ret_path = 'D:\PyCharm_projects\demostrator\output\george\lely\2019-12-10-131939\2019-12-10-131939.nc';
%     ret_path = 'D:\PyCharm_projects\demostrator\output\george\vre\2019-12-14-024659';
    
    if nargin == 0
        im_path = '../exercise/images';
        outdir = '../exercise';
        name = 'full_set';
    else
        [~, name, ~] = fileparts(im_path);
    end
    csv_out = fullfile(outdir, sprintf('%s.csv', name));
    var_names = {'Cab', 'LAI'};  % or dir(images) if only .tif
    fprintf('reading images from `%s`\n', im_path)
    
    [~, ~, ext] = fileparts(im_path);
    if strcmp(ext, '.nc')
        get_val = @(x) ncread(im_path, x);
    else
        get_val = @(x) imread(fullfile(im_path, sprintf('%s.tif', x)));
    end

    [r, c] = size(get_val(var_names{1}));
    
    %% subset
%     ir = get80_central(r);
%     ic = get80_central(c);
%     r = length(ir);
%     c = length(ic);
%     
%     figure
%     imagesc(get_val(var_names{1}))
% %     colorbar
%     title(colorbar,'[\mug cm^{-2}]');
%     title(var_names{1})
%     rectangle('Position',[min(ic)-1, min(ir)-1, c, r], 'LineWidth',2,'EdgeColor','r')
    
    %%
    fprintf('Flattening %d rows, %d columns\n', r, c)
    vals = nan(r*c, length(var_names));
    for i=1:length(var_names)
        var = var_names{i};
        v = get_val(var);
%         v = v(ir, ic);
        vals(:, i) = v(:);  % flattening
    end
    
    i_nans = any(isnan(vals), 2);
    df = array2table(vals, 'VariableNames', var_names);
    df_clean = df(~i_nans, :);
    df_clean.(['ind_' num2str(r), '_', num2str(c)]) = find(~i_nans);  % or fprintf()
    
    fprintf('found %d not-nan pixels\n', sum(~i_nans))
    
    %% saving
%     out = struct();
%     out.i_nans = i_nans;
%     out.df_clean = df_clean;
%     out.r = r;
%     out.c = c;
% %     reshape(i_nans, r, c);
%     save(fullfile(outdir, sprintf('%s.mat', name)), 'out')
    
    writetable(df_clean, csv_out)
    fprintf('saved to `%s`\n', csv_out)
end

function ind = get80_central(x, step)
    
    if nargin < 2
        step = 10;
    end

    mid = round(x / 2);
    lo = max(0, mid - step);
    hi = min(x, mid + step);
    ind = lo:hi;

end
