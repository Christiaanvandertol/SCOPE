function [res, res_std] = lut_search(params, lut_params, response)

% lut_params = readtable('lut_in.csv');
% params = readtable("D:\PyCharm_projects\SCOPE-CvT\data\input\dataset exercise\2019-12-14-024659.csv");

%% check validity: all lut params in params
p_names = params.Properties.VariableNames;
lut_names = lut_params.Properties.VariableNames;

% absent = setdiff(lut_names, p_names);
% assert(isempty(absent), '%s parameter must be in input table, because it is in lut\n', absent)
% 
% extra = setdiff(p_names, lut_names);
% if ~isempty(extra)
%     warning('%s parameters are not in lut, thus will not effect\n', strjoin(extra, ','))
% end

%% normalization to standardize RMSE, column alignment
p_ordered = nan(size(params, 1), length(lut_names));
for i=1:length(lut_names)
    v = lut_names{i};
    v_min = min(lut_params.(v));
    v_max = max(lut_params.(v));
    lut_params.(v) = (lut_params.(v) - v_min) / (v_max - v_min);
    p_ordered(:, i) = (params.(v) - v_min) / (v_max - v_min);
end
lut_arr = table2array(lut_params);

tic
ind = arrayfun(@(x) top_indices(x, p_ordered, lut_arr), 1:size(p_ordered, 1), 'UniformOutput', false);
% rmses = arrayfun(@(x) top_indices(x, p_ordered, lut_arr), 1:size(p_ordered, 1), 'UniformOutput', false);
ind = cell2mat(ind');
toc

% flu_path = "D:\PyCharm_projects\SCOPE-CvT\SCOPE_v1.73\fluxes.csv";
% opt = detectImportOptions(flu_path);
% flu = readtable(flu_path, opt);
% 
% actot = flu.Actot(ind);  % mean, median, std, replace nan?

actot = response(ind);

res = nanmedian(actot, 2);
res_std = nanstd(actot, 0, 2);  % (v, flag, dim)

% ori_path = "D:\PyCharm_projects\SCOPE-CvT\SCOPE_v1.73\output\exercise_2020-01-15-1724\fluxes.csv";
% opt = detectImportOptions(ori_path);
% ori = readtable(ori_path, opt);
% 
% figure
% errorbar(ori.Actot, res, res_std, 'o')
% title('LUT vs SCOPE')

end


function ind = top_indices(i, meas_arr, lut_arr)

    if mod(i, 100000) == 0
        fprintf('done %i / %i  pixels\n', i, length(meas_arr))
    end
    
    rmses = sqrt(mean((meas_arr(i, :) - lut_arr) .^ 2, 2));
    
    [~, I] = sort(rmses);
    ind = I(1:10)';  % top 10
%     rmses = rmses(ind);
end


