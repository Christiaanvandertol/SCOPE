function [res, res_std] = lut_search(params, lut_params, response)

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

    actot = response(ind);

    res = nanmedian(actot, 2);
    res_std = nanstd(actot, 0, 2);  % (v, flag, dim)

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


