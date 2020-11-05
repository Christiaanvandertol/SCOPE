function mly_ts = load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_)

%     mly_df = readtable("D:\PyCharm_projects\SCOPE2.0\input\dataset LatinHypercube\input_mSCOPE_timeseries.csv");
    mly_df = readtable(mSCOPE_ts_path);

    t_mly = mly_df.(t_column);
    if all(t_mly <= 367)  % doy is provided
        year_n = unique(year(t_));
        assert(length(year_n) == 1, 'Multiple seasons in mSCOPE are supported only if t_column is TIMESTAMP, not DOY')
        t_mly = datestr(datenum(year_n, 0, t_mly), 'yyyymmddHHMMSS.FFF');
    end
    t_mly = timestamp2datetime(t_mly);

    assert(min(t_mly) <= min(t_) && max(t_mly) >= max(t_), ['Interpolation of mSCOPE values is not possible, '...
        'time range from mSCOPE file does not fully cover range of TS data'])
    
    mly_df.(t_column) = [];
    variables = mly_df.Properties.VariableNames;
    param_names = variables(1:nly:length(variables));

    mly_in_t_ = interp1(t_mly, table2array(mly_df), t_);

    % length(t_) == size(mly_in_t_, 1)
    mly_split = mat2cell(mly_in_t_, length(t_), repmat(nly, 1, length(param_names)));

    mly_ts = cell2struct(mly_split, param_names, 2);
end

