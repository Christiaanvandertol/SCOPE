function [V, xyt, mly_ts, atmo_paths]  = load_timeseries(V, F, xyt, path_input)
    
    %% filenames
    Dataset_dir = ['dataset ' F(5).FileName];
    meteo_ec_csv = F(6).FileName;
    vegetation_retrieved_csv  = F(7).FileName;
    mSCOPE_csv = F(10).FileName;

    t_column = F(strcmp({F.FileID}, 't')).FileName;
    year_column = F(strcmp({F.FileID}, 'year')).FileName;
    atmo_column = F(strcmp({F.FileID}, 'atmos_names')).FileName;
    
    %% read berkeley format dataset
    df = readtable(fullfile(path_input, Dataset_dir, meteo_ec_csv), ...
        'TreatAsEmpty', {'.','NA','N/A'});
%     df = standardizeMissing(df, -9999); > 2013a
    t_ = df.(t_column);

    if iscell(t_)
        t_ = cellfun(@str2num, t_);
    end
    t_(t_ == -9999) = nan;
    
    if all(t_ <= 367)  % doy is provided
        %assert(~isempty(year_column), 'Please, provide year in your .csv')
        if(isempty(year_column)) 
            warning('t is DOY, converting to date with year = 2020, as `year` in .csv was empty')
            year_n = 2020; 
        else
            year_n = df.(year_column);
        end
        % then we calculate ts for you
        t_ = datestr(datenum(year_n, 0, t_), 'yyyymmddHHMMSS.FFF');
    end
    
    t_ = timestamp2datetime(t_);

    xyt.startDOY = timestamp2datetime(xyt.startDOY);
    xyt.endDOY = timestamp2datetime(xyt.endDOY);
    year_n = year(t_);
    
    %% filtering

    assert(~any(isnat(t_)), 'some t_ in timeseries contain NaN, please, remove or provide. You must know the date')
    
    % this filtering also removes NaNs in t_
    % may be not desirable if df_ts is going to be merged directly => assert
    time_i = (t_ >= xyt.startDOY) & (t_ <= xyt.endDOY);
    df_sub = df(time_i, :);

    %% time 
    t_ = t_(time_i);
    xyt.t = t_;
    xyt.year = year_n(time_i);  % for legacy and doy to date convertion

    %% optional interpolation_csv file
    interpolatable_cols = {};
    if ~isempty(vegetation_retrieved_csv)
        df_int = readtable(fullfile(path_input, Dataset_dir, vegetation_retrieved_csv), ...
            'TreatAsEmpty', {'.','NA','N/A'});
        t_int = df_int.(t_column);
        t_int = timestamp2datetime(t_int, year_n);
        if (min(t_int) > max(t_)) || (max(t_int) < min(t_))
            error('Correct the timestamp, please. The timestamp of vegetation_retrieved_csv is outside of the timestamp of meteo_ec_csv, all simulations would be 0 or NaNs')
        elseif (min(t_int) > min(t_)) && (max(t_int) < max(t_))
            warning('the timestamp of vegetation_retrieved_csv starts later and ends earlier than the timestamp of meteo_ec_csv, head and tail simulations will be 0 or NaNs')
        elseif min(t_int) > min(t_)
            warning('the timestamp of vegetation_retrieved_csv starts later than the timestamp of meteo_ec_csv, head simulations will be 0 or NaNs')
        elseif max(t_int) < max(t_)
            warning('the timestamp of vegetation_retrieved_csv ends earlier than the timestamp of meteo_ec_csv, tail simulations will be 0 or NaNs')
        end
        interpolatable_cols = df_int.Properties.VariableNames;
    end
    
    %% optional mSCOPE file
    mly_ts = struct();
    if ~isempty(mSCOPE_csv)
        mSCOPE_ts_path = fullfile(path_input, Dataset_dir, mSCOPE_csv);
        nly = str2double(F(11).FileName);
        mly_ts = load_mSCOPE_ts(mSCOPE_ts_path, nly, t_column, t_);
        mly_ts.nly = nly;
    end
    
    %% make correspondence: F.FileID : index in V struct
    i_empty = cellfun(@isempty, {F.FileName});
    f_ids = {F(~i_empty).FileID};
    f_names = {F(~i_empty).FileName};
    v_names = {V.Name};
    [~, iF, iV] = intersect(f_ids, v_names, 'stable');
    
    %% read fields that were provided (f_ids)
    for i = 1:length(iF)  % == length(iV)
        fi_i = iF(i);
        vi_i = iV(i);
        col_name = char(f_names(fi_i));
        % TODO replace missing by default?
        if any(strcmp(interpolatable_cols, col_name))
            V(vi_i).Val = interp1(t_int, df_int.(col_name), t_);
        else
            tmp = df_sub.(col_name);
            tmp(tmp == -9999) = nan;
            if all(isnan(tmp))
                warning('%s has NaNs along all timestamp. Calculations may fail', col_name)
            end
            V(vi_i).Val = tmp;
        end
    end

    %% special cases
    %% Irradiance files
    atmo_paths = [];
    if ~isempty(atmo_column)
        atmo_paths = fullfile(path_input, Dataset_dir, df_sub.(atmo_column));
    end
    %% tts calculation
    if ~any(strcmp(f_ids, 'tts'))  % tts wasn't read
        vi_tts = strcmp(v_names, 'tts');
        if isdatetime(t_)
            get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
            t_ = get_doy(t_);
        end
        DOY_  = floor(t_);
        time_ = 24*(t_-DOY_);
        if all(time_ == 0)
            time_ = 12;
			      warning('taking midday time to calculate tts')
        end
        ttsR  = calczenithangle(DOY_,time_ - xyt.timezn ,0,0,xyt.LON,xyt.LAT);     %sun zenith angle in rad
        V(vi_tts).Val = min(85, ttsR / pi * 180);     
    end

    %% ea calculation
    if ~any(strcmp(f_ids, 'ea')) && any(strcmp(f_ids, 'Ta'))  % ea wasn't read but Ta was
        ta = V(strcmp(v_names, 'Ta')).Val;
        es = satvap(ta);
        vi_ea = strcmp(v_names, 'ea');
        rh_column = F(strcmp({F.FileID}, 'RH')).FileName;
        vpd_column = F(strcmp({F.FileID}, 'VPD')).FileName;
        if ~isempty(rh_column)
            rh = df_sub.(rh_column);
            rh(rh == -9999) = nan;
            if any(rh > 10)
                rh = rh / 100;    % rh from [0 100] to [0 1]
                warning('converted relative hudimity from [0 100] to [0 1]')
            end
            ea = es .* rh;
            warning('calculated ea from Ta and RH')
            V(vi_ea).Val = ea;
        elseif ~isempty(vpd_column)
            vpd = df_sub.(vpd_column);
            vpd(vpd == -9999) = nan;
            ea = es - vpd;
            warning('calculated ea from Ta and VPD')
            if any(ea < 0)
                warning('some ea < 0, is your VPD in hPa?')
            end
            V(vi_ea).Val = ea;
        end 
    end

    %% units convertion
    %% p
    if any(strcmp(f_ids, 'p'))
        vi_p = strcmp(v_names, 'p');
        p = V(vi_p).Val;
        if any(p < 500)
            p = p * 10;
            warning('converted air pressure from kPa to hPa')
        end
        V(vi_p).Val = p;
    end

    %% smc
    if any(strcmp(f_ids, 'SMC'))
        vi_smc = strcmp(v_names, 'SMC');
        smc = V(vi_smc).Val;
        if any(smc > 1)
            smc = smc / 100;  % SMC from [0 100] to [0 1]
            warning('converted soil moisture content from from [0 100] to [0 1]')
        end     
        V(vi_smc).Val = smc;
    end
end
