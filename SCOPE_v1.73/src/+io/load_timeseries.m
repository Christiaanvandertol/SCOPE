function [V, xyt]  = load_timeseries(V, F, xyt, path_input, interpolatable)
    
    %% filenames
    Dataset_dir         = ['dataset ' F(5).FileName];
    berkeley_file       = F(6).FileName;
    interpolation_csv   = F(7).FileName;

    t_column = F(strcmp({F.FileID}, 't')).FileName;
    rh_column = F(strcmp({F.FileID}, 'RH')).FileName;

    %% read berkeley format dataset
    df = readtable(fullfile(path_input, Dataset_dir, berkeley_file), ...
        'TreatAsEmpty', {'.','NA','N/A'});
    t_ = df.(t_column);
    year_n = ones(size(t_)) * 2006;
    if any(t_ > 367)  % is timestamp
        % t = datestr(datenum(2006, 0, df.t), 'yyyymmddHHMMSS.FFF');
        t_ = io.timestamp2datetime(t_);
        xyt.startDOY = io.timestamp2datetime(xyt.startDOY);
        xyt.endDOY = io.timestamp2datetime(xyt.endDOY);
        year_n = year(t_);
    end
    
    %% filtering
    time_i = (t_ >= xyt.startDOY) & (t_ <= xyt.endDOY);   

    % SW < 0 is needed for soil heat model but may cause bad Actot flux
    % is_rin_high = true(size(df, 1), 1);
    % if ~isempty(Rin_file)
    %      rin       = replace_9999(df.(Rin_file));
    %      is_rin_high = rin >= 50;
    % end
    % i = time_i & is_rin_high;
    df_sub = df(time_i, :);
    if isempty(df_sub)
        error('The subset of the dataframe for timeseries is empty. Are startDOY, endDOY correct?')
    end

    %% time 
    t_ = t_(time_i);
    xyt.t = t_;
    xyt.year = year_n(time_i);  % for legacy and doy to date convertion

    %% optional interpolation_csv file
    interpolation_present = false;
    if ~isempty(interpolation_csv)
        interpolation_present = true;
        df_int = readtable(fullfile(path_input, Dataset_dir, interpolation_csv), ...
            'TreatAsEmpty', {'.','NA','N/A'});
        t_int = df_int.(t_column);
        if any(t_int > 367)
            t_int = io.timestamp2datetime(t_int);
        end
        assert(min(t_) >= min(t_int) & max(t_) <= max(t_int), '`interpolation_csv` timestamp is outside `ec_file_berkeley` timestamp')
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
        if interpolation_present && any(strcmp(interpolatable, f_ids(fi_i)))
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
    %% tts calculation
    if ~any(strcmp(f_ids, 'tts'))  % tts wasn't read
        vi_tts = strcmp(v_names, 'tts');
        if isdatetime(t_)
            get_doy = @(x) juliandate(x) - juliandate(datetime(year(x), 1, 0));
            t_ = get_doy(t_);
        end
        DOY_  = floor(t_);
        time_ = 24*(t_-DOY_);
        ttsR  = equations.calczenithangle(DOY_,time_ - xyt.timezn ,0,0,xyt.LON,xyt.LAT);     %sun zenith angle in rad
        V(vi_tts).Val = min(85, ttsR / pi * 180);     
    end

    %% ea calculation
    if ~any(strcmp(f_ids, 'ea')) && ~any(strcmp(f_ids, 'Ta'))  % ea wasn't read but Ta was
        ta = V(strcmp(v_names, 'Ta')).Val;
        vi_ea = strcmp(v_names, 'ea');
        if ~isempty(rh_column)
            rh = replace_9999(df_sub.(rh_file));
            if any(rh > 10)
                rh = rh / 100;    % rh from [0 100] to [0 1]
            end
            V(vi_ea).Val = equations.satvap(ta) .* rh;
            warning('calculated ea from Ta and RH')
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
