function dt = timestamp2datetime(ts, year_n)
    % ts : char or char array like ['20190101'; '20190102']
    %      int or int array [20190101; 20190102]
    % dt : datetime or datetime array

    if iscell(ts)
        ts = cellfun(@str2num, ts);
    end
    ts(ts == -9999) = nan;
    
    if isnumeric(ts) && all(ts <= 367)  % doy is provided
        warning('t is DOY, converting to date with year = %d, as `year` in .csv was empty', year_n(1))
        ts = datestr(datenum(year_n, 0, ts), 'yyyymmddHHMMSS.FFF');
    end
    
    if isnumeric(ts)
        ts = num2str(ts);
        % char or char array like ['20190101'; '20190102']
    end
    
    switch size(ts, 2)
        case 18  % SS.SSS
%             error('milliseconds timestamp is not implemented due to ambiguity of reading: incomplete ts silently become NaNs')
            format = 'yyyyMMddHHmmss.SSS';
        case 14  % SS
            format = 'yyyyMMddHHmmss';
        case 12  % HH, HR
            format = 'yyyyMMddHHmm';
        case 10  % HH (not standard)
            format = 'yyyyMMddHH';
        case 8  % DD
            % in this case HH == 0, MM == 0, SS == 0
            format = 'yyyyMMdd';
        case 6 % MM, WW
            error('weekly and monthly timestamps are not implemented')
        case 4 % YY
            error('yearly timestamp is not implemented')
        otherwise
            error('format of timestamp is not any truncation of `yyyyMMddHHmm`')
    end
    dt = datetime(ts, 'InputFormat', format);
end