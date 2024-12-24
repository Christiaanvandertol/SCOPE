function generate_lut_input(tab, n_spectra, outdir)
    % params = generate_lut_input(tab, n_spectra, outdir)
    if nargin == 0
        tab = readtable(fullfile('+lut', 'input_borders.csv'));
        n_spectra = 1000;
        outdir = fullfile('..', 'exercise');
    end
    out_file = fullfile(outdir, 'lut_in.csv');
    out_file_for_ts = fullfile(outdir, 'lut_in_for_ts.csv');
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete/rename it first', out_file)
    
    include = logical(tab.include);
    lb = tab.lower(include)';
    ub = tab.upper(include)';
    varnames = tab.variable(include)';

    % one row - one set of parameters
    lh = lhsdesign(n_spectra, sum(include));

    params = (ub-lb) .* lh + lb;

    if any(strcmp('LIDFa' , varnames))
        % abs(LIDFa + LIDFb) <= 1
        i_lidfa = strcmp('LIDFa', varnames);
        i_lidfb = strcmp('LIDFb', varnames);
        lidfa = params(:, i_lidfa);
        lidfb = params(:, i_lidfb);
        params(:, i_lidfa) = (lidfa + lidfb) / 2;
        params(:, i_lidfb) = (lidfa - lidfb) / 2;
    end

    t = array2table(params);
    t.Properties.VariableNames = varnames;
    writetable(t, out_file)

    t.t = datestr(datetime(2022, 7, 1, 0, 12, 1:n_spectra), 'yyyymmddHHMMSS.FFF');
    if ~any(strcmp('tts' , varnames))
        t.tts = repmat(30, [n_spectra, 1]);  % change here tts
    end
    writetable(t, out_file_for_ts)
    
    if verLessThan('matlab', '9.1')  % < 2016b
        varnames_in = '';
    else
        varnames_in = strjoin(varnames, ', ');
    end
    fprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in)
    fprintf('Saved lut input (parameters) in `%s`\n', out_file)
end