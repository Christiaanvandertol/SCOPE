function latin_hypercube_input(tab, n_spectra, outdir)
    if nargin == 0
        outdir = fullfile('input', 'dataset for_verification');
        tab = readtable(fullfile(outdir, 'input_borders.csv'));
        n_spectra = 30;
    end
    out_file = fullfile(outdir, 'lh_ts.csv');
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete it first', out_file)
    
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
    t.t = (1:size(t, 1))';
    writetable(t, out_file)
    
    if verLessThan('matlab', '9.1')  % < 2016b
        varnames_in = '';
    else
        varnames_in = strjoin(varnames, ', ');
    end
    fprintf('Sampled %i parameters: %s\n', length(varnames), varnames_in)
    fprintf('Saved lut input (parameters) in `%s`\n', out_file)
end