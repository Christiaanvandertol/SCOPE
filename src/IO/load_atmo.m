function atmo = load_atmo(atmfile, SCOPEspec)
    % default error is not clear enough
    assert(exist(atmfile, 'file') == 2, 'Atmospheric file `%s` does not exist', atmfile)
    [~, ~, ext] = fileparts(atmfile);
    if strcmp(ext, '.atm')
        atmo.M  = aggreg(atmfile, SCOPEspec);
    else
        raddata = load(atmfile);
        atmo.Esun_ = raddata(:,1);
        atmo.Esky_ = raddata(:,2);
    end
end