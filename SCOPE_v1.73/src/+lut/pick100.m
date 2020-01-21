function pick100(full_set_path, n_pixels)
    if nargin == 0
        full_set_path = '../exercise/full_set.csv';
        n_pixels = 100;
    end
    out_file = fullfile(fileparts(full_set_path), 'validation_in.csv');
    
    assert(exist(out_file, 'file') == 0, '`%s` file already exists, delete it first', out_file)
    
    df = readtable(full_set_path);
    ind = randi(size(df, 1), n_pixels, 1);
    df_out = df(ind, :);
    
    writetable(df_out, out_file)
    fprintf('saved 100 pixels to `%s`\n', out_file)
end