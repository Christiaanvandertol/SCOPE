path_x = fullfile('output', 'SCOPE_sparse_2020-06-07-1346', 'fluxes.csv');
path_y = fullfile('output', 'sunshade_2020-06-07-2234', 'fluxes.csv');
path_y = fullfile('output', 'bigleaf_2020-06-08-0022', 'fluxes.csv');
path_y = fullfile('output', 'bigleaf_my_2020-06-08-0030', 'fluxes.csv');

opts = detectImportOptions(path_x);
df_x = readtable(path_x, opts);
df_y = readtable(path_y, opts);


flu_names = {'Rnctot','lEctot','Hctot','Actot','Tcave', ...
    'Rnstot','lEstot','Hstot','Gtot','Tsave',...
    'Rntot','lEtot','Htot','rss'};

figure
for i=1:length(flu_names)
    subplot(3, 5, i)
    flux = flu_names{i};
    x = df_x.(flux);
    y = df_y.(flux);
    plot(x, y, 'o', 'MarkerFaceColor', 'r')
%     title(flux)
    hold on
    
    % liner model
    i_nans = isnan(x);
    lm_full = fitlm(x(~i_nans), y(~i_nans));   
    lm = polyfit(x(~i_nans), y(~i_nans), 1);
    predict_x = [min(x), max(x)];
    fit = polyval(lm, predict_x);
    plot(predict_x, fit, 'r:', 'LineWidth', 2.5)  % refline(lm(1), lm(2)) % lsline()
    % metrics
    rmse = sqrt(nanmean((x - y) .^ 2));
    bias = nanmean(x - y);
    title(sprintf('%s\n%.2g (rmse), %.2g (bias), \\color{red}%.1g (r^2adj)',...
        flux, rmse, bias, lm_full.Rsquared.Adjusted))
    xlabel('SCOPE')
    ylabel('bigleaf my')
    refline(1, 0)
end

sgtitle('Bigleaf ebal no Fc vs original SCOPE')
