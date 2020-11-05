function plot_1to1(meas, mod, what, out_path, var_name)
    if nargin < 5
        var_name = 'GPP';
    end

    figure
%     errorbar(val_out, res, res_std, 'o')
    units = '[W m^{-2}]';
    if any(strcmp(var_name, {'Actot', 'GPP'}))
        units = '[\mumol CO_2 m^{-2} s^{-1}]';
    end
    scatter(meas, mod)
    xlabel(sprintf('SCOPE %s %s', var_name, units))
    ylabel(sprintf('%s %s %s', what, var_name, units))
    l = refline(1, 0);
    l.Color = 'r';
    refline
    rmse = sqrt(mean((meas - mod) .^ 2));
    bias = mean(mod-meas);
    lm = fitlm(meas, mod); 
    title(sprintf('%s, RMSE = %.2f, bias=%.2f, R^2_{adj} = %.2f', var_name, rmse, bias, lm.Rsquared.Adjusted))
    legend('scatter', '1-to-1 line', 'regression', 'Location', 'southeast')
    saveas(gcf, out_path)
end