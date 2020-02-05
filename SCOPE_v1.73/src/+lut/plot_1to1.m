function plot_1to1(meas, mod, what, out_path)
    figure
%     errorbar(val_out, res, res_std, 'o')
    scatter(meas, mod)
    xlabel('SCOPE GPP [\mumol CO_2 m^{-2} s^{-1}]')
    ylabel([what ' GPP [\mumol CO_2 m^{-2} s^{-1}]'])
    l = refline(1, 0);
    l.Color = 'r';
    refline
    title([what ' vs SCOPE'])
    legend('scatter', '1-to-1 line', 'regression', 'Location', 'southeast')
    saveas(gcf, out_path)
end