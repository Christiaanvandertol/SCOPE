function plot_image(im, what, out_path)
    figure
    imagesc(im)
    cb = colorbar;
    title(cb, '[\mumol CO_2 m^{-2} s^{-1}]')
    title(['GPP with ' what])
    % if GPR is out of range (for SNAP), highly negative values appear
    caxis([quantile(im(:), 0.01), quantile(im(:), 0.99)])  
    saveas(gcf, out_path)
end