function plot_image(im, what, out_path)
    figure
    imagesc(im)
    cb = colorbar;
    title(cb, '[\mumol CO_2 cm^{-2} s^{-1}]')
    title([what ' GPP'])
    saveas(gcf, out_path)
end