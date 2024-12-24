function write_tiff(im, out_path, im_crs_path)
    
    copyfile(im_crs_path, out_path)
    t = Tiff(out_path, 'r+');
    if getTag(t, 'BitsPerSample') == 32
        im = single(im);
    end
    write(t, im)
    close(t)

    lut.compress_geotiff(out_path, im_crs_path)

end