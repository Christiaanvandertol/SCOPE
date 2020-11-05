function write_tiff(im, out_path)
    
    input_image_path = '../exercise/images/Cab.tif';
    copyfile(input_image_path, out_path)
    t = Tiff(out_path, 'r+');
    if getTag(t, 'BitsPerSample') == 32
        im = single(im);
    end
    write(t, im)
    close(t)

end