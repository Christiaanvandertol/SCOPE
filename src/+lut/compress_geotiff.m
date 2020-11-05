function compress_geotiff(path_in, path_out)

    ver_out = ver;
    toolboxes = {ver_out.Name};
    geo_tif = any(strcmp('Mapping Toolbox', toolboxes));
    
    if geo_tif
        fprintf('copying georeference tags from input image %s\n', path_in)
        geoinfo = geotiffinfo(path_in);
        key = geoinfo.GeoTIFFTags.GeoKeyDirectoryTag;
        R = geoinfo.SpatialRef;
    else
        warning(['Mapping Toolbox is not installed. Output .tifs can not be georeferenced.\n'...
            'Use gdal_translate to georeference from input image %s\n'...
            'Output .tifs are identical to input image'], path_in)
    end
    
    im = imread(path_in);
    comp.Compression = Tiff.Compression.PackBits;
    geotiffwrite(path_out, im, R, 'GeoKeyDirectoryTag', key, 'TiffTags', comp)
    
end