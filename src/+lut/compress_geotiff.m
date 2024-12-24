function compress_geotiff(im_io_path, im_crs_path)

    ver_out = ver;
    toolboxes = {ver_out.Name};
    geo_tif = any(strcmp('Mapping Toolbox', toolboxes));
    
    if geo_tif
        fprintf('copying georeference tags from input image %s\n', im_crs_path)
        geoinfo = geotiffinfo(im_crs_path);
        key = geoinfo.GeoTIFFTags.GeoKeyDirectoryTag;
        R = geoinfo.SpatialRef;
    else
        warning(['Mapping Toolbox is not installed. Output .tifs can not be georeferenced.\n'...
            'Use gdal_translate to georeference from input image %s\n'], im_crs_path)
        return
    end
    
    im = imread(im_io_path);
    comp.Compression = Tiff.Compression.PackBits;
    geotiffwrite(im_io_path, im, R, 'GeoKeyDirectoryTag', key, 'TiffTags', comp)
    
end