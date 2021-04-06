function savebrdfoutput(options,directional,angles,spectral,Output_dir)
%% directional (This is never saved as binary, but always as ASCII).

if ~exist([Output_dir '\directional'], 'dir')
    mkdir([Output_dir '\directional'])
end
    
Output_angle    =   [directional.tto';  directional.psi']; %#ok<*NASGU>
Output_refl     =   [spectral.wlS    directional.refl_];
Output_rso      =   [spectral.wlS     directional.rso_];
if options.calc_planck
    Output_rad      =   [spectral.wlT'   directional.Lot_];
end
if options.calc_fluor
    Output_fluor = [spectral.wlF'     directional.LoF_];
end

save([Output_dir,'Directional/',sprintf('refl (SunAngle %2.2f degrees).dat',angles.tts)],'Output_refl' ,'-ASCII','-TABS')
save([Output_dir,'Directional/',sprintf('rso (SunAngle %2.2f degrees).dat',angles.tts)],'Output_rso' ,'-ASCII','-TABS')
save([Output_dir,'Directional/',sprintf('Angles (SunAngle %2.2f degrees).dat',angles.tts)],'Output_angle','-ASCII','-TABS')
if options.calc_planck
    save([Output_dir,'Directional/',sprintf('Thermal radiances (SunAngle %2.2f degrees).dat',angles.tts)],'Output_rad','-ASCII','-TABS')
end
if options.calc_fluor
    save([Output_dir,'Directional/',sprintf('Fluorescence (SunAngle %2.2f degrees).dat',angles.tts)],'Output_fluor','-ASCII','-TABS')
end

fiddirtir       =   fopen([Output_dir,'Directional/','read me.txt'],'w');
fprintf(fiddirtir,'The Directional data is written in three files: \r\n');
fprintf(fiddirtir,'\r\n- Angles: contains the observatioin angles. \r\n');
fprintf(fiddirtir,'   * The 1st row gives the observation zenith  angles\r\n');
fprintf(fiddirtir,'   * The 2nd row gives the observation azimuth angles\r\n');
fprintf(fiddirtir,'\r\n- Thermal radiances: contains the radiance at spectral.wlT, emitted plus reflected incoming \r\n');
fprintf(fiddirtir,'   * The 1st column gives the wl values \r\n');
fprintf(fiddirtir,'   * The 2nd column gives the radiances corresponding to the directions given by first column in the Angles file\r\n');
fprintf(fiddirtir,'\r\n- refl and rso: contains the bidirectional distribution functions values, reflectance (based on rso and rdo) and rso. \r\n');
fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the BRDF values\r\n');
fprintf(fiddirtir,'   * The 2nd column gives the reflectance values corresponding to the directions given by first column in the Angles file\r\n');
fclose(fiddirtir);

