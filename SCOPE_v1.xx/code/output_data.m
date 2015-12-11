%% OUTPUT DATA
% author C. Van der tol
% modified:      31 Jun 2008: (CvdT) included Pntot in output fluxes.dat
% last modified: 04 Aug 2008: (JT)   included variable output directories
%                31 Jul 2008: (CvdT) added layer_pn.dat
%                19 Sep 2008: (CvdT) spectrum of outgoing radiation
%                19 Sep 2008: (CvdT) Pntot added to fluxes.dat
%                15 Apr 2009: (CvdT) Rn added to vertical profiles
%                03 Oct 2012: (CvdT) included boolean variabel calcebal
%                04 Oct 2012: (CvdT) included reflectance and fPAR
%                10 Mar 2013: (CvdT) major revision: introduced structures
%                22 Nov 2013: (CvdT) added additional outputs
%% Standard output
keyboard
% fluxes
%Output_dir='D:\Dropbox\ARTMOV3\models\3_Combined\2_SCOPE\SCOPE_v1.40_distr\output\';
fidf                = fopen([Output_dir,'fluxes.dat'],'a');
fprintf(fidf,'%9.0f %9.0f %9.0f %9.4f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f',...
    [k iter.counter xyt.year(k) xyt.t(k) fluxes.Rntot fluxes.lEtot fluxes.Htot fluxes.Rnctot fluxes.lEctot, ...
    fluxes.Hctot fluxes.Actot fluxes.Rnstot fluxes.lEstot fluxes.Hstot fluxes.Gtot fluxes.Resp 1E6*fluxes.aPAR  1E6*fluxes.aPAR_Cab fluxes.aPAR/rad.inPAR fluxes.aPAR_Wm2]);
if options.calc_fluor
    fprintf(fidf,'%9.4f %9.6f', rad.Eoutf,  rad.Eoutf./fluxes.aPAR_Wm2);
end
fprintf(fidf,'\r');

% surftemp
fidt                = fopen([Output_dir,'surftemp.dat'],'a');
fprintf(fidt,'%9.0f  %9.0f %9.4f % 9.2f %9.2f %9.2f %9.2f  %9.2f',...
    [k xyt.year(k) xyt.t(k) thermal.Ta thermal.Ts(1) thermal.Ts(2) thermal.Tcave thermal.Tsave]);
fprintf(fidt,'\r');

% aerodyn
fidra               = fopen([Output_dir,'aerodyn.dat'],'a');
fprintf(fidra,'%15.4f %15.4f %15.4f %15.4f',[thermal.raa, thermal.rawc, thermal.raws, thermal.ustar]);
fprintf(fidra,'\r');

% radiation
fidr                = fopen([Output_dir,'radiation.dat'],'a');
fprintf(fidr,'%9.0f  %9.0f %9.4f % 9.2f %9.2f  %9.2f',[k xyt.year(k) xyt.t(k) rad.Eouto rad.Eoutt + rad.Eoutte  rad.Eouto+rad.Eoutt + rad.Eoutte]);
fprintf(fidr,'\r');

% spectrum (added on 19 September 2008)
fidfho         = fopen([Output_dir,'spectrum_hemis_optical.dat'],'a');
fprintf(fidfho,'%9.5f ',rad.Eout_');
fprintf(fidfho,'\r');

fidfoo         = fopen([Output_dir,'spectrum_obsdir_optical.dat'],'a');
fprintf(fidfoo,'%9.5f ',rad.Lo_');
fprintf(fidfoo,'\r');

if options.calc_ebal
    fidto      = fopen([Output_dir,'spectrum_obsdir_BlackBody.dat'],'a');
    fprintf(fidto,'%9.2f', rad.LotBB_');
    fprintf(fidto,'\r');
    
    if options.calc_planck
        fidplanckh      = fopen([Output_dir,'spectrum_hemis_thermal.dat'],'a');
        fprintf(fidplanckh,'%9.2f', rad.Eoutte_');
        fprintf(fidplanckh,'\r');
        
        fidplancko      = fopen([Output_dir,'spectrum_obsdir_thermal.dat'],'a');
        fprintf(fidplancko,'%9.2f', rad.Lot_');
        fprintf(fidplancko,'\r');
    end
end

fidsi           = fopen([Output_dir,'irradiance_spectra.dat'],'a');
fprintf(fidsi,'%10.2f',meteo.Rin*(rad.fEsuno+rad.fEskyo)');
fprintf(fidsi,'\r');

fidref           = fopen([Output_dir,'reflectance.dat'],'a');
reflectance = pi*rad.Lo_./(rad.Esun_+rad.Esky_);

reflectance(spectral.wlS>3000) = NaN;
fprintf(fidref,'%9.5f',reflectance');
fprintf(fidref,'\r');

% wavelength (added on 23 February 2009)
% fidwl          = fopen([Output_dir,'wl.dat'],'a');
% fprintf(fidwl,'%9.5f ',spectral.wlS);
% fprintf(fidwl,' \r');

% input and parameter values (added June 2012)
fidv          = fopen([Output_dir,'pars_and_input.dat'],'a');
for i = 1:length(V)
    fprintf(fidv,'%12.3f',V(i).Val(vi(i)));
end
fprintf(fidv,'\r');

fidvs          = fopen([Output_dir,'pars_and_input_short.dat'],'a');
k2 = find(vmax>1);
for i = 1:length(k2)
    fprintf(fidvs,'%9.5f ',V(k2(i)).Val(vi(k2(i))));
end
fprintf(fidvs,' \r');


%% Optional Output

if options.calc_vert_profiles
    
    % gap
    fidgp               = fopen([Output_dir,'gap.dat'],'a');
    fprintf(fidgp,'%9.2f %9.2f %9.2f',[gap.Ps gap.Po gap.Pso]);
    fprintf(fidgp,'\r');
    
    fidpl           = fopen([Output_dir,'layer_aPAR.dat'],'a');
    fprintf(fidpl,'%9.2f',[1E6*profiles.Pn1d' 0]);
    fprintf(fidpl,'\r');
    
    fidplC          = fopen([Output_dir,'layer_aPAR_Cab.dat'],'a');
    fprintf(fidplC,'%9.2f',[1E6*profiles.Pn1d_Cab' 0]);
    fprintf(fidplC,'\r');
    
    if options.calc_ebal
        
        % leaftemp
        fidtc           = fopen([Output_dir,'leaftemp.dat'],'a');
        fprintf(fidtc,'%9.2f',[profiles.Tcu1d' profiles.Tch' profiles.Tc1d']);
        fprintf(fidtc,'\r');
        
        fidhl           = fopen([Output_dir,'layer_h.dat'],'a');
        fprintf(fidhl,'%9.2f',[profiles.Hc1d' fluxes.Hstot]);
        fprintf(fidhl,'\r');
        
        fidlel          = fopen([Output_dir,'layer_le.dat'],'a');
        fprintf(fidlel,'%9.2f',[profiles.lEc1d' fluxes.lEstot]);
        fprintf(fidlel,'\r');
        
        fidal           = fopen([Output_dir,'layer_a.dat'],'a');
        fprintf(fidal,'%9.2f',[profiles.A1d' fluxes.Resp]);
        fprintf(fidal,'\r');
        
        fidNPQ           = fopen([Output_dir,'layer_NPQ.dat'],'a');
        fprintf(fidNPQ,'%9.2f',[profiles.qE' 0]);
        fprintf(fidNPQ,'\r');
        
        fidrn           = fopen([Output_dir,'layer_rn.dat'],'a');
        fprintf(fidrn,'%9.2f',[profiles.Rn1d' fluxes.Rnstot]);
        fprintf(fidrn,'\r');
    end
    if options.calc_fluor
        fidfll          = fopen([Output_dir,'layer_fluorescence.dat'],'a');
        fprintf(fidfll,'%9.2f',profiles.fluorescence');       
        fprintf(fidfll,'\r');
    end
end

if options.calc_fluor% && options.calc_ebal
    fidfl          = fopen([Output_dir,'fluorescence.dat'],'a');
    fidfl1         = fopen([Output_dir,'fluorescencePSI.dat'],'a');
    fidfl2         = fopen([Output_dir,'fluorescencePSII.dat'],'a');
    fidflh         = fopen([Output_dir,'fluorescence_hemis.dat'],'a');
    for j=1:size(spectral.wlF,1)
        fprintf(fidfl,'%10.4f ',rad.LoF_);
        fprintf(fidfl1,'%10.4f ',rad.LoF1_);
        fprintf(fidfl2,'%10.4f ',rad.LoF2_);
        fprintf(fidflh,'%10.4f ',rad.Fhem_);
    end
    fprintf(fidfl,' \r');
    fprintf(fidfl1,' \r');
    fprintf(fidfl2,' \r');
    fprintf(fidflh,' \r');
end

fidp                = fopen([Output_dir,'BOC_irradiance.dat'],'a');
fprintf(fidp,'%9.0f', rad.Emin_(61,:)+(rad.Esun_*gap.Ps(61)')');
fprintf(fidp,'\r');

%%
if options.calc_directional && options.calc_ebal
    Output_angle    =   [directional.tto';  directional.psi'; angles.tts*ones(size(directional.psi'))];
    Output_brdf     =   [spectral.wlS'     directional.brdf_];
    if options.calc_planck
        Output_temp =   [spectral.wlT'    directional.Lot_];
    else
        Output_temp =   [directional.BrightnessT];
    end
    if options.calc_fluor
        Output_fluor = [spectral.wlF'     directional.LoF_];
    end
    
    save([Output_dir,'Directional/',sprintf('BRDF (SunAngle %2.2f degrees).dat',angles.tts)],'Output_brdf' ,'-ASCII','-TABS')
    save([Output_dir,'Directional/',sprintf('Angles (SunAngle %2.2f degrees).dat',angles.tts)],'Output_angle','-ASCII','-TABS')
    save([Output_dir,'Directional/',sprintf('Temperatures (SunAngle %2.2f degrees).dat',angles.tts)],'Output_temp','-ASCII','-TABS')
    
    if options.calc_fluor
        save([Output_dir,'Directional/',sprintf('Fluorescence (SunAngle %2.2f degrees).dat',angles.tts)],'Output_fluor','-ASCII','-TABS')
    end
    
    fiddirtir       =   fopen([Output_dir,'Directional/','read me.txt'],'w');
    fprintf(fiddirtir,'The Directional data is written in three files: \r\n');
    fprintf(fiddirtir,'\r\n- Angles: contains the directions. \r\n');
    fprintf(fiddirtir,'   * The 1st row gives the observation zenith  angles\r\n');
    fprintf(fiddirtir,'   * The 2nd row gives the observation azimuth angles\r\n');
    fprintf(fiddirtir,'   * The 3rd row gives the solar       zenith  angles\r\n');
    fprintf(fiddirtir,'\r\n- Temperatures: contains the directional brightness temperature. \r\n');
    fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the brightness temperature values (except for broadband)\r\n');
    fprintf(fiddirtir,'   * The 2nd column gives the Tb values corresponding to the directions given by first column in the Angles file\r\n');
    fprintf(fiddirtir,'\r\n- BRDF: contains the bidirectional distribution functions values. \r\n');
    fprintf(fiddirtir,'   * The 1st column gives the wl values corresponding to the BRDF values\r\n');
    fprintf(fiddirtir,'   * The 2nd column gives the BRDF values corresponding to the directions given by first column in the Angles file\r\n');
    fclose(fiddirtir);
end

%%
fclose all;