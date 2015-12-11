function [soil,leafbio,canopy,meteo,angles,xyt] = select_input(V,vi,canopy,options,xyt,soil)

soil.spectrum= V(16).Val(vi(16));
soil.rss = V(17).Val(vi(17));
soil.rs_thermal = V(18).Val(vi(18));
soil.cs = V(19).Val(vi(19));
soil.rhos = V(20).Val(vi(20));
soil.CSSOIL  = V(43).Val(vi(43));
soil.lambdas = V(21).Val(vi(21));
soil.rbs  = V(44).Val(vi(44));
soil.SMC  = V(54).Val(vi(54));

leafbio.Cab = V(1).Val(vi(1));
leafbio.Cca = V(2).Val(vi(2));
if options.Cca_function_of_Cab
    leafbio.Cca = 0.25*V(1).Val(vi(1));
end
leafbio.Cdm  = V(3).Val(vi(3));
leafbio.Cw = V(4).Val(vi(4));
leafbio.Cs   = V(5).Val(vi(5));
leafbio.N = V(6).Val(vi(6));
leafbio.Vcmo  = V(9).Val(vi(9));
leafbio.m  = V(10).Val(vi(10));
leafbio.Type = V(11).Val(vi(11));
leafbio.Tparam  = V(14).Val(:); % this is correct (: instead of 14)
leafbio.fqe(2)   = V(15).Val(vi(15));
leafbio.Rdparam  = V(13).Val(vi(13));

leafbio.rho_thermal   = V(7).Val(vi(7));
leafbio.tau_thermal  = V(8).Val(vi(8));

leafbio.Tyear = V(55).Val(vi(55));
leafbio.beta   = V(56).Val(vi(56));
leafbio.kNPQs   = V(57).Val(vi(57));
leafbio.qLs  = V(58).Val(vi(58));
leafbio.stressfactor  = V(59).Val(vi(59));


canopy.LAI  = V(22).Val(vi(22));
canopy.hc  = V(23).Val(vi(23));
canopy.LIDFa = V(26).Val(vi(26));
canopy.LIDFb  = V(27).Val(vi(26)); % this is correct (26 instead of 27)
canopy.leafwidth  = V(28).Val(vi(28));
canopy.rb   = V(38).Val(vi(38));
canopy.Cd  = V(39).Val(vi(39));
canopy.CR = V(40).Val(vi(40));
canopy.CD1  = V(41).Val(vi(41));
canopy.Psicor  = V(42).Val(vi(42));
canopy.rwc  = V(45).Val(vi(45));
canopy.kV = V(12).Val(vi(12));
canopy.zo  = V(24).Val(vi(24));
canopy.d = V(25).Val(vi(25));

meteo.z  = V(29).Val(vi(29));
meteo.Rin   = V(30).Val(vi(30));
meteo.Ta = V(31).Val(vi(31));
meteo.Rli  = V(32).Val(vi(32));
meteo.p  = V(33).Val(vi(33));
meteo.ea  = V(34).Val(vi(34));
meteo.u   = V(35).Val(vi(35));
meteo.Ca = V(36).Val(vi(36));
meteo.Oa  = V(37).Val(vi(37));

xyt.startDOY = V(46).Val(vi(46));
xyt.endDOY = V(47).Val(vi(47));
xyt.LAT = V(48).Val(vi(48));
xyt.LON = V(49).Val(vi(49));
xyt.timezn = V(50).Val(vi(50));
angles.tts = V(51).Val(vi(51));
angles.tto = V(52).Val(vi(52));
angles.psi = V(53).Val(vi(53));

%% derived input
if options.soil_heat_method ==1
    soil.GAM =  Soil_Inertia1(soil.SMC);
else
    soil.GAM  = Soil_Inertia0(soil.cs,soil.rhos,soil.lambdas);
end
if options.calc_rss_rbs
    [soil.rss,soil.rbs] = calc_rssrbs(soil.SMC,canopy.LAI,soil.rbs);
end

if leafbio.Type,
    leafbio.Type = 'C4';
else
    leafbio.Type = 'C3';
end
canopy.hot  = canopy.leafwidth/canopy.hc;
if options.calc_zo
    [canopy.zo,canopy.d ]  = zo_and_d(soil,canopy);
end
leafbio.fqe(1)   = leafbio.fqe(2)/5;

