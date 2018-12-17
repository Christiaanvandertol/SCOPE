function [V,xyt,canopy]  = load_timeseries(V,leafbio,soil,canopy,meteo,constants,F,xyt,path_input,options)

Dataset_dir         = ['dataset ' char(F(5).FileName)];
t_file              = char(F(6).FileName);
year_file           = char(F(7).FileName);
Rin_file            = char(F(8).FileName);
Rli_file            = char(F(9).FileName); 
p_file              = char(F(10).FileName); 
Ta_file             = char(F(11).FileName); 
ea_file             = char(F(12).FileName); 
u_file              = char(F(13).FileName); 
CO2_file            = char(F(14).FileName); 
z_file              = char(F(15).FileName); 
tts_file            = char(F(16).FileName);
LAI_file            = char(F(17).FileName);
hc_file             = char(F(18).FileName);
SMC_file            = char(F(19).FileName);
Vcmax_file          = char(F(20).FileName);
Cab_file            = char(F(21).FileName);

%% 1. Time and zenith angle
xyt.t               = load([path_input,Dataset_dir,'/' ,t_file] );
xyt.year            = load([path_input,Dataset_dir,'/',year_file]);
t_                  = xyt.t;

DOY_                = floor(t_);
time_               = 24*(t_-DOY_);

if ~isempty(tts_file)
    V(51).Val      = load([path_input,Dataset_dir,'/',tts_file]);
else
    ttsR            = calczenithangle(DOY_,time_ - xyt.timezn ,0,0,xyt.LON,xyt.LAT);     %sun zenith angle in rad
    V(51).Val       = min(85,ttsR/constants.deg2rad);                         %sun zenith angle in deg
end
%% 2. Radiation 
if ~isempty(Rin_file)
    V(30).Val           = load([path_input,Dataset_dir,'/',Rin_file]);
else 
    V(30).Val       = V(30).Val*ones(size(t_));
end
if ~isempty(Rli_file)
    V(32).Val           = load([path_input,Dataset_dir,'/',Rli_file]);
else
    V(32).Val       = V(32).Val*ones(size(t_));
end

%% 3. Windspeed, air temperature, humidity and air pressure
if ~isempty(u_file)% wind speed
    V(35).Val           = load([path_input,Dataset_dir,'/',u_file]); 
else
    V(35).Val           = V(35).Val*ones(size(t_));
end

if ~isempty(Ta_file)%air temperature
    V(31).Val           = load([path_input,Dataset_dir,'/',Ta_file]);
else
    V(31).Val           = V(31).Val*ones(size(t_));
end

if ~isempty(ea_file)%air temperature
    V(34).Val           = load([path_input,Dataset_dir,'/',ea_file]);
else
    V(34).Val           = V(34).Val*ones(size(t_));
end

if ~isempty(p_file)
    V(33).Val       = load([path_input,Dataset_dir,'/',p_file]);
else
    V(33).Val       = V(33).Val*ones(size(t_));
end

%% 4. Vegetation structure (measurement height, vegetation height and LAI)
if ~isempty(z_file)
    ztable          = load([path_input,Dataset_dir,'/',z_file]);
    V(29).Val       = interp1(ztable(:,1),ztable(:,2),t_);
else
    V(29).Val       = meteo.z*ones(size(t_));
end
if  ~isempty(LAI_file)
    LAItable        = load([path_input,Dataset_dir,'/',LAI_file]);
    V(22).Val         = interp1(LAItable(:,1),LAItable(:,2),t_);
else
    V(22).Val          = canopy.LAI*ones(size(time_));
end
if  ~isempty(hc_file)
    hctable         = load([path_input,Dataset_dir,'/',hc_file]);
    V(23).Val       = interp1(hctable(:,1),hctable(:,2),t_); 
    canopy.hc = V(23).Val;
    if options.calc_zo
        [V(24).Val ,V(25).Val ]  = zo_and_d(soil,canopy);
    else
        V(24).Val   = ones(size(t_))*V(24).Val;
        V(25).Val   = ones(size(t_))*V(25).Val;
    end
    
else
    V(23).Val        = canopy.hc*ones(size(t_)); 
    V(24).Val        = canopy.zo*ones(size(t_)); 
    V(25).Val        = canopy.d*ones(size(t_));
end


%% 5. Gas concentrations
if ~isempty(CO2_file)
    Ca_          = load([path_input,Dataset_dir,'/',CO2_file])*constants.Mair/constants.MCO2/constants.rhoa; % conversion from mg m-3 to ppm
    % mg(CO2)/m-3 * g(air)/mol(air) * mol(CO2)/g(CO2) * m3(air)/kg(air) * 10^-3 g(CO2)/mg(CO2) * 10^-3 kg(air)/g(air) * 10^6 ppm  
    jj              = isnan(Ca_);                           %find data with good quality  Ca data
    Ca_(jj)      = 380;
else
    Ca_ = ones(length(t_),1)* 380;
end
V(36).Val       = Ca_;

%% 6. Soil Moisture Content
if ~isempty(SMC_file)
    V(54).Val          = load([path_input,Dataset_dir,'/',SMC_file]);
end

%% 7. Leaf biochemical parameters
if ~isempty(Vcmax_file)
    Vcmaxtable      = load([path_input,Dataset_dir,'/',Vcmax_file]);
    V(9).Val         = interp1(Vcmaxtable(:,1),Vcmaxtable(:,2),t_);
else
    V(9).Val        = leafbio.Vcmo*ones(size(t_));
end

if ~isempty(Cab_file)
    Cabtable        = load([path_input,Dataset_dir,'/',Cab_file]);
    V(1).Val         = interp1(Cabtable(:,1),Cabtable(:,2),t_);
else V(1).Val        = leafbio.Cab*ones(size(t_));
end



