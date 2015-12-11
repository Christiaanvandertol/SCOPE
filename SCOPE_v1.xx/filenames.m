%The following three are always required,

X = {
'Simulation_Name'	, 'verificationdata';
'soil_file'		, 'soilnew.txt';
'leaf_file'		, 'optipar_fluspect_2014.txt'; 
'atmos_file' 		, 'FLEX-S3_std.atm';

%The following are only for the time series option!
'Dataset_dir'		, 'for_verification'; 
't_file'		, 't_.dat';
'year_file'		, 'year_.dat';
'Rin_file'		, 'Rin_.dat';
'Rli_file'		, 'Rli_.dat';
'p_file' 		, 'p_.dat';
'Ta_file'		, 'Ta_.dat';
'ea_file'		, 'ea_.dat';
'u_file'		, 'u_.dat';

%optional (leave empty for constant values From inputdata.TXT)
'CO2_file'		, '';
'z_file' 		, '';
'tts_file' 		, '';

%optional two column tables (first column DOY second column value)
'LAI_file'		, '';
'hc_file'		, '';
'SMC_file'		, '';
'Vcmax_file'		, '';
'Cab_file'		, ''};
