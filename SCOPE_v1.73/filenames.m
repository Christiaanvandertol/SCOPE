%The following three are always required,

X = {
'Simulation_Name'	, 'verification_run';
'soil_file'		, 'soilnew.txt';
'leaf_file'		, 'Optipar2017_ProspectD.mat'; 
'atmos_file' 		, 'FLEX-S3_std.atm';

%The following are only for the time series option!
'Dataset_dir'		, 'for_verification'; 
'verification_dir'	, 'verificationdata';
'meteo_ec_csv'	, 'ts_input.csv';

% values from `vegetation_retrieved_csv` will be linearly interpolated to timestamp (t) of meteo_ec_csv
'vegetation_retrieved_csv'	, '';

% Variables below are COLUMN NAMES in time_series_files (meteo_ec_csv, vegetation_retrieved_csv)
't'		, 'f_date';
'Rin'   , 'Rin';
'Rli'   , 'Rli';
'p' 	, 'p';
'Ta'	, 'Ta';
'u'		, 'u';
'ea'	, 'ea';
'RH'	, '';  % only for ea calculations ea = RH * satvap(Ta)

'tts' 	, '';  % if not present - calculated from t, timezn, LAT, LON

%optional (leave empty for constant values From inputdata.TXT)
% leaf
'Cab' 		, '';
'Cca'		, '';
'Cdm'		, '';
'Cw'	    , '';
'Cs'		, '';
'Cant'		, '';
'N'		    , '';

% BSM
'SMC'				, '';
'BSMBrightness'		, '';
'BSMlat'			, '';
'BSMlon'		    , '';

% canopy
'LAI'		, '';  
'hc'		, '';
'LIDFa'		, '';
'LIDFb'		, '';

% meteo
'z'		, '';
'Ca'	, '';

%biochemistry
'Vcmo'		, '';


%optional leaf inclination distribution file with 13 rows (see
%example). It MUST be located in ../data/leafangles/
'LIDF_file'     , ''};