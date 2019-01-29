N=[
1;		%calc_ebal		calculate the complete energy balance
0;		%calc_vert_profiles		calculate vertical profiles of fluxes and temperatures
1;		%calc_fluor		calculate chlorophyll fluorescence
0;		%calc_planck		calculate spectrum of thermal radiation with spectral emissivity instead of broadband
0;		%calc_directional		calculate BRDF and directional temperature for many angles specified in a file. Be patient, this takes some time
1;      %calc_xanthophyllabs		calculate dynamic xanthopyll absorption (zeaxanthin)
0;      %calc_PSI       0 (recommended): treat the whole fluorescence spectrum as one spectrum (new calibrated optipar), 1: differentiate PSI and PSII with Franck et al. spectra (of SCOPE 1.62 and older)
0;		%rt_thermal		0: provide emissivity values as input. 1: use values from fluspect and soil at 2400 nm for the TIR range
0;		%calc_zo		0: use the zo and d values provided in the inputdata, 1: calculate zo and d from the LAI, canopy height, CD1, CR, CSSOIL (recommended if LAI changes in time series)
0;      %0: use soil spectrum from a file, 1: simulate soil spectrum with the BSM model
0;		%SoilHeatMethod		0: standard calculation of thermal inertia from soil characteristics, 1: empiricaly calibrated formula (make function), 2: as constant fraction of soil net radiation
1;		%Fluorescence_model		0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
0;		%calc_rss_rbs		0: use resistance rss and rbs as provided in inputdata. 1:  calculate rss from soil moisture content and correct rbs for LAI (calc_rssrbs.m)
1;		%applTcorr		correct Vcmax and rate constants for temperature in biochemical.m
1;		%verify		verifiy the results (compare to saved 'standard' output) to test the code for the firstt ime
1;		%saveheaders		write header lines in output files
0;		%makeplots 		plot the results
1];		%simulation		0: individual runs. Specify one value for constant input, and an equal number (>1) of values for all input that varies between the runs.
		%		1: time series (uses text files with meteo input as time series)
		%		2: Lookup-Table (specify the values to be included. All possible combinations of inputs will be  used)
