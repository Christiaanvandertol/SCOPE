function [spectral] = define_bands

    % Define spectral regions for SCOPE v_1.40
    % All spectral regions are defined here as row vectors
    % WV Jan. 2013
    
    % 3 spectral regions for SCOPE
    
    reg1 =   400 :    1 :  2400;
    reg2 =  2500 :  100 : 15000;
    reg3 = 16000 : 1000 : 50000;
    
    spectral.wlS  = [reg1 reg2 reg3];
        
    % Other spectral (sub)regions
    
    spectral.wlP   = reg1;                            % PROSPECT data range
    spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
    spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
    spectral.wlO   = reg1;                            % optical part
    spectral.wlT   = [reg2 reg3];                     % thermal part
    spectral.wlZ   = 500:1:600;                       % xanthophyll region
    wlS            = spectral.wlS;
    spectral.wlPAR = wlS(wlS>=400 & wlS<=700);  % PAR range
    
    % Data used by aggreg routine to read in MODTRAN data
    
    spectral.SCOPEspec.nreg = 3;
    spectral.SCOPEspec.start = [ 400  2500  16000];
    spectral.SCOPEspec.end   = [2400 15000  50000];
    spectral.SCOPEspec.res   = [   1   100   1000];
    
end
