 function biochem_out = biochemical(biochem_in)
%
% Date: 	21 Sep 2012
% Update:   20 Feb 2013
% Update:      Aug 2013: correction of L171: Ci = Ci*1e6 ./ p .* 1E3;
%   
% Authors: 	Joe Berry and Christiaan van der Tol, contributions of others.
% Sources: 	
%           Farquhar et al. 1980, Collatz et al (1991, 1992).
%
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1)
%    - photosynthesis of a leaf or needle (umol m-2 s-1)
%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
%
% Usage:
% biochem_out = biochemical(biochem_in)
% the function was tested for Matlab 7.2.0.232 (R2006a)
%
% Calculates net assimilation rate A, fluorescence F using biochemical model
%
% Input (units are important):
% structure 'biochem_in' with the following elements:
%   Fluorescence_model          integer with
%                               0: with sustained quenching after drought, as in Lee et al. (2013) 
%                               1: calibrated for cotton data set: no drought
% Cs        % [umol m-3]            initial estimate of conc. of CO2 in the
%                                   ...bounary layer of the leaf
% Q         % [umol photons m-2 s-1]net radiation, PAR
% T         % [oC or K]             leaf temperature
% eb        % [hPa]                 intial estimate of the vapour pressure in leaf boundary layer
% O         % [mmol m-3]            concentration of O2 (in the boundary
%                                   ...layer, but no problem to use ambient)
% p         % [hPa]                 air pressure
% Vcmo      % [umol/m2/s]           maximum carboxylation capacity
% m         % []                    Ball-Berry coefficient 'm' for stomatal regulation
% Type      % []                    text parameter, either 'C3' for C3 or any
%                                   ...other text for C4
% tempcor   % []                    boolean (0 or 1) whether or not
%                                   ...temperature correction to Vcmax has to be applied.
% Tsparams  % [],[],[K],[K],[K]     vector of 5 temperature correction
%                                   parameters, look in spreadsheet of PFTs. Only if tempcor=1, otherwise use
%                                   dummy values
% Rdparam   % []                    respiration as fraction of Vcmax
% stressfactor []                   optional input: stress factor to reduce Vcmax (for
%                                   example soil moisture, leaf age). Default value = 1.

%
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional
% matrices
%
% Output:
% structure 'biochem_out' with the following elements:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Cs        % [umol/m3]             CO2 concentration in the boundary layer
% eta0      % []                    fluorescence as fraction of dark
%                                   ...adapted (fs/fo)
% rcw       % [s m-1]               stomatal resistance
% qE        % []                    non photochemical quenching
% fs        % []                    fluorescence as fraction of PAR
% Ci        % [umol/m3]             internal CO2 concentration
% Kn        % []                    rate constant for excess heat
% fo        % []                    dark adapted fluorescence (fraction of aPAR)
% fm        % []                    light saturated fluorescence (fraction of aPAR)
% qQ        % []                    photochemical quenching
% Vcmax     % [umol/m2/s]           carboxylation capacity after
%                                   ... temperature correction 

%% input
Rdparam       = biochem_in.Rdparam;                          %[umol/m2/s]        dark respiration rate at 25 oC as fraction of Vcmax
Tparams       = biochem_in.Tparams;       %[]                 temperature sensitivities of Vcmax, etc (dummy values of applTcorr was selected above)
Cs            = biochem_in.Cs;
Q             = biochem_in.Q;
T             = biochem_in.T;
eb            = biochem_in.eb;
O             = biochem_in.O;
p             = biochem_in.p;
Vcmo          = biochem_in.Vcmo;
m             = biochem_in.m;
Type          = biochem_in.Type;
tempcor       = biochem_in.tempcor;
stressfactor  = biochem_in.stressfactor;
model_choice  = biochem_in.Fluorescence_model;

T           = T+273.15*(T<100); % convert temperatures to K if not already

rhoa        = 1.2047;           % [kg m-3]       specific mass of air
Mair        = 28.96;            % [g mol-1]      molecular mass of dry air

%% parameters (at optimum temperature)
Kcopt       = 350;              % [ubar]        kinetic coefficient for CO2 (Von Caemmerer and Furbank, 1999)
Koopt       = 450;              % [mbar]        kinetic coeeficient for  O2 (Von Caemmerer and Furbank, 1999)
Kf          = 0.05;             % []            rate constant for fluorescence
%Kd          = 0.95;             % []            rate constant for thermal deactivation at Fm
Kp          = 4.0;              % []            rate constant for photochemisty
kpopt       = Vcmo/56*1E6;      % []            PEPcase rate constant for C02, used here: Collatz et al: Vcmo = 39 umol m-1 s-1; kp = 0.7 mol m-1 s-1.
if isfield(biochem_in,'atheta') 
    atheta = biochem_in.atheta;
else
    atheta      = 0.8;
end
if biochem_in.Fluorescence_model==0
    % default drought values: 
    Knparams = [5.01, 1.93, 10];
else
    % default general values (cotton dataset)
    Knparams = [2.48, 2.83, 0.114];    
end

%% temperature definitions
Tref        = 25+273.15;        % [K]           absolute temperature at 25 oC
slti        = Tparams(1);
%shti        = Tparams(2);
Thl         = Tparams(3);
%Thh         = Tparams(4);
Trdm        = Tparams(5);

%% convert all to bar
Cs          = Cs * 1e-6 .* p .*1E-3;
O           = O * 1e-3 .* p .*1E-3 * ~strcmp('C4',Type);       % forced to be zero for C4 vegetation (this is a trick to prevent oxygenase)
Kcopt       = Kcopt * 1e-6;
Koopt       = Koopt * 1e-3;

%% temperature corrections
qt          = 0.1 * (T-Tref) * tempcor;
%TH          = 1 + tempcor* exp(shti .* (T   -Thh));
TH          = 1 + tempcor* exp((-220E3+703*T)./(8.314*T));
TL          = 1 + tempcor* exp(slti .* (Thl -T));

Kc          = Kcopt * 2.1.^qt;
Ko          = Koopt * 1.2.^qt;
kp          = kpopt.* 1.8.^qt;
Kd          = max(0.0301*(T-273.15)+ 0.0773,0.8738);

Rd          = Rdparam * Vcmo .* 1.8.^qt ./(1+exp(1.3*(T-Trdm)));
switch Type
    case 'C3'
        Vcmax       =           Vcmo .* 2.1.^qt ./TH * stressfactor;
    otherwise 
        Vcmax       =           Vcmo .* 2.1.^qt ./(TL.*TH) * stressfactor;
end
spfy        = 2600 * 0.75 .^qt;            % This is, in theory, Vcmax/Vomax.*Ko./Kc, but used as a separate parameter

%% calculation of potential electron transport rate
po0         = Kp./(Kf+Kd+Kp);         % dark photochemistry fraction (Genty et al., 1989)
Je          = 0.5*po0 .* Q;          % electron transport rate

%% calculation of the intersection of enzyme and light limited curves
% this is the original Farquhar model
gam         = 0.5 ./spfy .*O; %[bar]       compensation point [bar]

%% calculation of internal CO2 concentration, photosynthesis
RH          = eb./satvap(T-273.15);  
% a = 5;
% D0 = 5;
% gs0 = 1E-6;
    
switch Type
    case 'C3'
        Ci          = max(.3*Cs,Cs.*(1-1./(m.*RH)));
        Vc          = Vcmax.*(Ci-gam)./((Kc .* (1+O./Ko)) + Ci);
        Vs          = Vcmo/2 .* 1.8.^qt;    %
        effcon      = 0.2;
    otherwise                   %C4
        %if biochem_in.A == -999
            Ci          = max(.1*Cs,Cs.*(1-1./(m.*RH)));
        %else
        %    gs   = gs0 + a*max(0,biochem_in.A)./((Cs-gam) .* (1+(1-RH).*satvap(T-273.15)./D0));
        %    Ci = Cs - biochem_in.A./gs;
        %end

        Vc          = Vcmax;
        Vs          = kp.*Ci;
        effcon      = 0.17;                    % Berry and Farquhar (1978): 1/0.167
end
Ve          = Je.*(Ci-gam)./(Ci+2*gam) .* effcon;

[a1,a2]     = abc(atheta,-(Vc+Ve),Vc.*Ve);
V           = min(a1,a2).*(Ci>gam) + max(a1,a2).*(Ci<=gam);
[a1,a2]     = abc(0.98,-(V+Vs),V.*Vs);
Ag          = min(a1,a2);
A           = Ag - Rd;
Ja          = Ag ./((Ci-gam)./(Ci+2*gam))./effcon;        % actual electron transport rate

rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3    * 1e6 ./ p .* 1E3;
rcw(A<=0)   = 0.625*1E6;

%% fluorescence (Replace this part by Magnani or other model if needed
ps          = po0.*Ja./Je;               % this is the photochemical yield
ps(isnan(ps)) = 0.8;
[eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = Fluorescencemodel(ps,Kp,Kf,Kd,Knparams);
Kpa         = ps./fs*Kf;

%% convert back to ppm
Ci          = Ci*1e6 ./ p .* 1E3;

%% Collect outputs

biochem_out.A       = A;
biochem_out.Ci      = Ci;
biochem_out.eta     = eta;
biochem_out.rcw     = rcw;
biochem_out.qE      = qE;
biochem_out.fs      = fs;
biochem_out.Kn      = Kn;
biochem_out.fo0     = fo0;
biochem_out.fm0     = fm0;
biochem_out.fo      = fo;
biochem_out.fm      = fm;
biochem_out.qQ      = qQ;
biochem_out.Vcmax   = Vcmax;
biochem_out.Kp      = Kpa;
biochem_out.ps      = ps;
biochem_out.Ja      = Ja;
biochem_in.A        = A;
return;
%%% end of function biochemical

%% abc formula
function [x2,x1] = abc(a,b,c)
if a == 0
    x1      = -c./b;
    x2      = x1;
else
    x1      = (-b+sqrt(b.^2-4.*a.*c))./(2.*a);
    x2      = (-b-sqrt(b.^2-4.*a.*c))./(2.*a);
end
return;
%%% end of abc formula

function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = Fluorescencemodel(ps,Kp,Kf,Kd,Knparams)
po0         = Kp./(Kf+Kd+Kp);         % dark photochemistry fraction (Genty et al., 1989)
x           = 1-ps./po0 ;                % degree of light saturation

Kno         = Knparams(1);
alpha       = Knparams(2);
beta        = Knparams(3);

% using exp(-beta) expands the interesting region between 0-1
%beta = exp(-beta);
x_alpha     = x.^alpha; % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
Kn          = Kno * (1+beta).* x_alpha./(beta + x_alpha);

fo0         = Kf./(Kf+Kp+Kd);           % dark adapted fluorescence yield Fo
fo          = Kf./(Kf+Kp+Kd+Kn);           % dark adapted fluorescence yield Fo
fm          = Kf./(Kf   +Kd+Kn);        % light adapted fluorescence yield Fm
fm0         = Kf./(Kf   +Kd);        % light adapted fluorescence yield Fm
fs          = fm.*(1-ps);
eta         = fs./fo0;
qQ          = 1-(fs-fo)./(fm-fo);       % photochemical quenching
qE          = 1-(fm-fo)./(fm0-fo0);     % non-photochemical quenching

eta         = (1+5)/5*eta-1/5;

return