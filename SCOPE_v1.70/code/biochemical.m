function biochem_out = biochemical(biochem_in,Ci_input)
%
% Date: 	21 Sep 2012
% Update:   20 Feb 2013
% Update:      Aug 2013: correction of L171: Ci = Ci*1e6 ./ p .* 1E3;
% Update:   2016-10 - (JAK) major rewrite to accomodate an iterative solution to the Ball-Berry equation
%                   - also allows for g_m to be specified for C3 plants, but only if Ci_input is provided.
% Authors: 	Joe Berry and Christiaan van der Tol, Ari Kornfeld, contributions of others.
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
% the function was tested for Matlab R2013b
%
% Calculates net assimilation rate A, fluorescence F using biochemical model
%
% Input (units are important):
% structure 'biochem_in' with the following elements:
% Knparams   % [], [], []           parameters for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha);
%       [Kno, Kn_alpha, Kn_beta]
%  or, better, as individual fields:
%   Kno                                     Kno - the maximum Kn value ("high light")
%   Kn_alpha, Kn_beta                      alpha, beta: curvature parameters  
%
% Cs        % [ppmV or umol mol]    initial estimate of conc. of CO2 in the
%                                   ...bounary layer of the leaf
% Q         % [umol photons m-2 s-1]net radiation, PAR
% fPAR     % [0-1]                 fraction of incident light that is absorbed by the leaf (default = 1, for compatibility)
% T         % [oC or K]             leaf temperature
% eb        % [hPa = mbar]          intial estimate of the vapour pressure in leaf boundary layer
% O         % [mmol/mol]            concentration of O2 (in the boundary
%                                   ...layer, but no problem to use ambient)
% p         % [hPa]                 air pressure
% Vcmax25 (Vcmo)  % [umol/m2/s]     maximum carboxylation capacity @ 25 degC
% BallBerrySlope (m) % []           Ball-Berry coefficient 'm' for stomatal regulation
% BallBerry0 % []              (OPTIONAL) Ball-Berry intercept term 'b' (if present, an iterative solution is used)
%                                     setting this to zeo disables iteration. Default = 0.01
%
% Type      % ['C3', 'C4']          text parameter, either 'C3' for C3 or any
%                                   ...other text for C4
% tempcor   % [0, 1]               boolean (0 or 1) whether or not
%                                   ...temperature correction to Vcmax has to be applied.
% Tparams  % [],[],[K],[K],[K]     vector of 5 temperature correction parameters, look in spreadsheet of PFTs.
%                                     Only if tempcor=1, otherwise use dummy values
% ...Or replace w/ individual values:
% slti        []              slope of cold temperature decline (C4 only)
% shti        []              slope of high temperature decline in photosynthesis
% Thl         [K]             T below which C4 photosynthesis is <= half that predicted by Q10
% Thh         [K]             T above which photosynthesis is <= half that predicted by Q10
% Trdm        [K]             T at which respiration is <= half that predicted by Q10

% effcon    [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4

% RdPerVcmax25 (Rdparam)  % []     respiration as fraction of Vcmax25
% stressfactor [0-1]               stress factor to reduce Vcmax (for
%                                   example soil moisture, leaf age). Use 1 to "disable" (1 = no stress)
%  OPTIONAL
% Kpep25 (kp)   % [umol/m2/s]         PEPcase activity at 25 deg C (defaults to Vcmax/56
% atheta      % [0-1]                  smoothing parameter for transition between Vc and Ve (light- and carboxylation-limited photosynthesis)
% useTLforC3  % boolean              whether to enable low-temperature attenuation of Vcmax in C3 plants (its always on for C4 plants)
% po0         %  double            Kp,0 (Kp,max) = Fv/Fm (for curve fitting)
% g_m         % mol/m2/s/bar      Mesophyll conductance (default: Infinity, i.e. no effect of g_m)

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

if nargin < 2
    Ci_input = [];
end
%% input
 % environmental
if isfield(biochem_in, 'Cs')
%   assert(all(biochem_in.Cs(:) >=0), 'Negative CO2 (Cs) is not allowed!');
   Cs         = max(0, biochem_in.Cs); % just make sure we don't deal with illegal values
else
    % if Cs is missing, Ci must have been supplied. Forcing Cs = NaN invalidates rcw & gs.
   Cs         = NaN; %biochem_in.Ci;
end
Q             = biochem_in.Q;
assert(all(Q(:) >=0), 'Negative light is not allowed!');

T             = biochem_in.T + 273.15*(biochem_in.T<200); % convert temperatures to K if not already
eb            = biochem_in.eb;
O             = biochem_in.O;
p             = biochem_in.p;

 % physiological
Type          = biochem_in.Type;
if isfield(biochem_in, 'Vcmax25')
    % new field names
    Vcmax25       = biochem_in.Vcmax25;
    BallBerrySlope = 9;
    if isfield(biochem_in, 'BallBerrySlope')  % with g_m and Ci specified, we don't pass BBslope
        BallBerrySlope = biochem_in.BallBerrySlope;
    end
    RdPerVcmax25       = biochem_in.RdPerVcmax25;
else
    % old field names: Vcmo, m, Rdparam
    Vcmax25       = biochem_in.Vcmo;
    BallBerrySlope = biochem_in.m;
    RdPerVcmax25       = biochem_in.Rdparam;
end
BallBerry0 = 0.01; % default value
if isfield(biochem_in, 'BallBerry0')
    BallBerry0 = biochem_in.BallBerry0;
end
if isfield(biochem_in, 'effcon')
    effcon        = biochem_in.effcon;
elseif strcmpi('C3', Type)
    effcon =  1/5;
else
    effcon = 1/6; % C4
end

% Mesophyll conductance: by default we ignore its effect
%  so Cc = Ci - A/gm = Ci
g_m = Inf;
if isfield(biochem_in, 'g_m')
    g_m = biochem_in.g_m * 1e6; % convert from mol to umol
end

% SCOPE provides PAR as APAR, so default (for SCOPE) = 1
%    The curve-fitting GUI may not be providing APAR and should therefore explicitly set fPAR
fPAR = 1;  % fraction of incident light that is absorbed by the leaf
if isfield(biochem_in, 'fPAR')
    fPAR = biochem_in.fPAR;
end

  % physiological options
tempcor       = biochem_in.tempcor;
stressfactor  = biochem_in.stressfactor;
%model_choice  = biochem_in.Fluorescence_model;
if isfield(biochem_in, 'useTLforC3')
    useTLforC3   = biochem_in.useTLforC3;
else
    useTLforC3 = false;
end
    
  % fluoeresence
if isfield(biochem_in, 'Knparams')
    Knparams      = biochem_in.Knparams;
elseif isfield( biochem_in, 'Kn0')
    Knparams = [biochem_in.Kn0, biochem_in.Kn_alpha, biochem_in.Kn_beta];
elseif isfield(biochem_in, 'Fluorescence_model') && biochem_in.Fluorescence_model==0
    % default drought values: 
    Knparams = [5.01, 1.93, 10];
else
    % default general values (cotton dataset)
    Knparams = [2.48, 2.83, 0.114];    
end

if isfield(biochem_in, 'po0')
    po0 = biochem_in.po0;
else
    po0 = [];
end

% physiological temperature parameters: temperature sensitivities of Vcmax, etc 
if isfield(biochem_in, 'Tparams')
   Tparams = biochem_in.Tparams;
    slti        = Tparams(1);
    shti        = Tparams(2);
    Thl         = Tparams(3);
    Thh         = Tparams(4);
    Trdm        = Tparams(5);
else
    slti        = biochem_in.slti;
    shti        = biochem_in.shti;
    Thl         = biochem_in.Thl;
    Thh         = biochem_in.Thh;
    Trdm        = biochem_in.Trdm;
end

%  NOTE: kpep (kp), atheta parameters in next section

%% parameters (at optimum temperature)
Tref        = 25+273.15;        % [K]           absolute temperature at 25 oC

Kc25       = 350;              % [ubar]        kinetic coefficient (Km) for CO2 (Von Caemmerer and Furbank, 1999)
Ko25       = 450;              % [mbar]        kinetic coeeficient (Km) for  O2 (Von Caemmerer and Furbank, 1999)
spfy25     = 2600;      %  specificity (tau in Collatz e.a. 1991)
                        %     This is, in theory, Vcmax/Vomax.*Ko./Kc, but used as a separate parameter

Kpep25       = (Vcmax25/56)*1E6;      % []      (C4) PEPcase rate constant for CO2, used here: Collatz et al: Vcmax25 = 39 umol m-1 s-1; kp = 0.7 mol m-1 s-1.
if isfield(biochem_in,'Kpep25')
    Kpep25 = biochem_in.kpep;
elseif isfield(biochem_in,'kp')
    Kpep25 = biochem_in.kp;
end
if isfield(biochem_in,'atheta') 
    atheta = biochem_in.atheta;
else
    atheta      = 0.8;
end

 % electron transport and fluorescence
Kf          = 0.05;             % []            rate constant for fluorescence
%Kd          = 0.95;             % []           rate constant for thermal deactivation at Fm
Kd          = max(0.8738,  0.0301*(T-273.15)+ 0.0773);
Kp          = 4.0;              % []            rate constant for photochemisty

% note:  rhoa/Mair = L/mol (with the current units) = 24.039 L/mol
%    and  V/n = RT/P ==>  T = 292.95 K @ 1 atm (using R_hPa = 83.144621; 1 atm = 1013.25 hPa)
%   ??!!  These values are used only for rcw, however.
rhoa        = 1.2047;           % [kg m-3]       specific mass of air
Mair        = 28.96;            % [g mol-1]      molecular mass of dry air


%% convert all to bar: CO2 was supplied in ppm, O2 in permil, and pressure in mBar
ppm2bar =  1e-6 .* (p .*1E-3);
Cs          = Cs .* ppm2bar;
O           = (O * 1e-3) .* (p .*1E-3) .* strcmp('C3',Type);    % force O to be zero for C4 vegetation (this is a trick to prevent oxygenase)
Kc25       = Kc25 * 1e-6;
Ko25       = Ko25 * 1e-3;

%% temperature corrections
qt          = 0.1 * (T-Tref) * tempcor;  % tempcorr = 0 or 1: this line dis/enables all Q10 operations
TH          = 1 + tempcor* exp(shti .* (T   -Thh));
%TH          = 1 + tempcor* exp((-220E3+703*T)./(8.314*T));
TL          = 1 + tempcor* exp(slti .* (Thl -T));

QTVc   =  2.1; % Q10 base for Vcmax and Kc
Kc     = Kc25 * exp(log(2.1).*qt);
Ko     = Ko25 * exp(log(1.2).*qt);
kpepcase = Kpep25.* exp(log(1.8).*qt);  % "pseudo first order rate constant for PEP carboxylase WRT pi (Collatz e.a. 1992)


% jak 2014-12-04: Add TL for C3 as well, works much better with our cotton temperature dataset (A-T)
if strcmpi(Type, 'C3') && ~useTLforC3
   Vcmax = Vcmax25 .* exp(log(QTVc).*qt) ./TH * stressfactor;
else
   Vcmax = Vcmax25 .* exp(log(QTVc).*qt) ./(TL.*TH) * stressfactor;
end

% specificity (tau in Collatz e.a. 1991)
spfy        = spfy25 * exp(log(0.75).*qt);

% "Dark" Respiration
Rd          = RdPerVcmax25 * Vcmax25 .* exp(log(1.8).*qt) ./(1+exp(1.3*(T-Trdm)));

%% calculation of potential electron transport rate
if isempty(po0)  % JAK 2015-12: User can specify po0 from measured data
    po0     = Kp./(Kf+Kd+Kp);         % maximum dark photochemistry fraction, i.e. Kn = 0 (Genty et al., 1989)
end
Je          = 0.5*po0 .* Q .* fPAR;          % potential electron transport rate (JAK: add fPAR)

%% calculation of the intersection of enzyme and light limited curves
% this is the original Farquhar model
Gamma_star         = 0.5 .*O ./spfy; %[bar]       compensation point in absence of Rd (i.e. gamma*) [bar]

% Don't bother with...
% Gamma: CO2 compensation point including Rd: solve Ci for 0 = An = Vc - Rd, 
%  assuming Vc dominates at CO2 compensation point according to Farquar 1980. (from Leuning 1990)
%  This gives a realistic value for C4 as well (in which O and Gamma_star = 0)
%Gamma = (Gamma_star .* Vcmax  +   Rd .* MM_consts)  ./ (Vcmax - Rd); % C3
if strcmp(Type, 'C3')
    MM_consts = (Kc .* (1+O./Ko)); % Michaelis-Menten constants
    Vs_C3 = (Vcmax25/2) .* exp(log(1.8).*qt);
    %  minimum Ci (as fraction of Cs) for BallBerry Ci. (If Ci_input is present we need this only as a placeholder for the function call)
    minCi = 0.3;
else
    % C4
    MM_consts = 0; % just for formality, so MM_consts is initialized
    Vs_C3 = 0;     %  the same
    minCi = 0.1;  % C4
end


%% calculation of Ci (internal CO2 concentration)
RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
warnings = [];

fcount = 0; % the number of times we called computeA()
if  ~isempty(Ci_input) 
    Ci = Ci_input; % in units of bar.
    if any(Ci_input > 1)
        % assume Ci_input is in units of ppm. Convert to bar
        Ci = Ci_input .* ppm2bar;
    end
    A =  computeA(Ci);
    
elseif all(BallBerry0 == 0)
    % b = 0: no need to iterate:
    Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);
    A =  computeA(Ci);
    
else
    % compute Ci using iteration (JAK)
    % it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    % Here I use a fully vectorized method based on Brent's method (like fzero) with some optimizations.
    tol = 1e-7;  % 0.1 ppm more-or-less
    % Setting the "corner" argument to Gamma may be useful for low Ci cases, but not very useful for atmospheric CO2, so it's ignored.
    %                     (fn,                           x0, corner, tolerance)
    Ci = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi), Cs, [], tol); % [] in place of Gamma: it didn't make much difference
    %NOTE: A is computed in Ci_next on the final returned Ci. fixedp_brent_ari() guarantees that it was done on the returned values.
    %A =  computeA(Ci);
end

%% Test-function for iteration
%   (note that it assigns A in the function's context.)  
%   As with the next section, this code can be read as if the function body executed at this point. 
%    (if iteration was used). In other words, A is assigned at this point in the file (when iterating).
    function [err, Ci_out] = Ci_next(Ci_in, Cs, RH, minCi)
        % compute the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
        A = computeA(Ci_in);
        A_bar = A .* ppm2bar;
        Ci_out = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi); %[Ci_out, gs]
       
        err = Ci_out - Ci_in; % f(x) - x
    end

%% Compute Assimilation.
%  Note: even though computeA() is written as a separate function, 
%    the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section
    function [A, biochem_out] = computeA(Ci)
        % global: Type, Vcmax, Gamma_star, MM_consts, Vs_C3, effcon, Je, atheta, Rd    %Kc, O, Ko, Vcmax25, qt

        if strcmpi('C3', Type)
           %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.3, Ci_input);
           %effcon      = 0.2;
           % without g_m:
           Vs          = Vs_C3; % = (Vcmax25/2) .* exp(log(1.8).*qt);    % doesn't change on iteration.
           if any(g_m < Inf)
               % with g_m:
               Vc = sel_root( 1./g_m, -(MM_consts + Ci +(Rd + Vcmax)./g_m), Vcmax.*(Ci - Gamma_star + Rd./g_m), -1);
               Ve = sel_root( 1./g_m, -(Ci + 2*Gamma_star +(Rd + Je .* effcon)./g_m), Je .* effcon.*(Ci - Gamma_star + Rd./g_m), -1);
               CO2_per_electron = Ve ./ Je;
           else
               Vc          = Vcmax.*(Ci-Gamma_star)./(MM_consts + Ci);  % MM_consts = (Kc .* (1+O./Ko)) % doesn't change on iteration.
               CO2_per_electron = (Ci-Gamma_star)./(Ci+2*Gamma_star) .* effcon;
               Ve          = Je .* CO2_per_electron;
           end
        else  %C4
            %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.1, Ci_input);
            Vc          = Vcmax;
            Vs          = kpepcase.*Ci;
            %effcon      = 0.17;                    % Berry and Farquhar (1978): 1/0.167 = 6
            CO2_per_electron = effcon; % note: (Ci-Gamma_star)./(Ci+2*Gamma_star) = 1 for C4 (since O = 0); this line avoids 0/0 when Ci = 0
            Ve          = Je .* CO2_per_electron;
        end

        % find the smoothed minimum of Ve, Vc = V, then V, Vs
%         [a1,a2]     = abc(atheta,-(Vc+Ve),Vc.*Ve);
%         % select the min or max  depending on the side of the CO2 compensation point
%         %  note that Vc, Ve < 0 when Ci < Gamma_star (as long as Q > 0; Q = 0 is also ok),
%         %     so the original construction selects the value closest to zero.
%         V           = min(a1,a2).*(Ci>Gamma_star) + max(a1,a2).*(Ci<=Gamma_star);
%         [a1,a2]     = abc(0.98,-(V+Vs),V.*Vs);
%         Ag          = min(a1,a2);
        V           = sel_root(atheta,-(Vc+Ve),Vc.*Ve, sign(-Vc) ); % i.e. sign(Gamma_star - Ci)
        Ag          = sel_root(0.98,-(V+Vs),V.*Vs, -1);
        A           = Ag - Rd;
        
        if nargout > 1
            biochem_out.A = A;
            biochem_out.Ag = Ag;
            biochem_out.Vc = Vc;
            biochem_out.Vs = Vs;
            biochem_out.Ve = Ve;
            biochem_out.CO2_per_electron = CO2_per_electron;
        end
        fcount = fcount + 1; % # of times we called computeA

    end

% (ppm2bar), A_bar
%tic;
%toc
%fprintf('Ball-Berry converged in %d steps (largest_diff = %.4g)\n', counter, largest_diff/ mean(ppm2bar));
%% Compute A, etc.

% note: the following sets a bunch of "global" values in the nested function. Prob better to use [A biochem_out] = ....
%A =  computeA(Ci);  % done above
    % % For debugging:
    % if any(A ~= computeA(Ci) & ~isnan(A))
    %     error('My algorithm didn''t work!');
    % end
%[~, gs1] = BallBerry(Cs, RH, A .* ppm2bar, BallBerrySlope, BallBerry0, minCi, Ci);
gs = 1.6 * A .* ppm2bar ./ (Cs-Ci);

Ja          = Ag ./ CO2_per_electron;        % actual electron transport rate

 % stomatal resistance 
%old: rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3  ./ ppm2bar; %  * 1e6 ./ p .* 1E3;
% if BallBerry0 == 0  %if B-B intercept was specified, then we computed gs "correctly" above and don't need this.
%     rcw(A<=0 & rcw~=0)   = 0.625*1E6;
% end
%rcw         = (1./gs) *rhoa/Mair*1E3  ./ ppm2bar; %  * 1e6 ./ p .* 1E3;
rcw      =  (rhoa./(Mair*1E-3))./gs;

%% fluorescence (Replace this part by Magnani or other model if needed)
ps          = po0.*Ja./Je;               % this is the photochemical yield
nanPs = isnan(ps);
if any(nanPs)
    if numel(po0) == 1
        ps(nanPs) = po0;
    else
        ps(nanPs) = po0(nanPs);  % happens when Q = 0, so ps = po0 (other cases of NaN have been resolved) 
    end
end
ps_rel   = max(0,  1-ps./po0);       % degree of light saturation: 'x' (van der Tol e.a. 2014)

[eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = Fluorescencemodel(ps, ps_rel, Kp,Kf,Kd,Knparams);
Kpa         = ps./fs*Kf;

%% convert back to ppm
Cc = [];
if ~isempty(g_m)
    Cc    = (Ci - A/g_m) ./ ppm2bar;
end
Ci          = Ci  ./ ppm2bar;
%Cs          = Cs  ./ ppm2bar;

%% Collect outputs

biochem_out.A       = A;
biochem_out.Ci      = Ci;
if ~isempty(Cc)
    biochem_out.Cc = Cc;
end
biochem_out.rcw     = rcw;
biochem_out.gs      =  gs;
%  this would be the same if we apply the rcw(A<=0) cutoff:
%biochem_out.gs      = BallBerrySlope.*A.*RH./Cs; % mol/m2/s no intercept term.
biochem_out.RH      =  RH;
biochem_out.warnings = warnings;
biochem_out.fcount = fcount;  % the number of times we called computeA()
%fprintf('fcount = %d\n', fcount);

biochem_out.Vcmax   = Vcmax;
biochem_out.Vc = Vc;  % export the components of A for diagnostic charts
biochem_out.Ve = Ve;
biochem_out.Vs = Vs;
biochem_out.Rd = Rd;

biochem_out.Ja      = Ja;
biochem_out.ps      = ps; % photochemical yield
biochem_out.ps_rel  = ps_rel;   % degree of ETR saturation 'x' (van der Tol e.a. 2014)

 % fluoresence outputs:
% note on Kn: technically, NPQ = (Fm - Fm')/Fm' = Kn/(Kf + Kd);
%     In this model Kf + Kd is close to but NOT equal to 1 @ 25C Kf + Kd = 0.8798
%     vdT 2013 fitted Kn assuming NPQ = Kn, but maybe we shouldn't?
biochem_out.Kd      = Kd;  % K_dark(T)
biochem_out.Kn      = Kn;  % K_n(x);  x = 1 - ps/p00 == 1 - Ja/Je
biochem_out.NPQ     = Kn ./ (Kf + Kd); % why not be honest!
biochem_out.Kf      = Kf;  % Kf = 0.05 (const)
biochem_out.Kp0     = Kp;  % Kp = 4.0 (const): Kp, max
biochem_out.Kp      = Kpa; % Kp,actual
biochem_out.eta     = eta;
biochem_out.qE      = qE;
biochem_out.fs      = fs;  % keep this for compatibility with SCOPE
biochem_out.ft      = fs;  % keep this for the GUI ft is a synonym for what we're calling fs
biochem_out.SIF     = fs .* Q;
biochem_out.fo0     = fo0;
biochem_out.fm0     = fm0;
biochem_out.fo      = fo;
biochem_out.fm      = fm;
biochem_out.Fm_Fo    = fm ./ fo;  % parameters used for curve fitting
biochem_out.Ft_Fo    = fs ./ fo;  % parameters used for curve fitting
biochem_out.qQ      = qQ;
return;

end  % end of function biochemical
 

%% quadratic formula, root of least magnitude
function x = sel_root(a,b,c, dsign)
%  sel_root - select a root based on the fourth arg (dsign = discriminant sign)
%    for the eqn ax^2 + bx + c, 
%    if dsign is:
%       -1, 0: choose the smaller root
%       +1: choose the larger root
%  NOTE: technically, we should check a, but in biochemical, a is always > 0
    if a == 0  % note: this works because 'a' is a scalar parameter!
        x      = -c./b;
    else
        if any(dsign == 0)
            dsign(dsign == 0) = -1; % technically, dsign==0 iff b = c = 0, so this isn't strictly necessary except, possibly for ill-formed cases)
        end
        %disc_root = sqrt(b.^2 - 4.*a.*c); % square root of the discriminant (doesn't need a separate line anymore)
        %  in MATLAB (2013b) assigning the intermediate variable actually slows down the code! (~25%)
        x = (-b + dsign.* sqrt(b.^2 - 4.*a.*c))./(2.*a);
    end
end %of min_root of quadratic formula


%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
% BallBerrySlope, BallBerry0, 
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 6 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A) && nargout > 1
        gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
end

end % function

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
% clean it up:
%gs( Cs == 0 ) = would need to be max gs here;  % eliminate infinities
gs( isnan(Cs) ) = NaN;  % max(NaN, X) = X  (MATLAB 2013b) so fix it here
end


%% Fluorescence model
function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = Fluorescencemodel(ps,x, Kp,Kf,Kd,Knparams)
  % note: x isn't strictly needed as an input parameter but it avoids code-duplication (of po0) and it's inherent risks.
  
    Kno = Knparams(1);
    alpha = Knparams(2);
    beta = Knparams(3);

    % switch model_choice
    %     case 0, % drought 
    %         Kno = 5.01;
    %         alpha = 1.93;
    %         beta = 10;
    %         %Kn          = (6.2473 * x - 0.5944).*x; % empirical fit to Flexas' data
    %         %Kn          = (3.9867 * x - 1.0589).*x;  % empirical fit to Flexas, Daumard, Rascher, Berry data
    %     case 1, healthy (cotton)
    %         Kno = 2.48;
    %         alpha = 2.83;
    %         beta = 0.114;
    %         %p = [4.5531;8.5595;1.8510];
    %         %Kn   = p(1)./(p(3)+exp(-p(2)*(x-.5)));
    % end

    % using exp(-beta) expands the interesting region between 0-1
    %beta = exp(-beta);
    x_alpha = exp(log(x).*alpha); % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
    Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha);

    %Kn          = Kn .* Kd/0.8738;          % temperature correction of Kn similar to that of Kd

    fo0         = Kf./(Kf+Kp+Kd);        % dark-adapted fluorescence yield Fo,0
    fo          = Kf./(Kf+Kp+Kd+Kn);     % light-adapted fluorescence yield in the dark Fo
    fm          = Kf./(Kf   +Kd+Kn);     % light-adapted fluorescence yield Fm
    fm0         = Kf./(Kf   +Kd);        % dark-adapted fluorescence yield Fm
    fs          = fm.*(1-ps);            % steady-state (light-adapted) yield Ft (aka Fs)
    eta         = fs./fo0;
    qQ          = 1-(fs-fo)./(fm-fo);    % photochemical quenching
    qE          = 1-(fm-fo)./(fm0-fo0);  % non-photochemical quenching

    %eta         = eta*(1+5)/5 - 1/5;     % this corrects for 29% PSI contribution in PAM data, but it is quick and dirty correction that needs to be improved in the next 

end


