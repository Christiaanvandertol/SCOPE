function biochem_out = biochemical(leafbio,meteo,options,constants,fV)
%
% Date: 	21 Sep 2012
% Update:   20 Feb 2013
% Update:      Aug 2013: correction of L171: Ci = Ci*1e6 ./ p .* 1E3;
% Update:   2016-10 - (JAK) major rewrite to accomodate an iterative solution to the Ball-Berry equation
%                   - also allows for g_m to be specified for C3 plants, but only if Ci_input is provided.
% Update:   25 Feb 2021: Temperature reponse functions by Dutta et al.
%                       implemented
% Authors: 	Joe Berry and Christiaan van der Tol, Ari Kornfeld, contributions of others.
% Sources:
%           Farquhar et al. 1980, Collatz et al (1991, 1992), and:
%
% Dutta, D., Schimel, D. S., Sun, Y., Tol, C. V. D., & Frankenberg, C. (2019). 
% Optimal inverse estimation of ecosystem parameters from observations of carbon and energy fluxes. 
% Biogeosciences, 16(1), 77-103.
%
% Van der Tol, C., Berry, J. A., Campbell, P. K. E., & Rascher, U. (2014). 
% Models of fluorescence and photosynthesis for interpreting measurements 
% of solar?induced chlorophyll fluorescence. 
% Journal of Geophysical Research: Biogeosciences, 119(12), 2312-2327.
%
% Bonan, G. B., Lawrence, P. J., Oleson, K. W., Levis, S., Jung, M., 
% Reichstein, M., ... & Swenson, S. C. (2011). 
% Improving canopy processes in the Community Land Model version 4 (CLM4) 
% using global flux fields empirically inferred from FLUXNET data. 
% Journal of Geophysical Research: Biogeosciences, 116(G2).
%
%
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1)
%    - photosynthesis of a leaf or needle (umol m-2 s-1)
%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
%
% Usage:
% biochem_out = biochemical(leafbio,meteo,options,constants,fV)
% the function was tested for Matlab R2017a
%
% Calculates net assimilation rate A, fluorescence F using biochemical model
%
% Input (units are important):
% structure 'leafbio' with the following elements:
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
%
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

%% physiological options
tempcor       = options.apply_T_corr;

%% input

% atmospheric physics
rhoa            = constants.rhoa;           % [kg m-3]       specific mass of air
Mair            = constants.Mair;           % [g mol-1]      molecular mass of dry air
R               = constants.R;              % [J mol-1K-1]   Molar gas constant
Q               = meteo.Q;                  % [umol m-2 s-1] absorbed PAR flux
Cs              = meteo.Cs;
T               = meteo.T + 273.15*(meteo.T<200); % convert temperatures to K if not already
eb              = meteo.eb;
O               = meteo.Oa;
p               = meteo.p;

%biochemical
Type            = leafbio.Type;
Vcmax25         = fV.*leafbio.Vcmax25;
BallBerrySlope  = leafbio.BallBerrySlope;
RdPerVcmax25    = leafbio.RdPerVcmax25;
BallBerry0      = leafbio.BallBerry0;
Tref            = 25+273.15;        % [K]           absolute temperature at 25 oC

%Jmax25      = 1.97*Vcmax25;
%TPU25           = 0.06*Jmax25;              % Triose Phosphate utilization rate

% All of the Following Values are Adopted from Bonan et. al 2011 paper and
% pertaings to C3 photosythesis

Kc25            = 405;                   % [umol mol-1]
Ko25            = 279;                   % [mmol mol-1]
spfy25          = 2444;                  % specificity (Computed from Bernacchhi et al 2001 paper)

% convert all to bar: CO2 was supplied in ppm, O2 in permil, and pressure in mBar
ppm2bar     =  1e-6 .* (p .*1E-3);
Cs          = Cs .* ppm2bar;
O           = (O * 1e-3) .* (p .*1E-3) .* strcmp('C3',Type);    % force O to be zero for C4 vegetation (this is a trick to prevent oxygenase)
Kc25        = Kc25 * 1e-6;
Ko25        = Ko25 * 1e-3;
Gamma_star25    = 0.5 .*O./spfy25;      % [ppm] compensation point in absence of Rd
Rd25            = RdPerVcmax25 * Vcmax25;
if strcmpi('C3', Type)
    effcon =  1/5;
else
    effcon = 1/6; % C4
end
atheta      = 0.8;

% Mesophyll conductance: by default we ignore its effect
%  so Cc = Ci - A/gm = Ci
g_m = Inf;
if isfield(leafbio, 'g_m')
    g_m = leafbio.g_m * 1e6; % convert from mol to umol
end
stressfactor  = leafbio.stressfactor;

% fluorescence
Knparams    = [leafbio.Kn0, leafbio.Knalpha, leafbio.Knbeta];
Kf          = 0.05;             % []            rate constant for fluorescence
Kd          = max(0.8738,  0.0301*(T-273.15)+ 0.0773);
Kp          = 4.0;              % []            rate constant for photochemisty

%% temperature corrections

[f.Vcmax,f.Rd,f.TPU,f.Kc,f.Ko,f.Gamma_star] = deal(1);
if tempcor
    if strcmpi('C4', Type)
        % RdPerVcmax25 = 0.025;  % Rd25 for C4 is different than C3
        %   Rd25 = RdPerVcmax25 * Vcmax25;
        % Constant parameters for temperature correction of Vcmax
        Q10 = leafbio.TDP.Q10;                           % Unit is  []
        s1  = leafbio.TDP.s1;                            % Unit is [K]
        s2  = leafbio.TDP.s2;                            % Unit is [K^-1]
        s3  = leafbio.TDP.s3;                            % Unit is [K]
        s4  = leafbio.TDP.s4;                            % Unit is [K^-1]
        
        % Constant parameters for temperature correction of Rd
        s5  = leafbio.TDP.s5;                            % Unit is [K]
        s6  = leafbio.TDP.s6;                            % Unit is [K^-1]
        
        fHTv = 1 + exp(s1.*(T - s2));
        fLTv = 1 + exp(s3.*(s4 - T));
        Vcmax = (Vcmax25 .* Q10.^(0.1.*(T-Tref)))./(fHTv .* fLTv); % Temp Corrected Vcmax
        
        % Temperature correction of Rd
        
        fHTv = 1 + exp(s5.*(T - s6));
        Rd = (Rd25 .* Q10.^(0.1.*(T-Tref)))./fHTv; % Temp Corrected Rd
        % Temperature correction of Ke
        Ke25 = 20000 .* Vcmax25 ;               % Unit is  []
        Ke = (Ke25 .* Q10.^(0.1.*(T-Tref)));    % Temp Corrected Ke
        
    else
        % temperature correction of Vcmax
        deltaHa     = leafbio.TDP.delHaV;                % Unit is  [J K^-1]
        deltaS      = leafbio.TDP.delSV;                 % unit is [J mol^-1 K^-1]
        deltaHd     = leafbio.TDP.delHdV;                % unit is [J mol^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
        f.Vcmax     = fTv .* fHTv;
        
%         % temperature correction for TPU
%         deltaHa     = leafbio.TDP.delHaP;                % Unit is  [J K^-1]
%         deltaS      = leafbio.TDP.delSP;                 % unit is [J mol^-1 K^-1]
%         deltaHd     = leafbio.TDP.delHdP;                % unit is [J mol^-1]
%         fTv         = temperature_functionC3(Tref,R,T,deltaHa);
%         fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
%         f.TPU       = fTv .* fHTv;
        
        % temperature correction for Rd
        deltaHa     = leafbio.TDP.delHaR;               % Unit is  [J K^-1]
        deltaS      = leafbio.TDP.delSR;                % unit is [J mol^-1 K^-1]
        deltaHd     = leafbio.TDP.delHdR;               % unit is [J mol^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        fHTv        = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd);
        f.Rd        = fTv .* fHTv;
        
        % temperature correction for Kc
        deltaHa     = leafbio.TDP.delHaKc;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Kc        = fTv;
        
        % temperature correction for Ko
        deltaHa     = leafbio.TDP.delHaKo;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Ko        = fTv;
        
        % temperature correction for Gamma_star
        deltaHa     = leafbio.TDP.delHaT;               % Unit is  [J K^-1]
        fTv         = temperature_functionC3(Tref,R,T,deltaHa);
        f.Gamma_star = fTv;
        
        Ke          = 1; % dummy value (only needed for C4)
    end
else
    Ke          = 1; % dummy value (only needed for C4)
end

if strcmp('C3', Type)
    Vcmax       = Vcmax25.*f.Vcmax.*stressfactor;
    Rd          = Rd25    .*f.Rd.* stressfactor;
    %TPU         = TPU25   .* f.TPU;
    Kc          = Kc25    .* f.Kc;
    Ko          = Ko25    .* f.Ko;
end
Gamma_star   = Gamma_star25 .* f.Gamma_star;

%% calculation of potential electron transport rate
po0         = Kp./(Kf+Kd+Kp);         % maximum dark photochemistry fraction, i.e. Kn = 0 (Genty et al., 1989)
Je          = 0.5*po0 .* Q;          % potential electron transport rate (JAK: add fPAR);
%Gamma_star  = 0.5 .*O ./spfy; %[bar]       compensation point in absence of Rd (i.e. gamma*) [bar]

if strcmp(Type, 'C3')
    MM_consts = (Kc .* (1+O./Ko)); % Michaelis-Menten constants
    Vs_C3 = (Vcmax/2);
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
computeA()  % clears persistent fcount
computeA_fun = @(x) computeA(x, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, Ke);

if all(BallBerry0 == 0)
    % b = 0: no need to iterate:
    Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);
    %     A =  computeA_fun(Ci);   
else
    % compute Ci using iteration (JAK)
    % it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    % Here I use a fully vectorized method based on Brent's method (like fzero) with some optimizations.
    tol = 1e-7;  % 0.1 ppm more-or-less
    % Setting the "corner" argument to Gamma may be useful for low Ci cases, but not very useful for atmospheric CO2, so it's ignored.
    %                     (fn,                           x0, corner, tolerance)
    [Ci] = fixedp_brent_ari(@(x) Ci_next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0, computeA_fun, ppm2bar), Cs, [], tol); % [] in place of Gamma: it didn't make much difference
    %NOTE: A is computed in Ci_next on the final returned Ci. fixedp_brent_ari() guarantees that it was done on the returned values.
    %     A =  computeA_fun(Ci);
end

[A, biochem_out]    = computeA_fun(Ci);
Ag                  = biochem_out.Ag;
CO2_per_electron    = biochem_out.CO2_per_electron;

%% Compute A, electron transport rate, and stomatal resistance
gs                  = max(0, 1.6 * A .* ppm2bar ./ (Cs-Ci));     % stomatal conductance
Ja                  = Ag ./ CO2_per_electron;   % actual electron transport rate
rcw                 = (rhoa./(Mair*1E-3))./gs;    % stomatal resistance

%% fluorescence
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
biochem_out.RH      =  RH;
biochem_out.Vcmax   = Vcmax;
biochem_out.Rd      = Rd;
biochem_out.Ja      = Ja;
biochem_out.ps      = ps; % photochemical yield
biochem_out.ps_rel  = ps_rel;   % degree of ETR saturation 'x' (van der Tol e.a. 2014)
biochem_out.Kd      = Kd;  % K_dark(T)
biochem_out.Kn      = Kn;  % K_n(x);  x = 1 - ps/p00 == 1 - Ja/Je
biochem_out.NPQ     = Kn ./ (Kf + Kd); % 
biochem_out.Kf      = Kf;  % Kf = 0.05 (const)
biochem_out.Kp0     = Kp;  % Kp = 4.0 (const): Kp, max
biochem_out.Kp      = Kpa; % Kp,actual
biochem_out.eta     = eta;
biochem_out.qE      = qE;
biochem_out.fs      = fs;  % keep this for compatibility with SCOPE
biochem_out.ft      = fs;  % keep this 
biochem_out.SIF     = fs .* Q;
biochem_out.fo0     = fo0;
biochem_out.fm0     = fm0;
biochem_out.fo      = fo;
biochem_out.fm      = fm;
biochem_out.Fm_Fo   = fm ./ fo;  % parameters used for curve fitting
biochem_out.Ft_Fo   = fs ./ fo;  % parameters used for curve fitting
biochem_out.qQ      = qQ;
biochem_out.Phi_N   = Kn./(Kn +Kp+Kf+Kd);
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

x_alpha = exp(log(x).*alpha); % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha);

fo0         = Kf./(Kf+Kp+Kd);        % dark-adapted fluorescence yield Fo,0
fo          = Kf./(Kf+Kp+Kd+Kn);     % light-adapted fluorescence yield in the dark Fo
fm          = Kf./(Kf   +Kd+Kn);     % light-adapted fluorescence yield Fm
fm0         = Kf./(Kf   +Kd);        % dark-adapted fluorescence yield Fm
fs          = fm.*(1-ps);            % steady-state (light-adapted) yield Ft (aka Fs)
eta         = fs./fo0;
qQ          = 1-(fs-fo)./(fm-fo);    % photochemical quenching
qE          = 1-(fm-fo)./(fm0-fo0);  % non-photochemical quenching

end

%% Test-function for iteration
%   (note that it assigns A in the function's context.)
%   As with the next section, this code can be read as if the function body executed at this point.
%    (if iteration was used). In other words, A is assigned at this point in the file (when iterating).
function [err, Ci_out] = Ci_next(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar)
% compute the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
A = A_fun(Ci_in);
A_bar = A .* ppm2bar;
Ci_out = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi); %[Ci_out, gs]

err = Ci_out - Ci_in; % f(x) - x
end

%% Compute Assimilation.
%  Note: even though computeA() is written as a separate function,
%    the code is, in fact, executed exactly this point in the file (i.e. between the previous if clause and the next section
function [A, biochem_out] = computeA(Ci, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, kpepcase)
% global: Type, Vcmax, Gamma_star, MM_consts, Vs_C3, effcon, Je, atheta, Rd    %Kc, O, Ko, Vcmax25, qt
persistent fcount
if nargin == 0
    fcount = 0;
    return
end
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
fcount = fcount + 1; % # of times we called computeA

if nargout > 1
    biochem_out.A = A;
    biochem_out.Ag = Ag;
    biochem_out.Vc = Vc;
    biochem_out.Vs = Vs;
    biochem_out.Ve = Ve;
    biochem_out.CO2_per_electron = CO2_per_electron;
    biochem_out.fcount = fcount;
end

end

%% Temperature Correction Functions
% The following two functions pertains to C3 photosynthesis
function [fTv] = temperature_functionC3(Tref,R,T,deltaHa)
% Temperature function
tempfunc1 = (1 - Tref./T);
fTv = exp(deltaHa/(Tref*R).*tempfunc1);
end

function [fHTv] = high_temp_inhibtionC3(Tref,R,T,deltaS,deltaHd)
% High Temperature Inhibition Function
hightempfunc_num = (1+exp((Tref*deltaS-deltaHd)/(Tref*R)));
hightempfunc_deno = (1+exp((deltaS.*T - deltaHd)./(R.*T)));
fHTv = hightempfunc_num ./ hightempfunc_deno;
end





