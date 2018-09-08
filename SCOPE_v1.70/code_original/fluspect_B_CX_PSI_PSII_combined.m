 function leafopt = fluspect_B_CX_PSI_PSII_combined(spectral,leafbio,optipar)
%
% function [leafopt] = fluspect(spectral,leafbio,optipar)
% calculates reflectance and transmittance spectra of a leaf using FLUSPECT, 
% plus four excitation-fluorescence matrices
%
% Authors: Wout Verhoef, Christiaan van der Tol (tol@itc.nl), Joris Timmermans, 
% Date: 2007
% Update from PROSPECT to FLUSPECT: January 2011 (CvdT)
%
%      Nov 2012 (CvdT) Output EF-matrices separately for PSI and PSII
%   31 Jan 2013 (WV)   Adapt to SCOPE v_1.40, using structures for I/O
%   30 May 2013 (WV)   Repair bug in s for non-conservative scattering
%   24 Nov 2013 (WV)   Simplified doubling routine
%   25 Nov 2013 (WV)   Restored piece of code that takes final refl and
%                      tran outputs as a basis for the doubling routine
%   03 Dec 2013 (WV)   Major upgrade. Border interfaces are removed before 
%                      the fluorescence calculation and later added again
%   23 Dec 2013 (WV)   Correct a problem with N = 1 when calculating k 
%                      and s; a test on a = Inf was included
%   01 Apr 2014 (WV)   Add carotenoid concentration (Cca and Kca)
%   19 Jan 2015 (WV)   First beta version for simulation of PRI effect
%   17 Mar 2017 (CT)   Added Anthocyanins according to Prospect-D
%
% usage:
% [leafopt] = fluspect_b(spectral,leafbio,optipar)
% 
% inputs:
% Cab         = leafbio.Cab;
% Cca         = leafbio.Cca;
% V2Z         = leafbio.V2Z;  % Violaxanthin - Zeaxanthin transition status
%                               [0-1]
% Cw          = leafbio.Cw;
% Cdm         = leafbio.Cdm;
% Cs          = leafbio.Cs;
% Cant 	      = leafbio.Cant;
% N           = leafbio.N; 
% fqe         = leafbio.fqe;
% 
% nr          = optipar.nr;
% Kdm         = optipar.Kdm;
% Kab         = optipar.Kab;
% Kca         = optipar.Kca;
% KcaV        = optipar.KcaV;
% KcaZ        = optipar.KcaZ;
% Kw          = optipar.Kw;
% Ks          = optipar.Ks;
% phi         = optipar.phi;
% outputs:
% refl          reflectance
% tran          transmittance
% Mb            backward scattering fluorescence matrix, I for PSI and II for PSII
% Mf            forward scattering fluorescence matrix,  I for PSI and II for PSII

%% parameters
% fixed parameters for the fluorescence module
ndub        = 15;           % number of doublings applied

% Fluspect parameters
Cab         = leafbio.Cab;
Cca         = leafbio.Cca;
V2Z         = leafbio.V2Z;
Cw          = leafbio.Cw;
Cdm         = leafbio.Cdm;
Cs          = leafbio.Cs;
Cant 	    = leafbio.Cant;
N           = leafbio.N;
fqe         = leafbio.fqe;

nr          = optipar.nr;
Kdm         = optipar.Kdm;
Kab         = optipar.Kab;

if V2Z == -999 
    % Use old Kca spectrum if this is given as input
    Kca     = optipar.Kca;
else
    % Otherwise make linear combination based on V2Z
    % For V2Z going from 0 to 1 we go from Viola to Zea
    Kca     = (1-V2Z) * optipar.KcaV + V2Z * optipar.KcaZ;    
end

Kw          = optipar.Kw;
Ks          = optipar.Ks;
Kant        = optipar.Kant;
phi         = optipar.phi;

%% PROSPECT calculations
Kall        = (Cab*Kab + Cca*Kca + Cdm*Kdm + Cw*Kw  + Cs*Ks + Cant*Kant)/N;   % Compact leaf layer

j           = find(Kall>0);               % Non-conservative scattering (normal case)
t1          = (1-Kall).*exp(-Kall);
t2          = Kall.^2.*expint(Kall);
tau         = ones(size(t1));
tau(j)      = t1(j)+t2(j);
kChlrel     = zeros(size(t1));
kChlrel(j)  = Cab*Kab(j)./(Kall(j)*N);

talf        = calctav(59,nr);
ralf        = 1-talf;
t12         = calctav(90,nr);
r12         = 1-t12;
t21         = t12./(nr.^2);
r21         = 1-t21;

% top surface side
denom       = 1-r21.*r21.*tau.^2;
Ta          = talf.*tau.*t21./denom;
Ra          = ralf+r21.*tau.*Ta;

% bottom surface side
t           = t12.*tau.*t21./denom;
r           = r12+r21.*tau.*t;

% Stokes equations to compute properties of next N-1 layers (N real)
% Normal case

D           = sqrt((1+r+t).*(1+r-t).*(1-r+t).*(1-r-t));
rq          = r.^2;
tq          = t.^2;
a           = (1+rq-tq+D)./(2*r);
b           = (1-rq+tq+D)./(2*t);

bNm1        = b.^(N-1);                  %
bN2         = bNm1.^2;
a2          = a.^2;
denom       = a2.*bN2-1;
Rsub        = a.*(bN2-1)./denom;
Tsub        = bNm1.*(a2-1)./denom;

%			Case of zero absorption
j           = find(r+t >= 1);
Tsub(j)     = t(j)./(t(j)+(1-t(j))*(N-1));
Rsub(j)	    = 1-Tsub(j);

% Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
denom       = 1-Rsub.*r;
tran        = Ta.*Tsub./denom;
refl        = Ra+Ta.*Rsub.*t./denom;

leafopt.refl = refl;
leafopt.tran = tran;
leafopt.kChlrel = kChlrel;

% From here a new path is taken: The doubling method used to calculate
% fluoresence is now only applied to the part of the leaf where absorption
% takes place, that is, the part exclusive of the leaf-air interfaces. The
% reflectance (rho) and transmittance (tau) of this part of the leaf are
% now determined by "subtracting" the interfaces

Rb  = (refl-ralf)./(talf.*t21+(refl-ralf).*r21);  % Remove the top interface
Z   = tran.*(1-Rb.*r21)./(talf.*t21);             % Derive Z from the transmittance

rho = (Rb-r21.*Z.^2)./(1-(r21.*Z).^2);    % Reflectance and transmittance 
tau = (1-Rb.*r21)./(1-(r21.*Z).^2).*Z;    % of the leaf mesophyll layer
t   =   tau;
r   =   max(rho,0);                       % Avoid negative r

% Derive Kubelka-Munk s and k

I_rt     =   (r+t)<1;
D(I_rt)  =   sqrt((1 + r(I_rt) + t(I_rt)) .* ...
                  (1 + r(I_rt) - t(I_rt)) .* ...
                  (1 - r(I_rt) + t(I_rt)) .* ...
                  (1 - r(I_rt) - t(I_rt)));
a(I_rt)  =   (1 + r(I_rt).^2 - t(I_rt).^2 + D(I_rt)) ./ (2*r(I_rt));
b(I_rt)  =   (1 - r(I_rt).^2 + t(I_rt).^2 + D(I_rt)) ./ (2*t(I_rt));    
a(~I_rt) =   1;
b(~I_rt) =   1;

s        =   r./t;
I_a      =   (a>1 & a~=Inf);
s(I_a)   =   2.*a(I_a) ./ (a(I_a).^2 - 1) .* log(b(I_a));

k        =    log(b);
k(I_a)   =   (a(I_a)-1) ./ (a(I_a)+1) .* log(b(I_a));
kChl     =   kChlrel .* k;

%% Fluorescence of the leaf mesophyll layer
% Fluorescence part is skipped for fqe = 0

if fqe > 0

    wle         = spectral.wlE';    % excitation wavelengths, transpose to column
    wlf         = spectral.wlF';    % fluorescence wavelengths, transpose to column
    wlp         = spectral.wlP;     % PROSPECT wavelengths, kept as a row vector

    minwle      = min(wle);
    maxwle      = max(wle);
    minwlf      = min(wlf);
    maxwlf      = max(wlf);

    % indices of wle and wlf within wlp

    Iwle        = find(wlp>=minwle & wlp<=maxwle);
    Iwlf        = find(wlp>=minwlf & wlp<=maxwlf);

    eps         = 2^(-ndub);

    % initialisations
    te          = 1-(k(Iwle)+s(Iwle)) * eps;   
    tf          = 1-(k(Iwlf)+s(Iwlf)) * eps;  
    re          = s(Iwle) * eps;
    rf          = s(Iwlf) * eps;

    sigmoid     = 1./(1+exp(-wlf/10)*exp(wle'/10));  % matrix computed as an outproduct
    
    [Mf, Mb] = deal(fqe(1) * ((.5*phi(Iwlf))*eps) * kChl(Iwle)'.*sigmoid);
    
    Ih          = ones(1,length(te));     % row of ones
    Iv          = ones(length(tf),1);     % column of ones

    % Doubling routine
    
    for i = 1:ndub
        
        xe = te./(1-re.*re);  ten = te.*xe;  ren = re.*(1+ten);  
        xf = tf./(1-rf.*rf);  tfn = tf.*xf;  rfn = rf.*(1+tfn);
              
        A11  = xf*Ih + Iv*xe';           A12 = (xf*xe').*(rf*Ih + Iv*re');
        A21  = 1+(xf*xe').*(1+rf*re');   A22 = (xf.*rf)*Ih+Iv*(xe.*re)';
        
        Mfn   = Mf  .* A11 + Mb  .* A12;
        Mbn   = Mb  .* A21 + Mf  .* A22;
        
        te   = ten;  re  = ren;   tf   = tfn;   rf   = rfn;
        Mf  = Mfn; Mb = Mbn;
    end
    
    % Here we add the leaf-air interfaces again for obtaining the final 
    % leaf level fluorescences.
    
    g = Mb; f = Mf;
    
    Rb = rho + tau.^2.*r21./(1-rho.*r21);
        
    Xe = Iv * (talf(Iwle)./(1-r21(Iwle).*Rb(Iwle)))';
    Xf = t21(Iwlf)./(1-r21(Iwlf).*Rb(Iwlf)) * Ih;
    Ye = Iv * (tau(Iwle).*r21(Iwle)./(1-rho(Iwle).*r21(Iwle)))';
    Yf = tau(Iwlf).*r21(Iwlf)./(1-rho(Iwlf).*r21(Iwlf)) * Ih;
    
    A = Xe .* (1 + Ye.*Yf) .* Xf;
    B = Xe .* (Ye + Yf) .* Xf;
    
    gn = A .* g + B .* f;
    fn = A .* f + B .* g;
    
    leafopt.Mb  = gn;
    leafopt.Mf  = fn;
    
    [leafopt.MbI_rc,leafopt.MfI_rc] = deal(0.5*fqe(1) * ((phi(Iwlf)))  * kChlrel(Iwle)'.*sigmoid);
    
end    

return;

function tav = calctav(alfa,nr)

    rd          = pi/180;
    n2          = nr.^2;
    np          = n2+1;
    nm          = n2-1;
    a           = (nr+1).*(nr+1)/2;
    k           = -(n2-1).*(n2-1)/4;
    sa          = sin(alfa.*rd);

    b1          = (alfa~=90)*sqrt((sa.^2-np/2).*(sa.^2-np/2)+k);
    b2          = sa.^2-np/2;
    b           = b1-b2;
    b3          = b.^3;
    a3          = a.^3;
    ts          = (k.^2./(6*b3)+k./b-b/2)-(k.^2./(6*a3)+k./a-a/2);

    tp1         = -2*n2.*(b-a)./(np.^2);
    tp2         = -2*n2.*np.*log(b./a)./(nm.^2);
    tp3         = n2.*(1./b-1./a)/2;
    tp4         = 16*n2.^2.*(n2.^2+1).*log((2*np.*b-nm.^2)./(2*np.*a-nm.^2))./(np.^3.*nm.^2);
    tp5         = 16*n2.^3.*(1./(2*np.*b-nm.^2)-1./(2*np.*a-nm.^2))./(np.^3);
    tp          = tp1+tp2+tp3+tp4+tp5;
    tav         = (ts+tp)./(2*sa.^2);

return;