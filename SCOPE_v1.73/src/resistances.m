function [resist_out] = resistances(resist_in)
%
%   function resistances calculates aerodynamic and boundary resistances
%   for soil and vegetation
%
%   Date:       01 Feb 2008
%   Authors:    Anne Verhoef            (a.verhoef@reading.ac.uk)
%               Christiaan van der Tol  (tol@itc.nl)
%               Joris Timmermans        (j_timmermans@itc.nl)
%   Source:     Wallace and Verhoef (2000) 'Modelling interactions in
%               mixed-plant communities: light, water and carbon dioxide', in: Bruce
%               Marshall, Jeremy A. Roberts (ed), 'Leaf Development and Canopy Growth',
%               Sheffield Academic Press, UK. ISBN 0849397693
%               
%               ustar:  Tennekes, H. (1973) 'The logaritmic wind profile', J.
%               Atmospheric Science, 30, 234-238
%               Psih:   Paulson, C.A. (1970), The mathematical
%               representation of wind speed and temperature in the
%               unstable atmospheric surface layer. J. Applied Meteorol. 9,
%               857-861
%       
% Note: Equation numbers refer to equation numbers in Wallace and Verhoef (2000)
% 
% Usage:
%   [resist_out] = resistances(resist_in)
%
% The input and output are structures. These structures are further
% specified in a readme file.
%
% Input: 
%   resist_in   aerodynamic resistance parameters and wind speed
%
% The strucutre resist_in contains the following elements:
% u         =   windspeed
% L         =   stability
% LAI       =   Leaf Area Index

% rbs       =   Boundary Resistance of soil                         [s m-1]
% rss       =   Surface resistance of soil for vapour transport     [s m-1]
% rwc       =   Within canopy Aerodynamic Resistance canopy         [s m-1]

% z0m       =   Roughness lenght for momentum for the vegetation    [m]
% d         =   Displacement height (Zero plane)                    [m]
% z         =   Measurement height                                  [m]
% h         =   Vegetation height                                   [m]

%
% Output:
%   resist_out  aeorodynamic resistances
%
% The strucutre resist_out contains the following elements:
% ustar     =   Friction velocity                                   [m s-1]
% raa       =   Aerodynamic resistance above the canopy             [s m-1]                     
% rawc      =   Total resistance within the canopy (canopy)         [s m-1]
% raws      =   Total resistance within the canopy (soil)           [s m-1]

% rai       =   Aerodynamic resistance in inertial sublayer         [s m-1]
% rar       =   Aerodynamic resistance in roughness sublayer        [s m-1]
% rac       =   Aerodynamic resistance in canopy layer (above z0+d) [s m-1]

% rbc       =   Boundary layer resistance (canopy)                  [s m-1]
% rwc       =   Aerodynamic Resistance within canopy(canopy)(Update)[s m-1]

% rbs       =   Boundary layer resistance (soil) (Update)           [s m-1]
% rws       =   Aerodynamic resistance within canopy(soil)          [s m-1] 

% rss       =   Surface resistance vapour transport(soil)(Update)   [s m-1]

% uz0       =   windspeed at z0                                     [m s-1]
% Kh        =   Diffusivity for heat                                [m2s-1]

%% parameters
global constants
kappa = constants. kappa;

Cd        =  resist_in.Cd;

u         =  resist_in.u;
L         =  resist_in.L;
LAI       =  resist_in.LAI;

rbs       =  resist_in.rbs;
%rss       =  resist_in.rss;
rwc       =  resist_in.rwc;

z0m       =  resist_in.zo;
d         =  resist_in.d;
z         =  resist_in.z;
h         =  resist_in.hc;
w         =  resist_in.w;

% derived parameters
%zr: top of roughness sublayer, bottom of intertial sublayer
zr			= 2.5*h;                   %                            [m]			
%n: dimensionless wind extinction coefficient                       W&V Eq 33
n			= Cd*LAI/(2*kappa^2);      %                            [] 

%% stability correction for non-neutral conditions
%neu		= find(L >= -.001 & L <= .001);
unst        = find(L <  -4);
st          = find(L >  4E3);

% stability correction functions, friction velocity and Kh=Km=Kv
pm_z    	= psim(z -d,L,unst,st);
ph_z    	= psih(z -d,L,unst,st);
pm_h        = psim(h -d,L,unst,st);
%ph_h       = psih(h -d,L,unst,st);
ph_zr       = psih(zr-d,L,unst,st).*(z>=zr) + ph_z.*(z<zr);
phs_zr      = phstar(zr,zr,d,L,st,unst);
phs_h		= phstar(h ,zr,d,L,st,unst);

ustar   	= max(.001,kappa*u./(log((z-d)/z0m) - pm_z));%          W&V Eq 30

Kh                  = kappa*ustar*(zr-d);                  %                W&V Eq 35
resist_out.Kh(unst)	= kappa*ustar(unst)*(zr-d).*(1-16*(h-d)./L(unst)).^.5;% W&V Eq 35
resist_out.Kh(st)   = kappa*ustar(st)  *(zr-d).*(1+ 5*(h-d)./L(st)  ).^-1;% W&V Eq 35

%% wind speed at height h and z0m
uh			= max(ustar/kappa .* (log((h-d)/z0m) - pm_h     ),.01);
uz0 		= uh*exp(n*((z0m+d)/h-1));                      %       W&V Eq 32

%% resistances

resist_out.uz0 = uz0; 
resist_out.ustar = ustar;
rai = (z>zr).*(1./(kappa*ustar).*(log((z-d) /(zr-d))  - ph_z   + ph_zr));% W&V Eq 41 
rar = 1./(kappa*ustar).*((zr-h)/(zr-d)) 	 - phs_zr + phs_h;% W&V Eq 39
rac = h*sinh(n)./(n*Kh)*(log((exp(n)-1)/(exp(n)+1)) - log((exp(n*(z0m+ d )/h)-1)/(exp(n*(z0m +d )/h)+1))); % W&V Eq 42
rws = h*sinh(n)./(n*Kh)*(log((exp(n*(z0m+d)/h)-1)/(exp(n*(z0m+d)/h)+1)) - log((exp(n*(.01    )/h)-1)/(exp(n*(.01    )/h)+1))); % W&V Eq 43
rbc = 70/LAI * sqrt(w./uz0);						%		W&V Eq 31, but slightly different

resist_out.rai = rai;
resist_out.rar = rar;
resist_out.rac = rac;
resist_out.rws = rws;
resist_out.rbc = rbc;

raa  = rai + rar + rac;
rawc = rwc + rbc;
raws = rws + rbs;

resist_out.raa  = raa;          % aerodynamic resistance above the canopy           W&V Figure 8.6
resist_out.rawc	= rawc;			% aerodynamic resistance within the canopy (canopy)
resist_out.raws	= raws;			% aerodynamic resistance within the canopy (soil)

resist_out.raa  = min(4E2,raa);         % to prevent unrealistically high resistances
resist_out.rawc = min(4E2,rawc);        % to prevent unrealistically high resistances
resist_out.raws = min(4E2,raws);        % to prevent unrealistically high resistances

return


%% subfunction pm for stability correction (eg. Paulson, 1970)
function pm = psim(z,L,unst,st)
x       	= (1-16*z./L(unst)).^(1/4);
pm      	= zeros(size(L));
pm(unst)	= 2*log((1+x)/2)+log((1+x.^2)/2) -2*atan(x)+pi/2;   %   unstable
pm(st)      = -5*z./L(st);                                      %   stable
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function ph = psih(z,L,unst,st)
x       	= (1-16*z./L(unst)).^(1/4);
ph      	= zeros(size(L));
ph(unst)	= 2*log((1+x.^2)/2);                                %   unstable
ph(st)      = -5*z./L(st);                                      %   stable
return

%% subfunction ph for stability correction (eg. Paulson, 1970)
function phs = phstar(z,zR,d,L,st,unst)
x			= (1-16*z./L(unst)).^0.25;
phs 		= zeros(size(L));
phs(unst)   = (z-d)/(zR-d)*(x.^2-1)./(x.^2+1);
phs(st)     = -5*z./L(st);
return