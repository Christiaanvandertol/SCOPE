function [zom,d] = zo_and_d (soil,canopy)

% function zom_and_d calculates roughness length for momentum and zero
% plane displacement from vegetation height and LAI
%
% Date:     17 November 2008
%           17 April 2013 (structures)
%
% Author:   A. Verhoef
%           implemented into Matlab by C. van der Tol (c.vandertol@utwente.nl)
%
% Source:   Verhoef, McNaughton & Jacobs (1997), HESS 1, 81-91
% 
% usage: 
%       zo_and_d (soil,canopy)
%
% canopy fields used as inpuyt:
%   LAI         one sided leaf area index
%   hc           vegetation height (m)
%
% soil fields used:
%   Cd          Averaged drag coefficient for the vegetation              
%   CR          Drag coefficient for isolated tree
%   CSSOIL      Drag coefficient for soil
%   CD1         Fitting parameter
%   Psicor      Roughness layer correction
%
% constants used (as global)
%   kappa       Von Karman's constant
%
% output:
%   zom         roughness lenght for momentum (m)
%   d           zero plane displacement (m)
% 

%% constants
global constants
kappa   = constants.kappa;

%% parameters
CR      = canopy.CR;
CSSOIL  = soil.CSSOIL;
CD1     = canopy.CD1;
Psicor  = canopy.Psicor;
LAI     = canopy.LAI;
h       = canopy.hc;

%% calculations
sq      = sqrt(CD1*LAI/2);
G1      = max(3.3, (CSSOIL + CR*LAI/2).^(-0.5));
d       = (LAI>1E-7 & h>1E-7).*h.*(1-(1-exp(-sq))./sq);          % Eq 12 in Verhoef et al (1997)
zom     = (h-d).*exp(-kappa*G1 + Psicor);