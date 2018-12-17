function [lE, H, ec, Cc]  = heatfluxes(ra,rs,Tc,ea,Ta,e_to_q,PSI,Ca,Ci)

global constants

rhoa = constants.rhoa;
cp   = constants.cp;
MH2O = constants.MH2O;
R    = constants.R;

% author: Dr. ir. Christiaan van der Tol (tol@itc.nl)
% date:     7 Dec 2007
% updated: 15 Apr 2009 CvdT     changed layout
% updated: 14 Sep 2012 CvdT     added ec and Cc to output
%
% parent: ebal.m
%
% usage:
% function [lE, H]  = heatfluxes(ra,rs,Tc,ea,Ta,e_to_q,PSI)
% 
% this function calculates latent and sensible heat flux
%
% input:
%   ra          aerodynamic resistance for heat         s m-1
%   rs          stomatal resistance                     s m-1
%   Tc          leaf temperature                        oC
%   ea          vapour pressure above canopy            hPa
%   Ta          air temperature above canopy            oC
%   e_to_q      conv. from vapour pressure to abs hum   hPa-1
%   PSI         leaf water potential                    J kg-1
%   Ca          ambient CO2 concentration               umol m-3
%   Ci          intercellular CO2 concentration         umol m-3
%
% output:
%   lEc         latent heat flux of a leaf              W m-2
%   Hc          sensible heat flux of a leaf            W m-2
%   ec          vapour pressure at the leaf surface     hPa
%   Cc          CO2 concentration at the leaf surface   umol m-3

lambda      = (2.501-0.002361*Tc)*1E6;  %      [J kg-1]  Evapor. heat (J kg-1)
ei          = equations.satvap(Tc).*exp(1E-3*PSI*MH2O/R./(Tc+273.15));  
qi          = ei .* e_to_q;
qa          = ea .* e_to_q;

lE          = rhoa./(ra+rs).*lambda.*(qi-qa);   % [W m-2]   Latent heat flux
H           = (rhoa*cp)./ra.*(Tc-Ta);           % [W m-2]   Sensible heat flux
ec          = ea + (ei-ea)*ra./(ra+rs);         % [W m-2] vapour pressure at the leaf surface
Cc          = Ca - (Ca-Ci).*ra./(ra+rs);        % [umol m-2 s-1] CO2 concentration at the leaf surface
