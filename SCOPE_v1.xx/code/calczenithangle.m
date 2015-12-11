function [Fi_s,Fi_gs,Fi_g,Omega_s] = calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat)
%
% author: Christiaan van der Tol (c.vandertol@utwente.nl)
% date:     Jan 2003
% update:   Oct 2008 by Joris Timmermans (j_timmermans@itc.nl): 
%               - corrected equation of time
%           Oct 2012 (CvdT) comment: input time is GMT, not local time!
%
% function [Fi_s,Fi_gs,Fi_g]= calczenithangle(Doy,t,Omega_g,Fi_gm,Long,Lat)
%
% calculates pi/2-the angle of the sun with the slope of the surface.
%
% input:
% Doy       day of the year
% t         time of the day (hours, GMT)
% Omega_g   slope azimuth angle (deg)
% Fi_gm     slope of the surface (deg)
% Long      Longitude (decimal)
% Lat       Latitude (decimal)
%
% output:
% Fi_s      'classic' zenith angle: perpendicular to horizontal plane
% Fi_gs     solar angle perpendicular to surface slope
% Fi_g      projected slope of the surface in the plane through the solar beam and the vertical
% 

%parameters (if not already supplied)
if nargin<6
    Long        =   13.75;                      % longitude
    Lat         =   45.5;                       % latitude
    if (nargin<4)
        Omega_g =   210;                        % aspect
        Fi_gm   =   30;                         % slope angle
    end
end

%convert angles into radials
G               =   (Doy-1)/365*2*pi;           % converts day of year to radials
Omega_g         =   Omega_g/180*pi;             % converts direction of slope to radials
Fi_gm           =   Fi_gm/180*pi;               % converts maximum slope to radials
Lat             =   Lat/180*pi;                 % converts latitude to radials

%computes the declination of the sun
d               =   0.006918-0.399912*cos(G  )+ 0.070247*sin(G  )- ...
                     0.006758*cos(2*G)+ 0.000907*sin(2*G)- ...
                     0.002697*cos(3*G)+ 0.00148*sin(3*G);
                                
%Equation of time
Et              =   0.017 + .4281 * cos(G) - 7.351 * sin(G) - 3.349 * cos(2*G) - 9.731 * sin(2*G);

%computes the time of the day when the sun reaches its highest angle                                
tm              =   12+(4*(-Long)-Et)/60;      % de Pury and Farquhar (1997), Iqbal (1983)

%computes the hour angle of the sun
Omega_s         =   (t-tm)/12*pi;

%computes the zenithangle (equation 3.28 in De Bruin)
Fi_s            =   acos(sin(d)*sin(Lat)+cos(d)*cos(Lat).*cos(Omega_s));

%computes the slope of the surface Fi_g in the same plane as the solar beam
Fi_g            =   atan(tan(Fi_gm).*cos(Omega_s-Omega_g));

%computes the angle of the sun with the vector perpendicular to the surface
Fi_gs           =   Fi_s + Fi_g;