function [es,s] = satvap(T)

%% function [es,s]= satvap(T) 
% Author: Dr. ir. Christiaan van der Tol
% Date: 2003
%
% calculates the saturated vapour pressure at 
% temperature T (degrees C)
% and the derivative of es to temperature s (kPa/C)
% the output is in mbar or hPa. The approximation formula that is used is:
% es(T) = es(0)*10^(aT/(b+T));
% where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
% and s(T) = es(T)*ln(10)*a*b/(b+T)^2

%% constants
a           = 7.5;
b           = 237.3;         %degrees C

%% calculations
es          = 6.107*10.^(7.5.*T./(b+T));
s           = es*log(10)*a*b./(b+T).^2;