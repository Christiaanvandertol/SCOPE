function [T_C]  =   Brightness_T(H)

global constants
sigmaSB = constants.sigmaSB;
C2K     = constants.C2K;

T_C             = (H/sigmaSB).^(1/4)-C2K;
end