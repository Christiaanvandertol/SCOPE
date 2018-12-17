function Tav = tav(alfa,nr)
n2                                  =   nr.^2;
np                                  =   n2 + 1;
nm                                  =   n2 - 1;

% Stern's formula in Lekner & Dorf (1988) gives reflectance for alfa = 90 degrees

% y1 = (3*n2+2*nr+1)./(3*(nr+1).^2);
% y2 = 2*nr.^3.*(nr.^2+2*nr-1)./(np.^2.*nm);
% y3 = n2.*np.*log(nr)./nm.^2;
% y4 = n2.*nm.^2./np.^3.*log((nr.*(nr+1)./(nr-1)));

% st = y1-y2+y3-y4;

a                                   =   +((nr+1).^2)/2;
k                                   =   -((n2-1).^2)/4;
sin_a                               =   sind(alfa);
%
if alfa~=0    
    B2                              =   sin_a^2 - np/2;
    B1                              =   (alfa~=90) * sqrt( B2.^2 + k );
   
    b                               =   B1 - B2;
    b3                              =   b.^3;
    a3                              =   a.^3;
    
    ts                              =   (k.^2./(6*b3) + k./b - b./2) - ...
                                        (k.^2./(6*a3) + k./a - a./2);
                                    
    tp1                             =   -2*n2.*    (   b  -  a   ) ./ (np.^2);
    tp2                             =   -2*n2.*np.*(  log(b./a)  ) ./ (nm.^2);
    tp3                             =      n2.*    ( 1./b - 1./a ) ./ 2; 
    
%     tp4                             =   16*n2.^2.* (n2.^2+1) .* ( log(2*np.*b - nm.^2) - log(2*np.*a - nm.^2) ) ./ (np.^3.*nm.^2);    
%     tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./ (np.^3       );

    tp4                             =	16*n2.^2.* (n2.^2+1) .* ( log((2*np.*b - nm.^2)./(2*np.*a - nm.^2))  ) ./(np.^3.*nm.^2);
    tp5                             =   16*n2.^2.* (n2     ) .* ( 1./(2*np.*b - nm.^2) - 1./(2*np.*a - nm.^2)) ./(np.^3);							 
    tp                              =   tp1 + tp2 + tp3 + tp4 + tp5;
    Tav                             =   (ts + tp) ./ (2*sin_a.^2);
else
    Tav                             =   4 *nr/((nr+1)*(nr+1));
end
return