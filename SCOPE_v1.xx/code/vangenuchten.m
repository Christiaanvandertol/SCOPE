function out = vangenuchten(input,thetares, thetasat, alpha,n,option)
%h = vangenuchten(input,thetares, thetasat, alpha,n,option);
%if option not specified, or option <>1, h = input, and theta is calculated, otherwise theta = input, and h is calculated


m = 1-1/n;
%m = 1;

if nargin>5
    if option ==1
        theta   = input;
        Se      = (theta - thetares)/(thetasat - thetares);
        out     = -1/alpha*(Se.^(-1/m)-1).^(1/n);
    end
else
    h       = input;
    out     = thetares + (thetasat-thetares)./(1+abs(alpha*h).^n).^m;
end