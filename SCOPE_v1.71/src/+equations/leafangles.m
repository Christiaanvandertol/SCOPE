function [lidf]=  leafangles(a,b)                                     
% Subroutine FluorSail_dladgen
% Version 2.3 
% For more information look to page 128 of "theory of radiative transfer models applied in optical remote sensing of
% vegetation canopies"
%
% FluorSail for Matlab
% FluorSail is created by Wout Verhoef, 
% National Aerospace Laboratory (NLR)
% Present e-mail: w.verhoef@utwente.nl
%
% This code was created by Joris Timmermans, 
% International institute for Geo-Information Science and Earth Observation. (ITC)
% Email: j.timmermans@utwente.nl
%
%% main function
F           =   zeros(1,13);
for i=1:8                                                               
    theta   =   i*10;                  %                theta_l =  10:80
    F(i)    =   dcum(a,b,theta);      %                 F(theta)
end

for i=9:12                                                              
    theta   =   80 + (i-8)*2;                         % theta_l = 82:88
    F(i)    =   dcum(a,b,theta);                     %  F(theta)
end

for i=13:13                                          %  theta_l = 90:90
    F(i) =   1;                                      %  F(theta)
end

lidf        =   zeros(13,1);
for i=13:-1:2                                                           
    lidf(i) =   F(i) -   F(i-1);                     %  lidf   =   dF/dtheta;
end
lidf(1) =   F(1);                                    %  Boundary condition

%% SubRoutines
function [F]   =  dcum(a,b,theta)
rd  =   pi/180;                                     %   Geometrical constant
if a>1 
    F    =   1 - cos(theta*rd);
else
    eps     =   1e-8;
    delx    =   1;
    
    x       =   2*rd *theta;
    theta2  =   x;
                                                                        %    
    while max(delx > eps)
        y   =   a*sin(x) + 0.5*b*sin(2*x);
        dx  =   0.5*(y - x + theta2);
        x   =   x + dx;
        delx=   abs(dx);
    end
    F    =   (2*y + theta2)/pi;                     %   Cumulative leaf inclination density function
    %pag 139 thesis says: F = 2*(y+p)/pi. 
    %Since theta2=theta*2 (in rad), this makes F=(2*y + 2*theta)/pi    
end