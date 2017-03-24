if ~exist('fieldsite','var'), fieldsite = 1; end

if fieldsite ==1, plotdir = 'plot1'; else plotdir = 'plot2'; end

MCO2      = 44;         % [g mol-1]     molecular mass of carbon dioxide


if fieldsite == 1;
   V       = load('verivarsplot1.out');
   Tveg       = load('plot1/Tveg1.dat');
else
	V       = load('verivarsplot2.out');
    Tveg       = load('plot2/Tveg2.dat');
end
%Tc_         = mean(Tveg(:,2:end),2);

V(V<-900)    = NaN;

% load validation data
t       = V(:,1);
Rso     = V(:,2);
Rlo     = V(:,3);
Rn4c	  = V(:,4);
RnLite  = V(:,5);
Rn_bottom1 = V(:,6);
Rn_bottom2 = V(:,7);
Rn_bottom3 = V(:,8);
TsIRT   = V(:,9);
H       = V(:,11);
lE      = V(:,12);
A       = -V(:,13)/MCO2*1E3;
G1      = V(:,14);
G_rem   = V(:,15);
Ta_can  = V(:,16);
ea_can  = V(:,17);

if fieldsite == 1;
    Ts11    = V(:,18);
    Ts12    = V(:,20);
    Ts13    = V(:,21);
    Ts14    = V(:,22);
    Ts01    = V(:,19);
    Ts02    = V(:,23);
    Ts03    = V(:,24);
    Ts04    = V(:,25);
else
    Ts11    = V(:,19);
    Ts12    = V(:,20);
    Ts13    = V(:,23);
    Ts14    = V(:,24);
    Ts01    = V(:,18);
    Ts02    = V(:,21);
    Ts03    = V(:,22);
    Ts04    = V(:,25);
end

if fieldsite ==2
    lE      = NaN*lE;
    H       = NaN*H;
    A       = NaN*A;
end
Ts0     = mean([Ts01, Ts02, Ts03, Ts04],2);
Ts1     = mean([Ts11, Ts12, Ts13, Ts14],2);
Tc      = mean(Tveg(:,2:end),2);
Rn_bot  = mean([Rn_bottom1 Rn_bottom2 Rn_bottom3],2);

cd(plotdir)

save t.dat          t           -ascii
save Rso.dat        Rso         -ascii
save Rlo.dat        Rlo         -ascii
save Rn4c.dat       Rn4c        -ascii
save RnLite.dat     RnLite      -ascii
save Rn_bottom1.dat Rn_bottom1  -ascii
save Rn_bottom2.dat Rn_bottom2  -ascii
save Rn_bottom3.dat Rn_bottom3  -ascii
save Rn_bottom.dat  Rn_bot      -ascii
save TsIRT.dat      TsIRT       -ascii
save H.dat          H           -ascii
save lE.dat         lE          -ascii
save A.dat          A           -ascii
save G.dat          G1          -ascii
save G_rem.dat      G_rem       -ascii
save Ta_can.dat     Ta_can      -ascii
save ea_can.dat     ea_can      -ascii
save Ts0.dat        Ts0         -ascii
save Ts1.dat        Ts1         -ascii
save Tc.dat         Tc          -ascii
cd ..