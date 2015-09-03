const = define_constants;
rhoa = const.rhoa;
Mair = const.Mair;

A = 10;
RH = [.01:.01:1];
gam = 0;
T = 20;

pmol2pm3 = rhoa/Mair*1E3;

Cs = 250 *pmol2pm3;

a = 3;
D0 = 10;
gs0 = 1E-6;


gs   = gs0 + a*A./(Cs-gam) .* (1+(1-RH).*satvap(T)./D0);

Ci   = Cs - A./gs;
close all

figure
plot(RH,gs)
figure
plot(RH,Ci/pmol2pm3)