fieldsite= 2;
if fieldsite == 1
    M               = load(['drivingvarsplot1.dat']);
else
    M               = load(['drivingvarsplot2.out']);
end
M(M<-900) = NaN;

t_                  = M(:,1);
Ta200_              = M(:,2);
Ta400_              = M(:,3);
ea200_              = M(:,4);
ea400_              = M(:,5);
Rin_                = M(:,6);
Rli_                = M(:,7);
u_                  = M(:,8);
SMC_                = M(:,9);
Ca_                 = M(:,10);

jj                  = ~isnan(Ta200_);                       %find data with good quality Ta200 data
[ea_,Ta_]           = deal(NaN*(ones(size(t_))));
ea_(jj)             = ea200_(jj);
Ta_(jj)             = Ta200_(jj);

jj                  = ~isnan(Ta400_);                       %find data with good quality Ta400 data
Ta_(jj)             = Ta400_(jj);
ea_(jj)             = ea400_(jj);
p_                  = 1E3*ones(size(Ta_));
year_               = 2002*ones(size(Ta_));

save    't_.dat'    t_ -ascii
save    'Ta_.dat'   Ta_ -ascii
save    'ea_.dat'   ea_ -ascii
save    'Rin_WRONG.dat'  Rin_ -ascii
save    'Rli_WRONG.dat'    Rli_ -ascii
save    'u_.dat'    u_ -ascii
save    'SMC_.dat'  SMC_ -ascii
save    'CO2_.dat'  Ca_ -ascii
save    'p_.dat'    p_ -ascii
save    'year_.dat'    year_ -ascii