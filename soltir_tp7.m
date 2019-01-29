function [wn,wl,Tall] = soltir_tp7(filename)

% soltir_tp7 Reads MODTRAN tp7 file and applies a new MIT algorithm to derive
% 18 spectral functions for atmospheric correction and simulations at BOA
% and TOA

% This function reads a MODTRAN tp7 output file that contains data for 4
% runs with different albedos a:
%
% For a = 0.5 at ground altitude
% For a = 1.0 at ground altitude
% For a = 0.0 at TOA
% For a = 1.0 at TOA
%
% From the MODTRAN outputs 18 atmospheric properties are derived, 
% named t1...t18
% The radiance data are converted into units of mW m-2 sr-1 nm-1
% 
% © Wouter Verhoef 2011-2014

fid     = fopen(filename,'r');
modname = filename(1:size(filename,2)-4);
outname = [modname '.atm'];

for i = 1:7, fgetl(fid); end

rline = str2num(fgetl(fid));
tts   = rline(4);
cts   = cosd(tts);

s     = fgetl(fid);
rline = sscanf(s,'%10f',3);
wns   = rline(1); wne = rline(2); wstep = rline(3);
nw    = int32((wne-wns)/wstep)+1;

datarec = zeros(nw,15,4);        % MODTRAN5.2.1 tp7 15 column output format
Tall    = zeros(nw,18);          % 18 output spectra

fgetl(fid); fgetl(fid);

for ipass = 1:4
    for il=1:nw
        s=fgetl(fid);
        dline=str2num(s);
        datarec(il,:,ipass)=dline;
    end
    for j = 1:12, fgetl(fid); end
end

wn  = datarec(:,1,1);
fac = wn.*wn;
wl  = 1e7./wn;
wls = wl(nw);
wle = wl(1);

% MIT algorithm supporting fluorescence retrievals
% Wout Verhoef Sept. 2011

% Support all applications, T-18 system
% OpT in heavy absorption bands now estimnated from Planck Tb at 6500 nm
% Wout Verhoef Oct. - Nov. 2012

% Constants of Planck function

c1 = 1.191066e-22;
c2 = 14388.33;

tran_boa    = datarec(:,2,1);
tran_toa    = datarec(:,2,3);

too         = tran_toa;

toasun      = datarec(:,14,4).*fac/pi*cts;  % = Eso cos(tts) / pi

%BOA

grfl50_boa  = datarec(:,7,1).*fac;
sfem50      = datarec(:,4,1).*fac;
sfem0_boa   = 2*sfem50;
grfl100_boa = datarec(:,7,2).*fac;
delgtot_boa = grfl100_boa-sfem0_boa;
crdd        = (grfl50_boa-sfem50)./(grfl100_boa-grfl50_boa-sfem50);

rdd         = max(0,1-crdd);

OpT         = crdd.*delgtot_boa./tran_boa;

% OpT at 6500 nm is used to estimate brightness temperature of La(b), which
% is used to get a minimum radiance where O is zero

wlp         = 6500;
Lp          = interp1(wl,OpT,wlp,'nearest');
Tp          = c2/(wlp*1e-3*log(1+c1*(wlp*1e-9)^(-5)/Lp));
Lmin        = c1*(wl*1e-9).^(-5)./(exp(c2./(wl*1e-3*Tp))-1);

%TOA

grfl100_toa = datarec(:,7,4).*fac;
sfem0       = datarec(:,4,3).*fac;
delgtot_toa = grfl100_toa - sfem0;
OpTtran     = crdd.*delgtot_toa;
path100     = datarec(:,5,4).*fac;
path0       = datarec(:,5,3).*fac;

rso         = path0./toasun;

delpath     = path100 - path0;
ptem100     = datarec(:,3,4).*fac;
ptem0       = datarec(:,3,3).*fac;
delptem     = max(0,ptem100-ptem0);
delatmo     = delpath + delptem;

iT          = (wl > 4600);
ia          = (~iT & delpath == 0) | (iT & delptem == 0);

fO          = delpath./(delatmo+1e-60).*~ia + ~iT.*ia;
fT          = delptem./(delatmo+1e-60).*~ia +  iT.*ia;

O           = fO.*OpT;
T           = fT.*OpT;

% Correct cases where T = 0

i0          = (T == 0);
T(i0)       = Lmin(i0);
fT(i0)      = T(i0)./OpT(i0);
Ttran       = fT.*OpTtran;
Otran       = OpTtran-Ttran;

tdo         = delatmo./(delgtot_toa+1e-6).*tran_toa;

tdo(ia)     = 0;

gsun100_toa = datarec(:,8,4).*fac;
gsun100_boa = datarec(:,8,2).*fac;

tsstoo      = gsun100_toa./toasun;
tss         = gsun100_boa./toasun;
tsdM        = O./toasun-tss;

tsdM(ia)    = 0;

% Apply log regression to predict tsd from tdo, too, tss, and wl

Ir          = (wl>1500 & wl<1700) | (wl>2100 & wl<2300) | ...
               (wl>3950 & wl<4100);
y           = log(tsdM(Ir)./tss(Ir))-log(tdo(Ir)./too(Ir));
x           = log(wl(Ir));
n           = size(x,1);
xm          = sum(x)/n;
ym          = sum(y)/n;
a           = (x-xm)'*(y-ym)/((x-xm)'*(x-xm));
b           = ym-a*xm;
p           = a*log(wl)+b;
tsdp        = tss.*tdo./(too+1e-60).*exp(p);

% weight proportional to delatmo squared

wgtM        = delatmo.^2;

tsd         = (wgtM.*tsdM+tsdp)./(wgtM+1);

fsun        = (tss+1e-60)./(tss+1e-60+tsd);
iH          = (wl>2600);

tsdtoo      = Otran./toasun-tsstoo;
tsdtoo(iH)  = tsd(iH).*too(iH);

tssrddtoo   = tsstoo.*rdd;
toordd      = too.*rdd;
tssrdd      = tss.*rdd;
tsstdo      = fsun.*delpath.*crdd./toasun;
tsdtdo      = (1-fsun).*delpath.*crdd./toasun;
Lat         = ptem0-sfem0_boa./crdd.*tdo;
Lab         = T+sfem0_boa.*crdd;
Labtoo      = Ttran+crdd.*sfem0_boa.*too;
Labtdo      = crdd.*(delptem+sfem0_boa.*tdo);

Tall        = [toasun rso rdd tss tsd too tdo tsstoo tsdtoo tsstdo tsdtdo ...
                  tssrdd toordd tssrddtoo Lat Lab Labtoo Labtdo];

fclose(fid);

% Verification against MODTRAN

% pt0   = Lat;
% pt100 = Lat+Labtdo./(1-rdd);
% pa0   = rso.*toasun;
% pa100 = toasun.*(rso+(tsstdo+tsdtdo)./(1-rdd));
% gb50  = (toasun*.5.*(tss+(2*tsd+tssrdd)./(2-rdd))+Lab./(2-rdd)).*tran_boa;
% gb100 = (toasun.*(tss+tsd)./(1-rdd)+Lab./(1-rdd)).*tran_boa;
% gt100 = toasun.*(tsstoo+(tsdtoo+tssrddtoo)./(1-rdd))+Labtoo./(1-rdd);
% gs100 = toasun.*tsstoo;

% Write data to output file

fid=fopen(outname,'w');
str1=' WN (cm-1) WL (nm)  '; 
str2='      T1       T2       T3        T4        T5        T6        T7      ';
str3='  T8        T9        T10       T11       T12       T13       T14      ';
str4=' T15          T16          T17          T18';
str5=['                       toasun     rso      rdd       tss       tsd       too' ...
 '       tdo      tsstoo    tsdtoo    tsstdo    tsdtdo    tssrdd'  ...
 '    toordd  tssrddtoo     Lat          Lab         Labtoo       Labtdo'];
str=[str1 str2 str3 str4];
fprintf(fid,'%s\r\n',str);
fprintf(fid,'%s\r\n\r\n',str5);

for i = 1:nw
    str = sprintf('%9.2f',wn(i));
    str = [str sprintf('%10.3f',wl(i))];
    str = [str sprintf('%10.5f',Tall(i,1))];
    str = [str sprintf('%10.6f',Tall(i,2:14));];
    a   = sprintf('%14.6e',Tall(i,15:18));
    str = [str a];
    fprintf(fid,'%s\r\n',str);
end

fclose('all');

end

