
%Cabscaling
%LAIscaling

Rin = 600.46;
roundRin = min(max(1,round(Rin)),1000);

alpha = yyalpha(roundRin);
beta  = yybeta(roundRin);


%% 
Direc = '';

LAI = dlmread([Direc 'LAI_formap.asc'],'',6,0);
Cab = dlmread([Direc 'Lcc_formap.asc'],'',6,0);

LAI(LAI==-9999) = NaN;
Cab(Cab==-9999) = NaN;

LAIX = 4;
CabX = 40;
AX = 30;

LT = 1-exp(-alpha*LAI);
CT = Cab*(1+beta)./(Cab+beta);

LTX = 1-exp(-alpha*LAIX);
CTX = CabX*(1+beta)./(CabX+beta);

A = AX*LT/LTX .* CT/CTX;

F3 = figure(3); clf, set(F3,'Position',[360 231 316 691])
s3(1) = subplot(311);
imagesc(Cab), title('Cab (\mug cm^{-2})'), colorbar

s3(2) = subplot(312);
imagesc(LAI), title('LAI '), colorbar

s3(3) = subplot(313);
imagesc(A), title('A (\mumol m^{-2} s^{-1})'), colorbar
