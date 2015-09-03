[alpha,epsmax] = deal(zeros(81,2));
Pv = zeros(6,2);
P = zeros(5,2);
fm = zeros(901,2);

for t = 1:2
    
    target          = [0,1]; %0: SIF, 1: A
    wl              = [685 -999];
    scaling         = [1,.73];

    
    [Cabi,epsmax(:,t),alpha(:,t)]   = Cab_LAI_scaling(target(t),wl(t),scaling(t));
    [Pv(:,t),~,~]          = Vcmo_scaling(target(t),wl(t),0);
    [mu,fm(:,t),s7]               = m_Ascaling(epsmax(:,t),alpha(:,t),Pv(:,t),target(t),wl(t));
    [P(:,t),~,s3]                 = tts_scaling(target(t),wl(t),0);
end
%Dir = '../output/Bridge_database_2015-05-15-1011/';
Dir = '../output/PARCS_database_2015-05-18-1534/'; %Cca varies with Cab
%Dir = '../output/PARCS_database_2015-05-18-1541/'; %Cca is kept constant
Q = dlmread([Dir 'fluxes.dat'],'',2,0);                 % fluxes
R = dlmread([Dir 'radiation.dat'],'',2,0);                 % fluxes
F = dlmread([Dir 'fluorescence.dat'],'',2,0);                 % fluxes

D = dlmread('fluxes_new.txt','',2,0);
p = dlmread('Database BRIDGE inputs.txt','',2,0);
%SIF = Q(:,end-1);
SIF = F(:,wl(1)-639);
A   = Q(:,11);
Cab = p(:,2);
m   = p(:,8);
LAI = p(:,9);
tts = p(:,18);
Rin = R(:,end);%645*ones(length(Cab),1);
Ta = p(:,12);
ea = p(:,14);
rh = ea./satvap(Ta);

epsmaxi  = epsmax(round(Cab)+1,1);
alphai   = alpha(round(Cab)+1,1);
ftts        = calc_ftts(P(:,1),tts);
index2  = min(901,max(1,round(m.*rh*100)-99));
fmi      = fm(index2,1);

fLAI    = 1-exp(-alphai.*LAI);
SIF0        = Rin.*epsmaxi.*fLAI.*fmi.*ftts;
%fVcmoi      = SIF./SIF0;

x = (0:.01:2)';
%%

Vcmoin = p(:,7);
% xi = Vcmo./Rin;
%index = min(201,max(1,round(sqrt(Vcmoin./(Rin*2.2039))*100+1)));
%fVcmoiS = fVcmo(index,1);
%fVcmoiA = fVcmo(index,2);
fVcmoiS = calcfVcmo(Pv(:,1),sqrt(Vcmoin./(Rin*2.2039)));
%fVcmoiA = calcfVcmo(Pv(:,2),sqrt(Vcmoin./(Rin*2.2039)));

SIFmod = SIF0.*fVcmoiS;

fVcmoiS = SIF./SIF0;
%Vcmo,fVcmoiA] = deal(0*Cab);

xi = zeros(length(A),1);
for k = 1:length(A)
    xi(k) = lsqnonlin(@(xa) calcfVcmo(Pv(:,1),xa)- fVcmoiS(k), .3 , .01, .5);
end
fVcmoiA = calcfVcmo(Pv(:,2),xi) ;
Vcmo = xi.^2.*Rin*2.2039;

% for k = 1:length(Cab)
%     [~,I] = min(abs(fVcmo(1:50,1)-fVcmoiS(k)));
%     Vcmo(k) = x(I).^2.*Rin(k);
%     fVcmoiA(k) = fVcmo(I,2);
%     
% end

epsmaxiA  = epsmax(round(Cab)+1,2);
alphaiA = alpha(round(Cab)+1,2);
fttsA = calc_ftts(P(:,2),tts);
fmiA      = fm(index2,2);
fLAIA   = 1-exp(-alphaiA.*LAI);

Amod       = Rin.*epsmaxiA.*fLAIA.*fVcmoiA.*fttsA.*fmiA;
SIFmod2 = SIF0.*fVcmoiS;

%%
% figure(1), hold on
% subplot(211), plot(SIF,'x'), hold on, plot(SIFmod,'c'), plot(SIFmod2,'m')
% subplot(212), plot(A,'x'), hold on, plot(Amod,'r')
% 
% figure(12), hold on
% plot(A,Amod,'rx'), hold on
% plot([0 40],[0,40],'k')

F9 = figure(9);clf
set(F9,'Position',[473 153 685 532])

subplot(311)
plot(SIF,'k'), hold on
plot(SIFmod,'bx-')
legend('SCOPE','Eqs 1&2','Orientation','horizontal','location','best')
ylabel('SIF760 (Wm^{-2}\mum^{-1}sr^{-1})')
set(gca,'xtick',[1:1:40],'FontSize',8)

subplot(312)
plot(Vcmoin,'k'), hold on
plot(Vcmo,'bx-')
legend('SCOPE','Retrieved','Orientation','horizontal','location','best')

I = [17,18,29,32,33]';
plot(I,Vcmo(I),'ro')
text(17,350,'Planophile')
text(18,20,'Erectophile')
text(25,130,'Tropical')
text(32,150,'Plagiophile')
text(33,100,'Erectophile')
set(gca,'xtick',[1:1:40],'FontSize',8)
ylabel('V_{cmo} (\mumol m^{-2}s^{-1})')

subplot(313)
plot(A,'k'), hold on
plot(Amod,'bx-')
set(gca,'xtick',[1:1:40],'FontSize',8)
ylabel('Photosynthesis (\mumol m^{-2}s^{-1})')
legend('SCOPE','Retrieved','Orientation','horizontal','location','best')