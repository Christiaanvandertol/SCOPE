%% this is the directory where the LUT is located
Dir = '../output/C3_Cab_Vcmo_Rin_sens_2015-05-06-1050/';


%% the files of the lookup table are loaded
Q = dlmread([Dir 'fluxes.dat'],'',2,0);                 % fluxes
p = dlmread([Dir 'pars_and_input_short.dat'],'',1,0);   % parameters
f = dlmread([Dir 'fluorescence.dat'],'',2,0);   % parameters

%%
Avec        = Q(:,11);              % photosynthesis

Cab         = p(:,1);               % leaf area index
Vcmo        = p(:,2);               % incident shortwave light (W m-2)
Rin         = p(:,3);
aPAR        = Q(:,17);
faPAR       = Q(:,19);
iPAR        = aPAR./faPAR;

Cabunique = unique(Cab);            % unique values in the LUT
Rinunique = unique(Rin);            % unique values in the LUT
Vcmounique = unique(Vcmo);

% place the photosynthesis data in a matrix instead of a vector
A = zeros(length(Cabunique),length(Rinunique).*length(Vcmounique)); 
A(:) = Avec(:)./Rin(:);

%Anorm = A./(ones(length(A),1)*max(A));
Anorm = Avec./Rin;
% Anorm2 = mean(A./(max(A,[],2)*ones(1,136)));


% initializatios
Amod = NaN*A;
xi = zeros(length(Rinunique),4);
gamma = .1;



Anorm2 = Anorm;

for k = 1:length(Cabunique)
    I = find(Cab==Cabunique(k));
    Anorm2(I) = Anorm(I)./max(Anorm(I));
end

F2 = figure(2); clf
set(F2,'Position',[360 461 615 461])

s2(1) = subplot(221);
plot(Vcmo,Avec,'kx');
xlabel('V_{cmo} (\mumol m^{-2}s^{-1})')
ylabel('SIF(740) (W m^{-2}\mum^{-1}sr^{-1})')
set(gca,'xlim',[0 220])

% s2(2) = subplot(132);
% plot(Vcmo,Anorm,'kx');
% xlabel('V_{cmo} (\mumol m^{-2}s^{-1})')
% ylabel('SIF(740) (W m^{-2}\mum^{-1}sr^{-1})')

s2(2) = subplot(222);
plot(Vcmo,Anorm2,'kx','MarkerSize',3)
ylabel('f_{Vcmo}')
xlabel('V_{cmo}')


s2(3) = subplot(223);
plot(Vcmo./iPAR,Anorm2,'kx','MarkerSize',3)
ylabel('f_{Vcmo}')
xlabel('V_{cmo}/iPAR')

s2(4) = subplot(224);
plot(log(Vcmo./iPAR),Anorm2,'kx','MarkerSize',3), hold on
[xc,yc] = consolidator(log(Vcmo./iPAR),Anorm2,@mean);

logVcmoi = (-6:.01:2)';
yy = interp1(xc,yc,logVcmoi);
yyVcmo = smooth(yy,30);
plot(logVcmoi,yyVcmo,'k','LineWidth',2)
plot([log(.22) log(.22)],[.5,1],'k--')
xlabel('log(V_{cmo}/iPAR)')
%set(gca,'yticklabel','')
text(-.9,.6,'light limited','FontSize',8)
text(-9,.6,'light saturated','FontSize',8)
ylabel('f_{Vcmo}')

resizefigure(s2,2,2,.1,.14,.08,.1);