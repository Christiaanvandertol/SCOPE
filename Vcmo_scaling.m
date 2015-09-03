function [P,fVcmo,x,y] = Vcmo_scaling(target,wl,doplot)
%% this is the directory where the LUT is located
%Dir = '../output/C3_Cab_Vcmo_Rin_sens_2015-05-06-1050/'; %Cca = 20
Dir = '../output/C3_Cab_Vcmo_Rin_sens_2015-05-18-1612/'; %Cca varies with Cab
%Dir = '../output/C3_Cab_Vcmo_Rin_sens_2015-05-11-1319/';

%% the files of the lookup table are loaded
Q = dlmread([Dir 'fluxes.dat'],'',2,0);                 % fluxes
p = dlmread([Dir 'pars_and_input_short.dat'],'',1,0);   % parameters
f = dlmread([Dir 'fluorescence.dat'],'',2,0);   % parameters
fhem = dlmread([Dir 'fluorescence_hemis.dat'],'',2,0);

%%
spi = 2*target +(wl==685);

switch target
    case 0
        if wl>0
            Avec        = f(:,wl-639);              % fluorescence
        else
            Avec = Q(:,21);
%             Avec = zeros(size(fhem,1),1);
%             for k = 1:size(fhem,1)
%                 Avec(k)        = 0.001 * Sint(fhem(k,:),(640:850));
%             end
        end
    otherwise
        Avec        = Q(:,11);              % photosynthesis
end

Cab         = p(:,1);               % leaf area index
Vcmo        = p(:,2);               % incident shortwave light (W m-2)
%Rin         = p(:,3);
aPAR        = Q(:,17);
faPAR       = Q(:,19);
iPAR        = aPAR./faPAR;

Cabunique = unique(Cab);            % unique values in the LUT
Anorm = Avec./iPAR;
Anorm2 = Anorm;
for k = 1:length(Cabunique)
    I = find(Cab==Cabunique(k));
    Anorm2(I) = Anorm(I)./max(Anorm(I));
end

%x = sqrt(Vcmo./iPAR);
x = (Vcmo./iPAR);
P = lsqnonlin(@(p) Anorm2 - (p(1)+(p(2)*x.^p(6)).*(1./(p(4)+p(3)*x.^p(5)))),[1 1 1 1 3 2],[-1 -1 -1 -1 -1 -2],[10 10 10 10 5 5] );
fVcmo = P(1)+(P(2)*x.^P(6)).*(1./(P(4)+P(3)*x.^P(5)));
x = (0:.01:4)';
y = P(1)+(P(2)*x.^P(6)).*(1./(P(4)+P(3)*x.^P(5)));

%%
if nargin>2
    if doplot
        F4 = figure(4+spi); clf
        set(F4,'Position',[360 461 615 461])
        
        s2(1) = subplot(221);
        plot(Vcmo,Avec,'kx');
        xlabel('V_{cmo} (\mumol m^{-2}s^{-1})')
        if target == 0
            switch wl
                case 685, ylabel('SIF685 (W m^{-2}\mum^{-1}sr^{-1})')
                case 760, ylabel('SIF760 (W m^{-2}\mum^{-1}sr^{-1})')
            end
        else
            ylabel('A (\mumol m^{-2} s^{-1})')
        end
        set(gca,'xlim',[0 220])
        
        s2(2) = subplot(222);
        plot(Vcmo,Anorm2,'kx','MarkerSize',3)
        ylabel('f_{Vcmo}')
        xlabel('V_{cmo}')
        
        s2(3) = subplot(223);
        plot(Vcmo./iPAR,Anorm2,'kx','MarkerSize',3), hold on
        ylabel('f_{Vcmo}')
        xlabel('V_{cmo}/iPAR')
        set(gca,'ylim',[0 1])
        plot(x,y,'k')
        
        s2(4) = subplot(224);
        %plot(sqrt(Vcmo./iPAR),Anorm2,'kx','MarkerSize',3), hold on
        plot(log(Vcmo./iPAR),Anorm2,'kx','MarkerSize',3), hold on
        
        plot(log(x),y,'k')
        
        %xlabel('sqrt(V_{cmo}/iPAR)')
        xlabel('log(V_{cmo}/iPAR)')
        
        ylabel('f_{Vcmo}')
        %set(gca,'xlim',[0,2.5],'ylim',[0,1])
        set(gca,'xlim',[-4,2],'ylim',[0,1])
        
        resizefigure(s2,2,2,.1,.14,.08,.1);
    end
end