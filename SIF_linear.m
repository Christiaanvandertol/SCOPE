% Cabrand = 70*rand(1,100);
% LAIrand = 7*rand(1,100);
% Rinrand = 800*rand(1,100);
% Vcmorand = 150*rand(1,100);
% 
% Cdmrand = .01+.01*rand(1,100);
% Nrand = 1+2*rand(1,100);
% Trand = 10+20*rand(1,100);
% mrand = 2+10*rand(1,100);
% ttsrand = 30+20*rand(1,100);
% earand = .5*satvap(Trand);
% 
% 
% xlswrite('../input_data.xlsx',Cabrand,'inputdata','B9');
% xlswrite('../input_data.xlsx',LAIrand,'inputdata','B45');
% xlswrite('../input_data.xlsx',Rinrand,'inputdata','B53');
% xlswrite('../input_data.xlsx',Vcmorand,'inputdata','B19');
% 
% xlswrite('../input_data.xlsx',Cdmrand,'inputdata','B11');
% xlswrite('../input_data.xlsx',Nrand,'inputdata','B14');
% xlswrite('../input_data.xlsx',Trand,'inputdata','B54');
% xlswrite('../input_data.xlsx',mrand,'inputdata','B20');
% xlswrite('../input_data.xlsx',ttsrand,'inputdata','B82');
% xlswrite('../input_data.xlsx',earand,'inputdata','B57');
% 
% %xlswrite('../input_data.xlsx','forward_random_input' ,'filenames','B4');
% 
% %%
% 
% SCOPE

%%
% Vcmo_SIFscaling
% Cab_LAI_SIFscaling

Vcmo_Ascaling
Cab_LAI_Ascaling
m_Ascaling

%%
F5 = figure(5); clf, hold on, z = zeros(2,1);
set(F5,'Position',[375 380 321 241])
for k = 1:2
    switch k
        case 1, outputdir = '../output/C3_Cab_LAI_sens_2015-05-06-1426/';
        case 2, outputdir = '../output/forward_random_input_2015-05-06-1537/';
    end
    
    
    Q = dlmread([outputdir 'fluxes.dat'],'',2,0);                 % fluxes
    p = dlmread([outputdir 'pars_and_input_short.dat'],'',1,0);   % parameters
    f = dlmread([outputdir 'fluorescence.dat'],'',2,0);   % parameters
    SIF740 = Q(:,11);%f(:,740-649);
    
    
    Cab     = p(:,1);
    Vcmo	= p(:,2+(k==2)*2);%p(:,2);
    LAI     = p(:,3+(k==2)*3);%p(:,3);
    Rin     = p(:,4+(k==2)*3);%p(:,4);
    
    aPAR        = Q(:,17);
    faPAR       = Q(:,19);
    iPAR        = aPAR./faPAR;
    
    
    %
    % Cab = 30;
    % LAI = 3;
    % Rin = 400;
    % Vcmo = 20;
    %
    %
    index = min(800,max(1,round(log(Vcmo./Rin)*100)+600));
    epsmax  = yyepsmax(round(Cab)+1);
    alpha   = yyalpha(round(Cab)+1);
    fLAI    = 1-exp(-alpha.*LAI);
    fVcmo   = yyVcmo(index);
    
    SIF740_l = Rin.*epsmax.*fLAI.*fVcmo;
    
    z(k) = plot(SIF740_l,SIF740,'ko','MarkerFaceColor',[.5 .5 .5],'MarkerSize',4);

end
plot([0 2.2],[0,2.2],'k')
set(z(2),'MarkerEdgeColor','k','MarkerFaceColor','none')
xlabel('SIF(740) - SCOPE_{reduced}')
ylabel('SIF(740) - SCOPE_{full}')