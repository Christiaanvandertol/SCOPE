% infor   = envihdrread('../Fluorescence/Fs_FLD3_20120823-SEL-1350-0600-L2-FLUO_subset1_destriped_radiance_georef_new method.hdr');
% 
% Fm     = envidataread('../Fluorescence/Fs_FLD3_20120823-SEL-1350-0600-L2-FLUO_subset1_destriped_radiance_georef_new method.bsq',infor);
% Fm = Fm*1E3;
% F0 = load('F0.dat');
% F0 = F0(1:6221,92:end);
% 
% LC = dlmread('../Land use class/mnfsel1350FLUO_roi1_potato0_4_corn0_4_ENVIstd.txt','',5,0);
% LC = LC(1:6221,1:485);
% 
% 
% [x,x2,c] = deal(zeros(8,1));
% for j = 1:length(Fm(:))
%     x(LC(j)+1)  = x(LC(j)+1) + Fm(j) * (Fm(j)>-1 & Fm(j)<10);
%     x2(LC(j)+1) = x2(LC(j)+1) + Fm(j).*Fm(j)* (Fm(j)>-1 & Fm(j)<10);
%     c(LC(j)+1)  = c(LC(j)+1) + 1* (Fm(j)>-1 & Fm(j)<10);
% end
% % 
% Fcl = x./c;
% stdFcl = 1./(c.*(c-1)) .* (c.* x2 - x.^2);
% 
% infor   = envihdrread('i:/HyPlantData/LCC and LAI/20120823-SEL-1350-0600-L2-DUAL_radiance_03aero_1h2o_FINAL_REFL_rect_resampled_HyMap_40_bands_LCC_wrap.hdr');
% LCC     = envidataread('i:/HyPlantData/LCC and LAI/20120823-SEL-1350-0600-L2-DUAL_radiance_03aero_1h2o_FINAL_REFL_rect_resampled_HyMap_40_bands_LCC_wrap',infor);
% infor   = envihdrread('i:/HyPlantData/LCC and LAI/20120823-SEL-1350-0600-L2-DUAL_radiance_03aero_1h2o_FINAL_REFL_rect_resampled_HyMap_40_bands_LAI_wrap.hdr');
% LAI     = envidataread('i:/HyPlantData/LCC and LAI/20120823-SEL-1350-0600-L2-DUAL_radiance_03aero_1h2o_FINAL_REFL_rect_resampled_HyMap_40_bands_LAI_wrap',infor);
% 
% Cab     = LCC(:,:,1);
% LAI     = LAI(:,:,1);
% 
% Cab     = Cab(1:6221,92:end);
% LAI     = LAI(1:6221,92:end);
% % 
% % % 
% % % GR = Fm(2100:2300,170:250);
% % % SO = Fm(2600:2800,80:130);
% 
% 
% SO = Fm(3800:4000,295:335);
% GR = Fm(1830:1930,90:110);
% FO = Fm(250:650,340:440);
% SB = Fm(5150:5300,50:160);
% CO = Fm(4580:4720,50:70);
% PO  = Fm(4930:5100,70:110);
% 
% %GRs0 = F0(2100:2300,170:250);
% SOs0 = F0(3800:4000,295:335);
% GRs0 = F0(1830:1930,90:110);
% FOs0 = F0(250:650,340:440);
% SBs0 = F0(5150:5300,50:160);
% COs0 = F0(4580:4720,50:70);
% POs0 = F0(4930:5100,70:110);
% 
% GRs = GRs0;
% SOs = SOs0;
% FOs = FOs0;%*0.92;
% SBs = SBs0*1.4;
% COs = COs0*1.05;
% POs = POs0;
% 
% meanC = [mean(SO(:)), mean(GR(:)),mean(FO(:)),mean(SB(:)),mean(CO(:)),mean(PO(:))]';
% meanCs0 = [nanmean(SOs0(:)), nanmean(GRs0(:)),nanmean(FOs0(:)),nanmean(SBs0(:)),nanmean(COs0(:)),nanmean(POs0(:))]';
% meanCs = [nanmean(SOs(:)), nanmean(GRs(:)),nanmean(FOs(:)),nanmean(SBs(:)),nanmean(COs(:)),nanmean(POs(:))]';
% stdC = [std(SO(:)), std(GR(:)),std(FO(:)),std(SB(:)),std(CO(:)),std(PO(:))]';
% stdCs = [nanstd(SOs(:)), nanstd(GRs(:)),nanstd(FOs(:)),nanstd(SBs(:)),nanstd(COs(:)),nanstd(POs(:))]';
% 
% R = Fm./F0;
% R(R>2) = NaN;
% R(R<-.5) = NaN;
% 
% figure(1), clf
% labels = {'unclassified','trees','soil','grass','sugar beet','maize','potato','urban'}';
% 
% subplot(144), cd('ojwoodford-sc-b84c84b'), sc(LC,'jet', [0.5 0.5 0.5]), cd .., z=colorbar; set(z,'yticklabel',labels), title('Land Cover')
% subplot(142), cd('ojwoodford-sc-b84c84b'), sc(Fm,[0 2.3],'jet', [0.5 0.5 0.5]), cd ..,  colorbar, title('SIF_{HyPlant}')
% subplot(141), cd('ojwoodford-sc-b84c84b'), sc(F0,'jet', [0.5 0.5 0.5]), cd ..,  colorbar, title('SIF_{SCOPE,0}')
% subplot(143), cd('ojwoodford-sc-b84c84b'),  sc(R,[-.5 1.3], 'jet', [0.5 0.5 0.5]); cd .., colorbar, title('SIF_{HyPlant}/SIF_{SCOPE,0}')
% 
% %caxis([.4 .7]), colorbar
% figure(2), clf
% crops = {'soil','grass','forest','sugar beet','maize','potato'};
% 
% %subplot(211)
% plot(meanC,meanCs0,'ko'), hold on
% plot([0,2.5],[0,2.5],'k')
% for i = 1:6
%     plot([meanC(i),meanC(i)],[meanCs0(i)-stdCs(i),meanCs0(i)+stdCs(i)])
%     plot([meanC(i)-stdC(i),meanC(i)+stdC(i)],[meanCs0(i),meanCs0(i)])
%     text(meanC(i),meanCs0(i)-.1,crops(i))
% end
% set(gca,'xlim',[-1 3.5],'ylim',[-1,3.5])
% ylabel('SCOPE SIF (mW m^{-2} nm^{-2} sr^{-1}')
% title('without physiological response')
% xlabel('HyPlant SIF (mW m^{-2} nm^{-2} sr^{-1}')

%%

[Cabi,epsmax,alpha]   = Cab_LAI_scaling(0,740);
%%
Rin = 655.8605;  % from EC data spreadsheet

J = find(~isnan(Cab.*LAI.*Fm) & Fm>0);
epsmaxi  = epsmax(round(Cab(J))+1);
alphai   = alpha(round(Cab(J))+1);

Fmf     = Fm(J)/Rin;
fLAI    = 1-exp(-alphai.*LAI(J));

Fmn     = Fmf./epsmaxi;

y = zeros(20,1);
figure(6), clf
for k = 1:20
    I = find(fLAI>((k-1)/20) & fLAI<(k/20));
    y(k) = median(Fmn(I));
end
plot((0:.05:.95),y)
%plot(fLAI,Fm(~isnan(Cab.*LAI.*Fm))./epsmaxi,'.')

