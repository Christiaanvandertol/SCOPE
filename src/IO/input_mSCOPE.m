function mly=input_mSCOPE(parameter_file)
    mlayer        = xlsread(parameter_file);
    mly.nly         = mlayer(1,1);
    mly.pLAI        = mlayer(2,:);
    mly.pCab        = mlayer(3,:);
    mly.pCca        = mlayer(4,:);
    mly.pCdm        = mlayer(5,:);
    mly.pCw         = mlayer(6,:);
    mly.pCs         = mlayer(7,:);
    mly.pN          = mlayer(8,:);
    totLAI          = sum(mly.pLAI);
    mly.totLAI      = totLAI;
%     else
%         vi = ones(length(V),1);
%         [~,leafbio,canopy,~,~,~]  = select_input(V,vi,canopy,options,constants);
%         mly.pLAI = canopy.LAI;
%         mly.totLAI = canopy.LAI; 
%         mly.mCab = leafbio.Cab;
%         mly.mCca = leafbio.Cca;
%         mly.mCdm = leafbio.Cdm;
%         mly.mCw = leafbio.Cw;
%         mly.mCs = leafbio.Cs;
%         mly.mN = leafbio.N;
%     end

%% just in case averages
%     mly.mCab        = (mly.pCab*mly.pLAI')./totLAI;
%     mly.mCca        = (mly.pCca*mly.pLAI')./totLAI;
%     mly.mCdm        = (mly.pCdm*mly.pLAI')./totLAI;
%     mly.mCw         = (mly.pCw*mly.pLAI')./totLAI;
%     mly.mCs         = (mly.pCs*mly.pLAI')./totLAI;
%     mly.mN          = (mly.pN*mly.pLAI')./totLAI;
end