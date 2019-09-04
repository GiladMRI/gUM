HankelTemporalLen=2;
WhichInnIdxs=1:nEchos;
[~, ~, ~,H]=ghankel(numel(WhichInnIdxs),HankelTemporalLen,Sz);
%%
[ Ux_LLR, sx_LLR, Vx_LLR ] = batch_svd(H*X);
R1r=Vx_LLR(:,:,2,1)./Vx_LLR(:,:,1,1); % R1 is simply the decay
UpdatedT2SMap_msr=-TimeBetweenEchos_ms./log(abs(R1r));
UpdatedB0Mapr=-(angle(R1r)/(2*pi))/(TimeBetweenEchos_ms/1e3); % in Hz
figure;subplot(1,2,1);gmontage(UpdatedB0Mapr,[-100 100]);title('full data B_0');removeTicks;colorbar;subplot(1,2,2);gmontage(UpdatedT2SMap_msr,[0 100])
%%
% [ U_LLR, s_LLR, V_LLR ] = batch_svd(H*X);
[ U_LLR, s_LLR, V_LLR ] = batch_svd(H*EstXPadded);

R1=V_LLR(:,:,2,1)./V_LLR(:,:,1,1); % R1 is simply the decay

UpdatedT2SMap_ms=-TimeBetweenEchos_ms./log(abs(R1));
UpdatedB0Map=-(angle(R1)/(2*pi))/(TimeBetweenEchos_ms/1e3); % in Hz

% figure;subplot(1,2,1);gmontage(UpdatedB0Map,[-100 100]);title('calib data B_0');removeTicks;colorbar;subplot(1,2,2);gmontage(UpdatedT2SMap_ms,[0 100])

% grmss(sum(V_LLR(:,:,1,:).*V_LLR(:,:,2,:),4))

VH_LLR=(permute43(V_LLR));
GoodDirV=VH_LLR(:,:,1,:);
BadDirV=VH_LLR(:,:,2,:);
% BadDirVb=cat(4,GoodDirV(:,:,1,2),-GoodDirV(:,:,1,1));
BadDirVb=cat(4,GoodDirV(:,:,1,1),-GoodDirV(:,:,1,2));

disp('Found initial BadDir');
%%
figure;subplot(1,2,1);
gmontage(UpdatedB0Map,[-100 100]);
title('calib B_0');removeTicks;colorbar
subplot(1,2,2);
gmontage(UpdatedT2SMap_ms,[0 100])
%% Test TH annihilation
TH=H*X;
% GoodDir=sum(TH.*VH_LLR(:,:,1,:),4);
% BadDir=sum(TH.*VH_LLR(:,:,2,:),4);
GoodDir=sum(TH.*GoodDirV,4);
BadDir=sum(TH.*BadDirV,4);
BadDirb=sum(TH.*BadDirVb,4);

BothDir=cat(4,GoodDir,BadDir,BadDirb);
% ShowAbsAngle(BothDir)

ShowAbsAngle(BothDir)
% %% Test TH annihilation, with U
% TH=H*X;
% UTH=U_LLR(:,:,:,:).*TH;figure;subplot(1,2,1);gmontage(UpdatedB0Map,[-100 100]);title('Innershot temporal hankel based initial B_0');removeTicks;colorbar;subplot(1,2,2);gmontage(UpdatedT2SMap_ms,[0 100])
% GoodDir=sum(UTH.*conj(V_LLR(:,:,1,:)),4);
% BadDir=sum(UTH.*conj(V_LLR(:,:,2,:)),4);
% BothDir=cat(4,GoodDir,BadDir);
% ShowAbsAngle(BothDir)
%% Solve given Bad dir
THDir_FN='THDir.txt';

ImSz=[Sz nEchos];
ImSzF=FillOnesTo16(ImSz);

Main_Ops={'fftc 3','fmac 0 0'};
THDir_Ops={'hankel 2 3 2','fmac 1 8'};

WriteLinopToFile(THDir_FN,{Main_Ops,THDir_Ops});

% test
% BART_TH=bart(['linopScript -L 1 ' THDir_FN],ImSzF,X,Msk,BadDirVb);

% Best: no W, BadDirLambda 1e5, U 3
BadDirLambda=1e5;
WLambda=1e-6;
% THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  ' ' THDir_FN],ImSzF,Sig,Msk,BadDirV);
% THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  '  -R W:3:0:' num2str(WLambda) ' ' THDir_FN],ImSzF,Sig,Msk,BadDirV);
Uval=3;
THA_Rec=bart(['picsS -m -u ' num2str(Uval) ' -R 2:' num2str(BadDirLambda)  ' ' THDir_FN],ImSzF,Sig,Msk,BadDirV);
% THA_Rec=bart(['picsS -m -u ' num2str(Uval) ' -R 2:' num2str(BadDirLambda)  '  -R W:3:0:' num2str(WLambda) ' ' THDir_FN],ImSzF,Sig,Msk,BadDirV);

NTHA_Rec=THA_Rec.*grmss(X)./grmss(THA_Rec);
NTHA_RecC{1}=NTHA_Rec;
THAs=ssimBySlices(NTHA_Rec,X);
fgmontage(NTHA_Rec);ylabel(BadDirLambda);xlabel(num2str([WLambda Uval]));title(mean(THAs));
%% Warm start?
Fac1=grmss(X)./grmss(THA_Rec);
WarmSartFN='WarmStart';
writecfl(WarmSartFN,X/Fac1);

WS_THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  '  -R W:3:0:' num2str(WLambda) ' -W ' WarmSartFN ' ' THDir_FN],ImSzF,Sig,Msk,BadDirVb);

NWS_THA_Rec=WS_THA_Rec.*grmss(X)./grmss(WS_THA_Rec);
WS_THAs=ssimBySlices(NWS_THA_Rec,X);
fgmontage(NWS_THA_Rec);ylabel(BadDirLambda);xlabel(WLambda);title(mean(THAs));
%% Iterate (get B,D using SVD)
for iter=24
    disp(['------ ' num2str(iter) ' -----------']);
    [ Uc_LLR, sc_LLR, Vc_LLR ] = batch_svd(H*THA_Rec);
    
    R1=Vc_LLR(:,:,2,1)./Vc_LLR(:,:,1,1); % R1 is simply the decay
    UpdatedT2SMap_ms=-TimeBetweenEchos_ms./log(abs(R1));
    UpdatedB0Map=-(angle(R1)/(2*pi))/(TimeBetweenEchos_ms/1e3); % in Hz
    
%     UpdatedT2SMap_ms=UpdatedT2SMap_msr;
%     UpdatedB0Map=UpdatedB0Mapr;
%     figure;subplot(1,2,1);gmontage(UpdatedB0Map,[-100 100]);title('Innershot temporal hankel based initial B_0');removeTicks;colorbar;subplot(1,2,2);gmontage(UpdatedT2SMap_ms,[0 100])
    
    UpdatedT2SMap_ms=max(UpdatedT2SMap_ms,4);
    R1M=exp(-TimeBetweenEchos_ms./UpdatedT2SMap_ms);
    R1P=-UpdatedB0Map*2*pi*TimeBetweenEchos_ms/1e3;
    R1x=R1M.*exp(1i*R1P);
    Vc_LLRx=cat(4,ones(Sz),R1x);
    Vc_LLRx=Vc_LLRx./grmss(Vc_LLRx,4)/sqrt(2); % like GoodDirV up to unimportant phase coefficient
    BadDirVb=cat(4,Vc_LLRx(:,:,1,1),-Vc_LLRx(:,:,1,2));

    VH_LLR=(permute43(Vc_LLR));
    BadDirV=VH_LLR(:,:,2,:);
    GoodDirV=VH_LLR(:,:,1,:);
%     BadDirVb=cat(4,GoodDirV(:,:,1,1),-GoodDirV(:,:,1,2));

%     THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  ' ' THDir_FN],ImSzF,Sig,Msk,BadDirVb);
%     THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  '  -R W:3:0:' num2str(WLambda) ' ' THDir_FN],ImSzF,Sig,Msk,BadDirV);
    THA_Rec=bart(['picsS -m -R 2:' num2str(BadDirLambda)  '  -R W:3:0:' num2str(WLambda) ' ' THDir_FN],ImSzF,Sig,Msk,BadDirVb);
    NTHA_Rec=THA_Rec.*grmss(X)./grmss(THA_Rec);
    NTHA_RecC{iter}=NTHA_Rec;
    THAs(iter,:)=ssimBySlices(NTHA_Rec,X);
end
disp('Finished iterations');
%%
figure;plot(THAs)