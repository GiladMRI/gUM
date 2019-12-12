BaseFP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';

EPTIRecP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Recon/Human/';
EBaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
R=load([EBaseP 'meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib_FirstRep_Raw.mat']);

OutP=[EBaseP 'meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns/'];

parameters=load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Data_acq/meas_prot_Human_SMS1_1p9_3SHOT_GE.mat');
parameters=parameters.parameters;

ES_ms=R.meas.prot.iEffectiveEpiEchoSpacing/1000;

FirstTE_ms=9;

sig=readcfl([OutP 'Sli' num2str(1) '_Sig']);

nEchos=size(sig,7);
Sz=gsize(sig,1:2);

TEs_ms=(FirstTE_ms+(0:(nEchos-1))*ES_ms);
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms3=TEs_ms3;
disp('ok');
%%
nccToUse=15;

isig=ifft1cg(sig,1);
isigP=permute(isig(:,:,1,1:nccToUse,1,1,:),[2 1 3:7]); % kPE RO 1 CC 1 1 Echos
isigPS=zeros([Sz 1 nccToUse 1 1 nEchos nShots]);
for i=1:nShots
    isigPS(ShotsLocs{i},:,:,:,:,:,:,i)=isigP(ShotsLocs{i},:,:,:,:,:,:);
end
    
nShots=3;
ShotsLocs={1:40,41:80,81:120};

MaskSamplesS=abs(isigPS)>0;
FMaskSamplesS=fft1cg(single(MaskSamplesS),1);

for i=1:nEchos
    for s=1:nShots
        kLoc(i,s)=find(MaskSamplesS(:,1,1,1,1,1,i,s));
    end
end
disp('got kLoc');
%%
T2svalues_ms=linspace(5,300,200);
Decays=exp(-TEs_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

Ops={'fmac 0 1'};
FMScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots.txt';
WriteLinopToFile(FMScriptFN,Ops);

Ops_Subspace=['fmac 1 32',Ops];
FMScriptFN_Subspace='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots_Subspace.txt';
WriteLinopToFile(FMScriptFN_Subspace,Ops_Subspace);

ImSz16=FillOnesTo16([Sz 1 1 1 1 nEchos]);

CompsPFN=[BaseFP 'CompsP'];

nComponentsToUse=4;
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);

writecfl(CompsPFN,CompsP);

ImSz16Comps=FillOnesTo16([Sz 1 1 1 nComponentsToUse 1]);
disp('Saved comps');
%%
GRE_Ref_FN='meas_MID04675_FID21637_gre_te4_9';
RefLocs=load([EBaseP GRE_Ref_FN '/Locs.mat']);
RefLocs=RefLocs.RotatedLocs;

sTwixX=load([EBaseP GRE_Ref_FN '/sTwixX.mat']);
sTwixX=sTwixX.sTwixX;

sTwixX.hdr.MeasYaps.sSliceArray.asSlice
sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition % 220 x 213
% Here: 228 x 228
RefS=load([EBaseP GRE_Ref_FN '/B0T2S.mat']);
FOVRef=[220 220];
FOVHere=[228 228];

RefSvals1=gflip(rot90(padarray(imresize(padarray(squeeze(RefS.s_vals(:,:,1,:)),[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefSvals1=circshift(circshift(RefSvals1,2,1),4,2); % positive is downware
RefSvals1r=permute(imresize(RefSvals1,Sz),[2 1 3]);
fgmontagex(RefSvals1r(:,:,6))

Msk=imfillholesBySlices(RefSvals1r>0.0001);

RefPDBase0=gflip(rot90(padarray(imresize(padarray(squeeze(RefS.PDBase0),[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefPDBase0=circshift(circshift(RefPDBase0,2,1),4,2); % positive is downware
RefPDBase0r=permute(imresize(RefPDBase0,Sz),[2 1 3]);

RefB0=gflip(rot90(padarray(imresize(padarray(squeeze(RefS.B0M_Hz),[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefB0=circshift(circshift(RefB0,2,1),4,2); % positive is downware
RefB0r=permute(imresize(RefB0,Sz),[2 1 3]);

SSigs=[0.0001 0.1:0.1:29];
for i=1:size(RefB0r,3)
    disp(i);
    SRefB0r(:,:,i)=SmoothByW(RefB0r(:,:,i),RefSvals1r(:,:,i),[21 21],SSigs);
end
save([OutP 'SRefB0r.mat'],'SRefB0r');

nSlices=size(SRefB0r,3);

B0MapToUse_Hz=SRefB0r;
Thickness_mm=3;
DistBetweenSlices_mm=3;
%%
dB0dx=symD(B0MapToUse_Hz,1);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dz=symD(B0MapToUse_Hz,3)*Thickness_mm/DistBetweenSlices_mm;
% dB0dz=symD(B0MapToUse_Hz,3);

ThroughPlanePhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dz)/1000;
ThroughPlaneDecayS=sinc(ThroughPlanePhaseDiffS/(2*pi));
%%
PEPhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dx)/1000;
%
PhaseDiffDueToImaging= 2*pi*(kLoc-1-60)/120;
PhaseDiffDueToImaging78=permute(PhaseDiffDueToImaging,[3:8 1 2]);

PEPhaseDiffTotal=perm73(PEPhaseDiffS)+PhaseDiffDueToImaging78;
PEPhaseDiffTotalS=sinc(PEPhaseDiffTotal/(2*pi));

ROPhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dy)/1000;
RODecayS=sinc(ROPhaseDiffS/(2*pi));
disp('ok');
%%
PEPhaseDiffTotalS2=PEPhaseDiffTotalS*0;
for SliI=1:nSlices
    disp(SliI);
    PEOversampleFac=31;
    
    tmp=B0MapToUse_Hz(:,:,SliI);
    tmpr=imresize(tmp,[Sz(1)*PEOversampleFac Sz(2)]);
    tmpr2=kron(tmp,ones([PEOversampleFac 1]));
    
    dtmpr=tmpr-tmpr2;
    
    PhaseDiff=2*pi.*EchoTimes_ms3.*dtmpr/1000;
    
    PhaseDiffDueToImaging= 2*pi*(kLoc-1-60)/120;
    
    Patchbase=linspace(-.5,.5,31).';
    PhaseDueToImaging78x=repmat(Patchbase.*PhaseDiffDueToImaging78,Sz);
    
    PEPhaseTotalx=PhaseDueToImaging78x+perm73(PhaseDiff);
    PEPhaseTotalx2=PartitionDim(PEPhaseTotalx,1,Sz(1));
    PEPhaseDiffTotal=permute(sum(cos(PEPhaseTotalx2),1),[9 2:8 1])/PEOversampleFac;
    
    PEPhaseDiffTotalS2(:,:,:,SliI,:,:,:,:)=PEPhaseDiffTotal;
end
disp('ok');
%%
SliI=16;

for SliI=[1:8 10:32]
    disp(['Slice ' num2str(SliI) ' ' datestr(now)]);
    
    SigFN=[BaseFP 'sisigPS' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesS_FN=[BaseFP 'TSCPSensFMaskSamplesS' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTD_FN=[BaseFP 'TSCPSensFMaskSamplesSTD' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTDPE_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPE' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTDPERO_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPERO' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTDPE2RO_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPE2RO' '_Sli' num2str(SliI)];
    
    %     TSCPSensFMaskSamplesSTDPE=readcfl(TSCPSensFMaskSamplesSTDPE_FN);
    %     TSCPSensFMaskSamplesSTDPERO=TSCPSensFMaskSamplesSTDPE.*perm73(RODecayS(:,:,:,SliI));
    %     writecfl(TSCPSensFMaskSamplesSTDPERO_FN,TSCPSensFMaskSamplesSTDPERO);
    %     
    %     disp(['Saved for subspace, slice #' num2str(SliI)]);
    % end

    if(~exist(SigFN,'file'))
        sig=readcfl([OutP 'Sli' num2str(SliI) '_Sig']);
        CurSens=readcfl([OutP 'Sli' num2str(SliI) '_SensCC']);

        % sig kRO kPE 1 CC 1 1 Echos
        isig=ifft1cg(sig,1);
        isigP=permute(isig(:,:,1,1:nccToUse,1,1,:),[2 1 3:7]); % kPE RO 1 CC 1 1 Echos
        isigPS=zeros([Sz 1 nccToUse 1 1 nEchos nShots]);
        for i=1:nShots
            isigPS(ShotsLocs{i},:,:,:,:,:,:,i)=isigP(ShotsLocs{i},:,:,:,:,:,:);
        end
        sisigPS=sum(isigPS,1);


        SensP=permute(CurSens(:,:,:,1:nccToUse),[2 1 3 4]);
        SensFMaskSamplesS=SensP.*FMaskSamplesS;

        writecfl(SigFN,sisigPS);

        SensFMaskSamplesS_FN=[BaseFP 'SensFMaskSamplesS'];
        writecfl(SensFMaskSamplesS_FN,SensFMaskSamplesS);
    end
    %
    % AHsig=bart(['linopScript -d 5 -A ' FMScriptFN],ImSz16,SigFN,SensFMaskSamplesS_FN);
    %
    % EstLLR=permute(S6.im_recon,[2 1 4 5 6 7 3]);
    % AHA_EstLLR=bart(['linopScript -d 5 -N ' FMScriptFN],ImSz16,EstLLR,SensFMaskSamplesS_FN);
    % disp('Saved sig and stuff');
    %
    % TSCP=exp(-1i*2*pi*SB0x.'.*TEs_ms3/1000);
    % TSCP=exp(-1i*2*pi*UpdatedB0Map_Hz_LLR.*TEs_ms3/1000);
    TSCP=exp(-1i*2*pi*SRefB0r(:,:,SliI).*TEs_ms3/1000);
    TSCP7=perm73(TSCP);
    TSCPSensFMaskSamplesS=Msk(:,:,SliI).*TSCP7.*SensFMaskSamplesS;
    TSCPSensFMaskSamplesSTD=Msk(:,:,SliI).*TSCP7.*SensFMaskSamplesS.*perm73(ThroughPlaneDecayS(:,:,:,SliI));
    TSCPSensFMaskSamplesSTDPE=Msk(:,:,SliI).*TSCP7.*SensFMaskSamplesS.*perm73(ThroughPlaneDecayS(:,:,:,SliI)).*PEPhaseDiffTotalS(:,:,:,SliI,:,:,:,:);
    
    TSCPSensFMaskSamplesSTDPE2=Msk(:,:,SliI).*TSCP7.*SensFMaskSamplesS.*perm73(ThroughPlaneDecayS(:,:,:,SliI)).*PEPhaseDiffTotalS2(:,:,:,SliI,:,:,:,:);
    
    TSCPSensFMaskSamplesSTDPERO=TSCPSensFMaskSamplesSTDPE.*perm73(RODecayS(:,:,:,SliI));
    TSCPSensFMaskSamplesSTDPE2RO=TSCPSensFMaskSamplesSTDPE2.*perm73(RODecayS(:,:,:,SliI));
    
    writecfl(TSCPSensFMaskSamplesS_FN,TSCPSensFMaskSamplesS);
    writecfl(TSCPSensFMaskSamplesSTD_FN,TSCPSensFMaskSamplesSTD);
    writecfl(TSCPSensFMaskSamplesSTDPE_FN,TSCPSensFMaskSamplesSTDPE);
    writecfl(TSCPSensFMaskSamplesSTDPERO_FN,TSCPSensFMaskSamplesSTDPERO);
    writecfl(TSCPSensFMaskSamplesSTDPE2RO_FN,TSCPSensFMaskSamplesSTDPE2RO);
    
    disp(['Saved for subspace, slice #' num2str(SliI)]);
end
%%
WhichEchosToUse=15:65;

SliI=3;
% for SliI=[1:8 10:32]
for SliI=[1:8 10:32]
    SigFN=[BaseFP 'sisigPS' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesS_FN=[BaseFP 'TSCPSensFMaskSamplesS' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTD_FN=[BaseFP 'TSCPSensFMaskSamplesSTD' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTDPE_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPE' '_Sli' num2str(SliI)];
    
    TSCPSensFMaskSamplesSTDPERO_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPERO' '_Sli' num2str(SliI)];
    TSCPSensFMaskSamplesSTDPE2RO_FN=[BaseFP 'TSCPSensFMaskSamplesSTDPE2RO' '_Sli' num2str(SliI)];
    
    TSCP=exp(-1i*2*pi*SRefB0r(:,:,SliI).*TEs_ms3/1000);
    TSCP7=perm73(TSCP);
    
    RecStr=' -S ';
%     RecStr=' -S -R W:3:32:0.01 ';

    RecSubspace=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesS_FN,CompsPFN);
    RecSubspaceX=sum(RecSubspace.*CompsP,6).*TSCP7;
    
    RecSubspaceTD=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesSTD_FN,CompsPFN);
    RecSubspaceTDX=sum(RecSubspaceTD.*CompsP,6).*TSCP7;
    
    RecSubspaceTDPE=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesSTDPE_FN,CompsPFN);
    RecSubspaceTDPEX=sum(RecSubspaceTDPE.*CompsP,6).*TSCP7;

    RecSubspaceTDPERO=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesSTDPERO_FN,CompsPFN);
    RecSubspaceTDPEROX=sum(RecSubspaceTDPERO.*CompsP,6).*TSCP7;
    
%     RecSubspaceTDPE2RO=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesSTDPE2RO_FN,CompsPFN);
%     RecSubspaceTDPE2ROX=sum(RecSubspaceTDPE2RO.*CompsP,6).*TSCP7;

    % [PDBase_Sub, UpdatedB0Map_Hz_Sub, UpdatedT2SMap_ms_Sub, s_vals_Sub, Fitted0_Sub, PDBase0_Sub]=...
    %     FitToModel_MPBD1CSf(squeeze(RecSubspaceX),WhichEchosToUse,ES_ms,FirstTE_ms);
    % 
    % [PDBase_SubTD, UpdatedB0Map_Hz_SubTD, UpdatedT2SMap_ms_SubTD, s_vals_SubTD, Fitted0_SubTD, PDBase0_SubTD]=...
    %     FitToModel_MPBD1CSf(squeeze(RecSubspaceTDX),WhichEchosToUse,ES_ms,FirstTE_ms);
    % 
    % [PDBase_SubTDPE, UpdatedB0Map_Hz_SubTDPE, UpdatedT2SMap_ms_SubTDPE, s_vals_SubTDPE, Fitted0_SubTDPE, PDBase0_SubTDPE]=...
    %     FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEX),WhichEchosToUse,ES_ms,FirstTE_ms);

    [PDBase_SubTDPERO, UpdatedB0Map_Hz_SubTDPERO, UpdatedT2SMap_ms_SubTDPERO, s_vals_SubTDPERO, Fitted0_SubTDPERO, PDBase0_SubTDPERO]=...
        FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEROX),WhichEchosToUse,ES_ms,FirstTE_ms);
    
%     [PDBase_SubTDPE2RO, UpdatedB0Map_Hz_SubTDPE2RO, UpdatedT2SMap_ms_SubTDPE2RO, s_vals_SubTDPE2RO, Fitted0_SubTDPE2RO, PDBase0_SubTDPE2RO]=...
%         FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPE2ROX),WhichEchosToUse,ES_ms,FirstTE_ms);

    disp(['Rec and fitted slice ' num2str(SliI)]);

    % RecSubspaceS(:,:,:,SliI)=squeeze(RecSubspaceX);
    % RecSubspaceTDS(:,:,:,SliI)=squeeze(RecSubspaceTDX);
    % RecSubspaceTDPES(:,:,:,SliI)=squeeze(RecSubspaceTDPEX);
    RecSubspaceTDPEROS(:,:,:,SliI)=squeeze(RecSubspaceTDPEROX);
    RecSubspaceTDPE2ROS(:,:,:,SliI)=squeeze(RecSubspaceTDPE2ROX);

    % Fitted0_SubS(:,:,:,SliI)=Fitted0_Sub;
    % Fitted0_SubTDS(:,:,:,SliI)=Fitted0_SubTD;
    % Fitted0_SubTDPES(:,:,:,SliI)=Fitted0_SubTDPE;
    Fitted0_SubTDPEROS(:,:,:,SliI)=Fitted0_SubTDPERO;
    Fitted0_SubTDPE2ROS(:,:,:,SliI)=Fitted0_SubTDPE2RO;
end
%%
save([EBaseP 'RecSubspaceTDPEROS.mat'],'RecSubspaceTDPEROS');
save([EBaseP 'RecSubspaceTDPE2ROS.mat'],'RecSubspaceTDPE2ROS');
%%
AllRecs4=cat(5,AllRecs,Fitted0_SubTDPEROS);
save([EBaseP 'AllRecs4.mat'],'AllRecs4');
%%
AllRecsA=cat(5,Fitted0_SubS,Fitted0_SubTDS,Fitted0_SubTDPES);
save([EBaseP 'AllRecs.mat'],'AllRecs');
%%
AllRecs=cat(5,Fitted0_SubS,Fitted0_SubTDS,Fitted0_SubTDPES);

save([EBaseP 'AllRecs.mat'],'AllRecs');
%%
AllRecsF=cat(5,Fitted0_SubS,Fitted0_SubTDS,Fitted0_SubTDPES,Fitted0_SubTDPEROS,Fitted0_SubTDPE2ROS);
save([EBaseP 'AllRecsF.mat'],'AllRecsF');

QQa=abs(perm43(squeeze(abs(AllRecs4(:,:,55,:,[1 2 4]))-abs(AllRecs4(:,:,55,:,3)))));

SliX=8;
fgmontagex(AllRecs4(:,:,55,SliX,:));caxis(caxis/2)
figure;imagesc(QQa(:,:,:,SliX));removeTicks;
fgmontagex(SRefB0r(:,:,SliX));
%%
exIdx=65;
NPD0Ref=RefPDBase0r(:,:,SliI)./grmss(RefPDBase0r(:,:,SliI)).*grmss(Fitted0_Sub(:,:,exIdx));
Eex=cat(3,Fitted0_Sub(:,:,exIdx),NPD0Ref,Fitted0_SubTD(:,:,exIdx),Fitted0_SubTDPE(:,:,exIdx));
fgmontagex(Eex);caxis(caxis/1.8);
%%
RecSubspaceTDPEX=RecSubspaceTDPEX2;
[PDBase_SubTDPE, UpdatedB0Map_Hz_SubTDPE, UpdatedT2SMap_ms_SubTDPE, s_vals_SubTDPE, Fitted0_SubTDPE, PDBase0_SubTDPE]=...
FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEX),WhichEchosToUse,ES_ms,FirstTE_ms);
UpdatedT2SMap_ms_SubTDPE2=UpdatedT2SMap_ms_SubTDPE;Fitted0_SubTDPE2=Fitted0_SubTDPE;
RecSubspaceTDPEX=RecSubspaceTDPEX3;
[PDBase_SubTDPE, UpdatedB0Map_Hz_SubTDPE, UpdatedT2SMap_ms_SubTDPE, s_vals_SubTDPE, Fitted0_SubTDPE, PDBase0_SubTDPE]=...
FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEX),WhichEchosToUse,ES_ms,FirstTE_ms);
UpdatedT2SMap_ms_SubTDPE3=UpdatedT2SMap_ms_SubTDPE;Fitted0_SubTDPE3=Fitted0_SubTDPE;
RecSubspaceTDPEX=RecSubspaceTDPEX4;
[PDBase_SubTDPE, UpdatedB0Map_Hz_SubTDPE, UpdatedT2SMap_ms_SubTDPE, s_vals_SubTDPE, Fitted0_SubTDPE, PDBase0_SubTDPE]=...
FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEX),WhichEchosToUse,ES_ms,FirstTE_ms);
UpdatedT2SMap_ms_SubTDPE4=UpdatedT2SMap_ms_SubTDPE;Fitted0_SubTDPE4=Fitted0_SubTDPE;

RecSubspaceTDPE_Comps=cat(3,RecSubspaceTDPEX2(:,:,exIdx),RecSubspaceTDPEX3(:,:,exIdx),RecSubspaceTDPEX4(:,:,exIdx));;
RecSubspaceTDPE_CompsF=cat(3,Fitted0_SubTDPE2(:,:,exIdx),Fitted0_SubTDPE3(:,:,exIdx),Fitted0_SubTDPE4(:,:,exIdx));

fgmontagex(perm43(cat(4,RecSubspaceTDPE_Comps,RecSubspaceTDPE_CompsF)));caxis(caxis/1.8);

fgmontagex(cat(3,Fitted0_SubTDPE2(:,:,exIdx),Fitted0_SubTDPE3(:,:,exIdx),Fitted0_SubTDPE4(:,:,exIdx)),'Size',[1 3]);
%%
NPD0Ref=RefPDBase0r(:,:,SliI)./grmss(RefPDBase0r(:,:,SliI)).*grmss(PDBase0_Sub);
PD0All=cat(3,PDBase0_Sub,NPD0Ref,PDBase0_SubTD,PDBase0_SubTDPE);
%%
fgmontagex(S16.im_recon(:,:,exIdx).');title('S16');
fgmontagex(RecSubspaceTDPEX(:,:,exIdx));title('TDPE');
fgmontagex(Fitted0_SubTDPE(:,:,exIdx));title('TDPE Fitted');

%%
SmallVarT=50e6;exIdx=65;
NPD0Ref=RefPDBase0r(:,:,SliI)./grmss(RefPDBase0r(:,:,SliI)).*grmss(Fitted0_Sub(:,:,exIdx));
Eex=cat(3,Fitted0_Sub(:,:,exIdx),NPD0Ref,Fitted0_SubTD(:,:,exIdx),Fitted0_SubTDPE(:,:,exIdx));
fgmontagex(Eex);caxis(caxis/1.8);

getSmallVars;
save([EBaseP 'EPTI_dB0_SmallVars.mat'],SmallVars{:});
%%
SpiP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1/';
for s=14 % 1:16
    % s=14;
    for i=1:4
%         Maps{i}=readcfl([SpiP 'gB0xSli' num2str(s) '_RA1_13_25A_Lm4c_Elem' num2str(i-1) '_iter7000']);
        Maps{i}=readcfl([SpiP 'gB0xSli' num2str(s) '_RA1_13_25A_Lm5b_Elem' num2str(i-1) '_iter7000']);
        if(i==4), Maps{i}=1./Maps{i}; end
        ElemL{i}=readcfl([SpiP 'gB0xSli' num2str(s) '_RA1_13_25A_ElemsL_' num2str(i-1)]);
        MElemL(:,:,:,i)=Maps{i}.*squeeze(ElemL{i});
        if(i>1), MElemL(:,:,:,i)=exp(MElemL(:,:,:,i));end
    end
    SpiRec2(:,:,:,s)=prod(MElemL,4);
end
%%