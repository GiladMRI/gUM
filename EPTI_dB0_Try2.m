EPTIRecP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Recon/Human/';
EBaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
R=load([EBaseP 'meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib_FirstRep_Raw.mat']);

OutP=[EBaseP 'meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns/'];

parameters=load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Data_acq/meas_prot_Human_SMS1_1p9_3SHOT_GE.mat');
parameters=parameters.parameters;

ES_ms=R.meas.prot.iEffectiveEpiEchoSpacing/1000;

FirstTE_ms=9;

TEs_ms=(FirstTE_ms+(0:(nEchos-1))*ES_ms);
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms3=TEs_ms3;
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
B0MapToUse_Hz=SRefB0r;
Thickness_mm=3;
DistBetweenSlices_mm=3;
%%
dB0dx=symD(B0MapToUse_Hz,1);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dz=symD(B0MapToUse_Hz,3)*Thickness_mm/DistBetweenSlices_mm;
% dB0dz=symD(B0MapToUse_Hz,3);

ThroughPlanePhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dz)/1000;
ThroughPlaneDecayS=sinc(ThroughPlanePhaseDiffS/(2*pi));
%%
PEPhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dx)/1000;

PEPhaseDiffSCurSli=perm73(PEPhaseDiffS(:,:,:,6));
%%
PhaseDiffDueToImaging= 2*pi*(kLoc-1-60)/120;
PhaseDiffDueToImaging78=permute(PhaseDiffDueToImaging,[3:8 1 2]);

PEPhaseDiffTotal=PEPhaseDiffSCurSli+PhaseDiffDueToImaging78;
PEPhaseDiffTotalS=sinc(PEPhaseDiffTotal/(2*pi));
%%
ITS_TSC_Cmnds={'fmac 2 0',['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_TSC_Cmnds);
BARTS_Aopx.Others{3}=perm73(ThroughPlaneDecay);
BARTS_Aopx.Others{1}=SensCSMap(:,:,:,1:nccToUse);
BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
%%
nShots=3;
ShotsLocs={1:40,41:80,81:120};
% sig kRO kPE 1 CC 1 1 Echos
isig=ifft1cg(sig,1);
nccToUse=15;
isigP=permute(isig(:,:,1,1:nccToUse,1,1,:),[2 1 3:7]); % kPE RO 1 CC 1 1 Echos
isigPS=zeros([Sz 1 nccToUse 1 1 nEchos nShots]);
for i=1:nShots
    isigPS(ShotsLocs{i},:,:,:,:,:,:,i)=isigP(ShotsLocs{i},:,:,:,:,:,:);
end
MaskSamplesS=abs(isigPS)>0;
FMaskSamplesS=fft1cg(single(MaskSamplesS),1);
sisigPS=sum(isigPS,1);

for i=1:nEchos
    for s=1:nShots
        kLoc(i,s)=find(MaskSamplesS(:,1,1,1,1,1,i,s));
    end
end
SensP=permute(SensCSMap(:,:,:,1:nccToUse),[2 1 3 4]);
SensFMaskSamplesS=SensP.*FMaskSamplesS;

% iisig=ifft1cg(isigP,1);
% QQ=sum(sum(iisig(:,:,1,1:nccToUse,1,1,:),7).*conj(SensP),4);

Ops={'fmac 0 1'};
FMScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots.txt';
WriteLinopToFile(FMScriptFN,Ops);

ImSz16=FillOnesTo16([Sz 1 1 1 1 nEchos]);

BaseFP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
SigFN=[BaseFP 'sisigPS'];
writecfl(SigFN,sisigPS);

SensFMaskSamplesS_FN=[BaseFP 'SensFMaskSamplesS'];
writecfl(SensFMaskSamplesS_FN,SensFMaskSamplesS);

AHsig=bart(['linopScript -d 5 -A ' FMScriptFN],ImSz16,SigFN,SensFMaskSamplesS_FN);

EstLLR=permute(S6.im_recon,[2 1 4 5 6 7 3]);
AHA_EstLLR=bart(['linopScript -d 5 -N ' FMScriptFN],ImSz16,EstLLR,SensFMaskSamplesS_FN);
%% Subspace
WhichEchosToUse=5:55;
[PDBase_LLR, UpdatedB0Map_Hz_LLR, UpdatedT2SMap_ms_LLR, s_vals_LLR, Fitted0_LLR, PDBase0_LLR]=...
    FitToModel_MPBD1CSf(EstLLR,WhichEchosToUse,ES_ms,FirstTE_ms);


T2svalues_ms=linspace(5,300,200);
Decays=exp(-TEs_ms.'./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

nComponentsToUse=2;
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);

Ops_Subspace=['fmac 1 32',Ops];
FMScriptFN_Subspace='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots_Subspace.txt';
WriteLinopToFile(FMScriptFN_Subspace,Ops_Subspace);

CompsPFN=[BaseFP 'CompsP'];
TSCPSensFMaskSamplesS_FN=[BaseFP 'TSCPSensFMaskSamplesS'];

writecfl(CompsPFN,CompsP);

ImSz16Comps=FillOnesTo16([Sz 1 1 1 nComponentsToUse 1]);

% TSCP=exp(-1i*2*pi*SB0x.'.*TEs_ms3/1000);
% TSCP=exp(-1i*2*pi*UpdatedB0Map_Hz_LLR.*TEs_ms3/1000);
TSCP=exp(-1i*2*pi*SRefB0r(:,:,6).*TEs_ms3/1000);

TSCP7=perm73(TSCP);
TSCPSensFMaskSamplesS=TSCP7.*SensFMaskSamplesS;

writecfl(TSCPSensFMaskSamplesS_FN,TSCPSensFMaskSamplesS);
%%
RecSubspace=bart(['picsS -d 5  ' FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesS_FN,CompsPFN);
RecSubspaceX=sum(RecSubspace.*CompsP,6).*TSCP7;
%%
TSCPSensFMaskSamplesS=TSCP7.*SensFMaskSamplesS.*perm73(ThroughPlaneDecayS(:,:,:,6));
writecfl(TSCPSensFMaskSamplesS_FN,TSCPSensFMaskSamplesS);

RecSubspaceTD=bart(['picsS -d 5  ' FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesS_FN,CompsPFN);
RecSubspaceTDX=sum(RecSubspaceTD.*CompsP,6).*TSCP7;

[PDBase_SubTD, UpdatedB0Map_Hz_SubTD, UpdatedT2SMap_ms_SubTD, s_vals_SubTD, Fitted0_SubTD, PDBase0_SubTD]=...
    FitToModel_MPBD1CSf(squeeze(RecSubspaceTDX),WhichEchosToUse,ES_ms,FirstTE_ms);
%%
SmallVarT=50e6;
getSmallVars;
save([EBaseP 'EPTI_dB0_SmallVars.mat'],SmallVars{:});
%%
TSCPSensFMaskSamplesS=TSCP7.*SensFMaskSamplesS.*perm73(ThroughPlaneDecayS(:,:,:,6)).*PEPhaseDiffTotalS;
writecfl(TSCPSensFMaskSamplesS_FN,TSCPSensFMaskSamplesS);

RecSubspaceTDPE=bart(['picsS -d 5  ' FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesS_FN,CompsPFN);
RecSubspaceTDPEX=sum(RecSubspaceTDPE.*CompsP,6).*TSCP7;

[PDBase_SubTDPE, UpdatedB0Map_Hz_SubTDPE, UpdatedT2SMap_ms_SubTDPE, s_vals_SubTDPE, Fitted0_SubTDPE, PDBase0_SubTDPE]=...
    FitToModel_MPBD1CSf(squeeze(RecSubspaceTDPEX),WhichEchosToUse,ES_ms,FirstTE_ms);
%%
NPD0Ref=RefPDBase0r(:,:,6)./grmss(RefPDBase0r(:,:,6)).*grmss(PDBase0_LLR);
PD0All=cat(3,PDBase0_LLR,NPD0Ref,PDBase0_SubTD,PDBase0_SubTDPE);

exIdx=65;
NPD0Ref=RefPDBase0r(:,:,6)./grmss(RefPDBase0r(:,:,6)).*grmss(Fitted0_LLR(:,:,exIdx));
Eex=cat(3,Fitted0_LLR(:,:,exIdx),NPD0Ref,Fitted0_SubTD(:,:,exIdx),Fitted0_SubTDPE(:,:,exIdx));
fgmontagex(Eex);caxis(caxis/1.8);
