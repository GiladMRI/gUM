MIDStr=BaseFN(6:11);
FN=[ScanP BaseFN '.dat'];
disp('ok');
%%
mainP=[ScanP BaseFN];
mkdir(mainP);

system(['chmod +777 -R ' mainP]);
disp([mainP ' Created']);

setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/bart-0.4.04b')

BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';
nccToUse=31;

Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;
%% Read raw
AData = mapVBVD(FN);
if(iscell(AData))
    ADatax=AData{end};
else
    ADatax=AData;
end
asSlice=ADatax.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[ADatax.hdr.Phoenix.sSliceArray.asSlice{:}];
end

nSlices=numel(ADatax.hdr.Phoenix.sSliceArray.asSlice);
try
    WipMemBlock=ADatax.hdr.MeasYaps.sWiPMemBlock;
catch
    WipMemBlock=ADatax.hdr.MeasYaps.sWipMemBlock;
end
MB=1; %WipMemBlock.alFree{9};
% MB
nSlicesNoMB=nSlices/MB;

% nRepsHdr=1+ADatax.hdr.Meas.lRepetitions;
nRepsHdr=1+ADatax.hdr.MeasYaps.lRepetitions;
disp('Read raw');
%%
TrajType=WipMemBlock.adFree{12};
ResType=floor((TrajType-10)/2)+1; % 2,1.5,1
TimingType=mod(TrajType-10,2)+1; % 5ms, 10ms
load('GAll5ms10ms.mat');
GTrajaCBase=GAll(:,TimingType,ResType);

if(TimingType==1)
    nInnerShots=10;
else
    nInnerShots=5;
end
%% Do something with the noise!
% noise=AData{1}.noise();
%%
for s=1:nSlices
    try
        SlbLoc(1,s)=asSlice(s).sPosition.dSag;
    catch
        SlbLoc(1,s)=0;
    end
    try
        SlbLoc(2,s)=asSlice(s).sPosition.dCor;
    catchSlicesToRead
        SlbLoc(2,s)=0;
    end
    try
        SlbLoc(3,s)=asSlice(s).sPosition.dTra;
    catch
        SlbLoc(3,s)=0;
    end
end

RotMat = transpose(Quat2RotMat(ADatax.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;

[U IA IB]=unique(ADatax.image.slicePos(3,:));
Qoffset=ADatax.image.slicePos(1:3,IA);

Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

try
    FOVx=ADatax.hdr.Meas.ReadFOV;
catch
    FOVx=ADatax.hdr.Config.ReadFoV;
end

FOVx=220;

dFOV=FOVx/1000;

FOV_mm=FOVx;


disp('ok');

GradReduceFac=WipMemBlock.adFree{4};
PhiRotPerRep=WipMemBlock.adFree{8};
AcqDwellTime_us=WipMemBlock.adFree{13}/1000;
%%
% ADCsToRead=1:10;
% ADCsToRead=1;
SlicesToRead=[10];
for SlicesToRead=1:nSlices
RepsToRead=1:(nRepsHdr-1);
% ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,ADCsToRead,:,:,:,:,:,:);
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,:,:,:,:,:,:,:);

ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);


dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;

disp('Read data');
%%
GTrajaC=GTrajaCBase/GradReduceFac;

g=GTrajaC;
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

kK=k*FOV_mm/1000/2/pi;

Kmax=ceil(max(abs(kK)));
NTrg=Kmax*2;
res_mm=FOVx/(max(abs(kK))*2);
%
figure;
plot(kK);
axis square;
axis equal;
ylabel('k');
%%
ChRMS=grmss(ADataIsL,[1 3]);
[SChRMS ChOrd]=sort(ChRMS,'descend');
% Data=squeeze(ADataIsL(:,ChOrd(1),1:nRepsToUse));

Ch2D=CombineDims(ADataIsL,[3 1]);
[~,S,sccmtx] = svd(Ch2D(1:10:end,:),'econ');
clear Ch2D 

ncc=31;
ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
for i=1:ncc
    ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
end
DataC=permute(ADataIsLCC,[1 3 2]);
disp('ok cc');
%%
TrgSz=[NTrg NTrg];
OnesSens=ones(TrgSz);
%%
Extra=4;
kKA=[zeros(4,1); kK; repmat(kK(end),[10 1])]; % in us

% Delay_us=-2.5; % positive data rightward %-6; % clockwise
Delays_us=-2:0.1:4;

Delays_us=1.6;
dlyI=1;
% for dlyI=1:numel(Delays_us)
    disp('-------');
    disp(dlyI);
    Delay_us=Delays_us(dlyI);
    AcqTimePoints_us=Extra*GradDwellTime_us+Delay_us+(0:AcqDwellTime_us:50000);

    Traj=interp1((0:(numel(kKA)-1))*GradDwellTime_us,kKA,AcqTimePoints_us);
    disp('ok delay');
    % plotTrajCrossPoints;
    
    nRepsToUse=nRepsHdr-1;

    TrajM=Traj.*exp(-1i*2*pi*PhiRotPerRep/360*(0:nRepsToUse-1)');

    DataC=permute(ADataIsLCC,[1 3 2]);
%     TrajM=TrajMa*exp(1i*2*pi*0/360)/1.01; % positive counter clockwize

    BARTTraj=cat(3,real(TrajM),imag(TrajM),imag(TrajM)*0);
    BARTTrajP=permute(BARTTraj,[3 2 1]);

    kx=BARTTrajP(1,:,:)*2*pi;
    ky=BARTTrajP(2,:,:)*2*pi;

    modx=exp(-1i.*(dx*kx+dy*ky));
    disp('ok mod');
    %% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
TrajPartToUse=0000+(1:2000);
RepsToUse=1:(nRepsHdr-1);

DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 3 5 6 7 8 2]).*modx;
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 size(DataPC,8)]);

Rec1p=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
% fgmontage(Rec1cp)
fgmontage(grmss(Rec1p,8));removeTicks;
title([num2str(numel(RepsToUse)) ' shots data x ' num2str(numel(TrajPartToUse)) 'points, pre channel recon, RMS']);
%%
Rec1pRMS=grmss(Rec1p,8);
[Out B BN]=CalcSlicesSNR(Rec1pRMS,true,7);
Msk=imfillholesBySlices(~BN);
Msk=getLargestComponent(Msk);
% Msk=imfillholesBySlices(Rec1pRMS>0.001);
se = strel('disk', 3);
DMsk=imdilate(Msk,se,'same');
SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p).*DMsk,0,TrgSz);
SelfSens1=SelfSens(:,:,:,1);
disp('ok SelfSens');
%% Recon each channel separately, CC
DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 ncc]);

disp('ok rec pre channel cc');
%%
SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
disp('ok SensCC');
    %% Per innershot recon (with added B0 variation)
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;

    InnerShotLen=floor(nTrajToUse/nInnerShots);
    EffTimeBetweenInnershots_ms=AcqDwellTime_us*InnerShotLen/1000;

    RepsToUse=1:(nRepsHdr-1);

    DataCCP=permute(DataPC,[1:3 8 4:7]);
    disp('ok x');
    %% Per innershot recon, several together
    clear InnerShotTraj InnerShotData
    for i=1:nInnerShots
        TrajPartToUseC{i}=(i-1)*InnerShotLen+(1:2000);
        InnerShotTraj{i}=BARTTrajP(:,TrajPartToUseC{i},RepsToUse);
        InnerShotData{i}=DataCCP(:,TrajPartToUseC{i},RepsToUse,1:nccToUse);
    end
    InnerShotTraj=cat(8,InnerShotTraj{:});
    InnerShotData=cat(8,InnerShotData{:});
    Rec1ccpMt=bart('pics -m -R W:3:0:0.0000001 -t ',InnerShotTraj,InnerShotData,DMsk.*SensCC(:,:,:,1:nccToUse));
    Rec1ccpMt=squeeze(Rec1ccpMt);
    disp('ok per innershot together');
    %% Spiral-ins
    clear InnerShotTrajIn InnerShotDataIn
    for i=1:size(Rec1ccpMt,3)
        TrajPartToUseCIn{i}=i*InnerShotLen-2000+(1:2000);
        InnerShotTrajIn{i}=BARTTrajP(:,TrajPartToUseCIn{i},RepsToUse);
        InnerShotDataIn{i}=DataCCP(:,TrajPartToUseCIn{i},RepsToUse,1:nccToUse);
    end
    InnerShotTrajIn=cat(8,InnerShotTrajIn{:});
    InnerShotDataIn=cat(8,InnerShotDataIn{:});
    Rec1ccpMtIn=bart('pics -m -R W:3:0:0.0000001 -t ',InnerShotTrajIn,InnerShotDataIn,DMsk.*SensCC(:,:,:,1:nccToUse));
    Rec1ccpMtIn=squeeze(Rec1ccpMtIn);
    disp('ok per innershot together spiral in');
    %% Spiral out-ins
    clear InnerShotTrajOutIn InnerShotDataOutIn
    for i=1:size(Rec1ccpMt,3)-1
        TrajPartToUseCOutIn{i}=i*InnerShotLen+(-1000:1000);
        InnerShotTrajOutIn{i}=BARTTrajP(:,TrajPartToUseCOutIn{i},RepsToUse);
        InnerShotDataOutIn{i}=DataCCP(:,TrajPartToUseCOutIn{i},RepsToUse,1:nccToUse);
    end
    InnerShotTrajOutIn=cat(8,InnerShotTrajOutIn{:});
    InnerShotDataOutIn=cat(8,InnerShotDataOutIn{:});
    Rec1ccpMtOutIn=bart('pics -m -R W:3:0:0.0000001 -t ',InnerShotTrajOutIn,InnerShotDataOutIn,DMsk.*SensCC(:,:,:,1:nccToUse));
    Rec1ccpMtOutIn=squeeze(Rec1ccpMtOutIn);
    disp('ok per innershot together spiral out-in');
    
    ResByDelay{dlyI}={Rec1ccpMt,Rec1ccpMtIn,Rec1ccpMtOutIn};
% end
%%
% for dlyI=1:numel(ResByDelay)
%     ResByDelayM(:,:,:,:,dlyI)=cat(4,ResByDelay{dlyI}{1}(:,:,2:end),ResByDelay{dlyI}{2}(:,:,1:end-1),ResByDelay{dlyI}{3});
% end
% %%
% ResByDelayM=single(ResByDelayM);
% save('ForDelayCalibByOutInRot.mat','ResByDelayM');
%%
% fgmontage(ResByDelayM(:,:,5,1,1:3:end),'Size',[3 7]);MaximizeFig
% fgmontage(ResByDelayM(:,:,5,2,1:3:end),'Size',[3 7]);MaximizeFig
%% Chosen: 1.6us
% fgmontage(ResByDelayM(:,:,5,1,37+(-2:2)),'Size',[1 5]);
% fgmontage(ResByDelayM(:,:,5,2,37+(-2:2)),'Size',[1 5]);
%%
ShowAbsAngle(Rec1ccpMtOutIn);subplot(1,2,2);title('per innershot, spiral-out-ins');
ShowAbsAngle(Rec1ccpMt(:,:,2:end));subplot(1,2,2);title('per innershot, spiral-outs');
ShowAbsAngle(Rec1ccpMtIn(:,:,1:end-1));subplot(1,2,2);title('per innershot, spiral-ins');
%%
OutGifFN=[mainP filesep 'SpiralOut_' num2str(Delay_us) 'us_Sli' num2str(SlicesToRead) '.gif'];
delete(OutGifFN);
figure;
pause(1);
imagesc(abs(Rec1ccpMt(:,:,1)));colormap gray
gif(OutGifFN)
for i=2:size(Rec1ccpMt,3)
    imagesc(abs(Rec1ccpMt(:,:,i)));colormap gray
    gif
end
disp(OutGifFN);
%%
HankelTemporalLen=2;
WhichInnIdxs=1:size(Rec1ccpMt,3);
[~, ~, ~,HInnerShot]=ghankel(numel(WhichInnIdxs),HankelTemporalLen,TrgSz);
%%
[ U_LLR, s_LLR, V_LLR ] = batch_svd(HInnerShot*(Rec1ccpMt(:,:,WhichInnIdxs)));

R1=V_LLR(:,:,2,1)./V_LLR(:,:,1,1); % R1 is simply the decay

InnerShotDiff_ms=InnerShotLen*AcqDwellTime_us/1e3;

UpdatedT2SMap_ms=-InnerShotDiff_ms./log(abs(R1));
UpdatedB0Map=-(angle(R1)/(2*pi))/(InnerShotDiff_ms/1e3); % in Hz

figure;subplot(1,2,1);
gmontage(UpdatedB0Map,[-100 100]);
title('Innershot temporal hankel based initial B_0');removeTicks;colorbar
subplot(1,2,2);
gmontage(UpdatedT2SMap_ms,[0 100]);removeTicks;
% Updated_TSC_L2=exp(-1i*angle(R1).*permute(0:(nTS_L2-1),[1 3 2]));
% % Updated_TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
SimSpiralIn=Rec1ccpMt.*(R1.^0.5);

UpdatedB0Map0=UpdatedB0Map;
UpdatedT2SMap_ms0=UpdatedT2SMap_ms;
%% Test TH annihilation
% TH=HInnerShot*Rec1ccpMt(:,:,WhichInnIdxs);
% VH_LLR=(permute43(V_LLR));
% GoodDir=sum(TH.*VH_LLR(:,:,1,:),4);
% BadDir=sum(TH.*VH_LLR(:,:,2,:),4);
% BothDir=cat(4,GoodDir,BadDir);
% ShowAbsAngle(BothDir)
%% Out+In gif
OutInGifFN=[mainP filesep 'SpiralOutandIn_' num2str(Delay_us) 'us_Sli' num2str(SlicesToRead) '.gif'];
delete(OutInGifFN);
figure;pause(1);
imagesc(cat(2,abs(Rec1ccpMt(:,:,1)),abs(Rec1ccpMtIn(:,:,1))));colormap gray;removeTicks;axis equal
gif(OutInGifFN)
for i=2:size(Rec1ccpMt,3)
    imagesc(cat(2,abs(Rec1ccpMt(:,:,i)),abs(Rec1ccpMtIn(:,:,i))));colormap gray;removeTicks;axis equal
    gif
end
disp(OutInGifFN);
%% out in rotation gif
OutInRotGifFN=[mainP filesep 'SpiralOutInRot_' num2str(Delay_us) 'us_Sli' num2str(SlicesToRead) '.gif'];
delete(OutInRotGifFN);
figure;
pause(1)
imagesc(abs(Rec1ccpMt(:,:,5)));colormap gray;removeTicks;axis equal
% pause(1)
gif(OutInRotGifFN);
imagesc(abs(Rec1ccpMtIn(:,:,5)));colormap gray;removeTicks;axis equal
gif
disp(OutInRotGifFN);
%%
figure;imshowpair(abs(Rec1ccpMt(:,:,3)),abs(Rec1ccpMtIn(:,:,3)),'checkerboard')
%%
nAcqPoints=numel(Traj);

% STraj=TrajQ.';
CurRep=1;
STraj=TrajM(CurRep,:);
STraj3=CTo3Rows(STraj);

% STraj3=BARTTrajP(:,:,CurRep);

Sz=TrgSz;
Sz16=FillOnesTo16(Sz);

SnufftStruct = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
disp('ok 0');
%%
nTS=20;
nTS=50;
TSB=GetTSCoeffsByLinear(nAcqPoints,nTS);

Kerns=NUFFT_to_Toep_2blocks(SnufftStruct,TSB);

% TimePoints_ms=linspace(0,nAcqPoints*1e3/SpBW,nTS);
% TimePoints_ms=(0:(nAcqPoints-1))*AcqDwellTime_us/1000;
TimePoints_ms=linspace(0,AcqTimePoints_us(end)/1000,nTS);
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TSC=exp(1i.*2*pi*UpdatedB0Map*(1e-3).*TimePoints_ms3);
% [x y z Channels 1 TS]

ScriptFN_TS=[BaseSP 'NuftTSC.txt'];
ScriptFN_TS=[BaseSP 'nuftTSC_N.txt'];

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # Data is [1 #Traj 1 Ch]
TSBP=permute(TSB,[3 1 4 5 6 7 2]);
TSCP=permute(TSC,[1 2 7 6 5 4 3]);
KernsP=permute(Kerns,[1 2 7 6 5 4 3]);
disp('ok a');
%% given B0
% nccToUse=7; % 7: 46 sec
% nccToUse=31;
% BigTrajPartToUse=1:45400;
BigTrajPartToUse=1:(floor(50000/AcqDwellTime_us/10)*10);
%% Given B0, T2S, W regularized
% 31Ch: 260sec.
CurRep=1;
R1x=min(max(abs(R1),0.3),1).*exp(-1i*angle(R1));
TSCx=R1x.^(TimePoints_ms3/InnerShotDiff_ms);
TSCxP=permute(TSCx,[1 2 7 6 5 4 3]);

WLambda=0.1;
% RegStr=['-R W:3:0:' num2str(WLambda)];
% tic;RecgD_W=bart(['picsS -m ' RegStr ' ' ScriptFN_TS],Sz16,DataCCP(:,BigTrajPartToUse,CurRep,1:nccToUse),...
%     DMsk.*SensCC(:,:,:,1:nccToUse),STraj3(:,BigTrajPartToUse),TSBP(:,BigTrajPartToUse,:,:,:,:,:,:,:),TSCxP,KernsP);
% t=toc;
% fgmontage(RecgD_W);MaximizeFig;xlabel(nccToUse);ylabel(t);title(['Single shot, given B_0,T_2^* estimate, ' RegStr]);
%%
% Rec1ccpMtC=Rec1ccpMt.*exp(-1i*angle(Rec1ccpMt(:,:,1)));
% Rec1ccpMtCM=Rec1ccpMt./Rec1ccpMt(:,:,1);
% fgmontage(Rec1ccpMt(:,:,1));title('39 shots 2ms, for ref');
%% Per rep/s, given B0,T2*, Short part of traj
% %% preparation
% WhichInnerShot=3;
% nTSShort=3;
% TSBShort=GetTSCoeffsByLinear(numel(TrajPartToUseC{WhichInnerShot}),nTSShort);
% TSBPShort=permute(TSBShort,[3 1 4 5 6 7 2]);
% 
% for CurRep=1:(nRepsHdr-1)
%     disp(CurRep);
%     STraj=TrajM(CurRep,TrajPartToUseC{WhichInnerShot});
%     STraj3MShort(:,:,CurRep)=CTo3Rows(STraj);
%     SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
%     KernsByRepShort{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRep,TSBShort);
%     KernsPShortC{CurRep}=permute(KernsByRepShort{CurRep},[1 2 7 6 5 4 3]);
% end
% KernsPMShort=cat(3,KernsPShortC{:});
% clear KernsPShortC
% 
% TimePointsShort_ms=linspace(0,AcqTimePoints_us(TrajPartToUseC{1}(end))/1000,nTSShort);
% TimePointsShort_ms3=permute(TimePointsShort_ms,[1 3 2]);
% TSCxShort=R1x.^(TimePointsShort_ms3/InnerShotDiff_ms);
% TSCxPShort=permute(TSCxShort,[1 2 7 6 5 4 3]);
% %% given B0,T2*, several reps together
% CurReps=1:39;
% % 14 sec
% 
% RecgD_SReps=bart(['picsS -m ' RegStr ' ' ScriptFN_TS],Sz16,DataCCP(:,TrajPartToUseC{WhichInnerShot},CurReps,1:nccToUse),...
%         SensCC(:,:,:,1:nccToUse).*TightMask,STraj3MShort(:,:,CurReps),TSBPShort,...
%         TSCxPShort,sum(KernsPMShort(:,:,CurReps,:,:,:,:),3));
% 
% fgmontage(RecgD_SReps);title('Multi-shot, short part of traj, given B0,T2*');removeTicks
% 
% fgmontage(Rec1ccpMt(:,:,3));title('39 shots 2ms, for ref');
%% Longer
TrajPartMed=1:2000;
TrajPartMed=2200+(1:2000);

TrajPartMed=1:4600; % kind of ok
TrajPartMed=1:40000;

TrajPartMed=1:(floor(50000/AcqDwellTime_us/10)*10);

nPointsMed=numel(TrajPartMed);
nTSMed=ceil((nPointsMed+1)/1000);
% nTSMed=8;
TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
% TimePointsMed_ms=0;
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
TSCxMed=R1x.^(TimePointsMed_ms3/InnerShotDiff_ms);
TSCxPMed=permute(TSCxMed,[1 2 7 6 5 4 3]);

clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
    SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
    KernsByRepMed{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRep,TSBMed);
    KernsPMedC{CurRep}=permute(KernsByRepMed{CurRep},[1 2 7 6 5 4 3]);
end
KernsPMMed=cat(3,KernsPMedC{:});
clear KernsPMedC
disp('Prepared Med');
%% given B0,T2*, several reps together
% CurReps=1:2;
% % 2000: 18v10 sec (Normal)
% 
% RegStrMed='-R W:3:0:0.0000001';
% 
% % ScriptFN_TS=[BaseSP 'NuftTSC.txt'];
% % ScriptFN_TS=[BaseSP 'nuftTSC_N.txt'];
% 
% % RecgD_SReps_Med=bart(['picsS -m ' RegStrMed ' ' ScriptFN_TS],Sz16,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
% %         SensCC(:,:,:,1:nccToUse).*TightMask,STraj3MMed(:,:,CurReps),TSBPMed,...
% %         TSCxPMed,sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
%     
% RecgD_SReps_Med=bart(['picsS -m ' RegStrMed ' ' ScriptFN_TS],Sz16,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
%         DMsk.*SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed,...
%         TSCxPMed,sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
% 
% fgmontage(RecgD_SReps_Med);title('Multi-shot, medium-length part of traj, given B0,T2*');removeTicks
% 
% fgmontage(Rec1ccpMt(:,:,1));title('39 shots 2ms, for ref');
% 
% fgmontage(Rec1ccpMtIn(:,:,1));title('39 shots 2ms, in, for ref');
%%
TightMask=DMsk;
%% Try to get all TSC, multi-shot, for B0,T2*
CurReps=1:39;

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC=FillOnesTo16(size(TSCxPMed));

% THLR_lambda=0.1;
THLR_lambda=10;
% RhoStr='';
RhoStr=[' -u ' num2str(1e-3) ' '];

THLRMultiShot=bart(['picsS -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16AllTSC,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed,1,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
%%
[~,~,~,H_AllTS]=ghankel(nTSMed,2,TrgSz);
[ ~, s_LLR_AllTS, V_LLR_AllTS] = batch_svd(H_AllTS*squeeze(THLRMultiShot));
R1ts=V_LLR_AllTS(:,:,2,1)./V_LLR_AllTS(:,:,1,1); % R1 is simply the decay
InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
UpdatedT2SMap_ms1=-InnerTSDiff_ms./log(abs(R1ts));
UpdatedB0Map1=-(angle(R1ts)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
%%
figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
%%
UpdatedB0Map=UpdatedB0Map1;
UpdatedT2SMap_ms=UpdatedT2SMap_ms1;
R1=R1ts;
R1x=min(max(abs(R1),0.3),1).*exp(-1i*angle(R1));
% R1x=R1;
InnerShotDiff_ms=InnerTSDiff_ms;
%%
DataForSlice{SlicesToRead}={UpdatedB0Map,UpdatedT2SMap_ms,UpdatedB0Map0,UpdatedT2SMap_ms0,THLRMultiShot,Rec1ccpMtOutIn,Rec1ccpMt,Rec1ccpMtIn,SelfSens1,sccmtx};

close all;

end % Slices loop
%%
save([mainP filesep 'DataForSlice.mat'],'DataForSlice');
% load([mainP filesep 'DataForSlice.mat']);
%%
clear KernsByRepMed              	1513758480    		1513.76MB
clear KernsPMMed                 	1513754112    		1513.75MB
clear ADataIsL                   	499200000     		499.20MB
clear DataForSlice               	492372000     		492.37MB
clear ADataIsLCC                 	483600000     		483.60MB
clear DataC                      	483600000     		483.60MB
clear DataCCP                    	439640760     		439.64MB
clear DataPC                     	439640760     		439.64MB
clear InnerShotData              	193440000     		193.44MB
clear InnerShotDataIn            	193440000     		193.44MB
save([mainP filesep 'CurStatus.mat']);