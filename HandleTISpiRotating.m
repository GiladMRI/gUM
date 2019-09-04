ScanP='/media/deni/bigdrive/SPIRAL_ASL/S03/RawData/';
BaseFN='meas_MID129_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_FID45302';
RefFldMapP=[ScanP 'meas_MID122_BP_fieldmap_5echosX_FID45295' filesep];

ScanP='/media/deni/bigdrive/SPIRAL_ASL/S04/RawData/';
BaseFN='meas_MID363_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_FID46893';
RefFldMapP=[ScanP 'meas_MID354_BP_fieldmap_9echos_2mm_Full_FID46884' filesep];

ScanP='/autofs/cluster/kawin/Gilad/PhantomBay4/';
BaseFN='meas_MID00656_FID03768_gSpi2d_T20_VD1_800k_d110';
RefFldMapP=[ScanP 'meas_MID00526_FID03646_gre_te7_40' filesep];

MIDStr=BaseFN(6:11);
FN=[ScanP BaseFN '.dat'];
disp('ok');
%%
WhichE=1;


mainP=[ScanP BaseFN '_E' num2str(WhichE)];
mkdir(mainP);

system(['chmod +777 -R ' mainP]);
%% Read raw
AData = mapVBVD(FN);
ADatax=AData{2};
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
%%
load('GAll1mm.mat');
GTrajaC=GAll(:,1,3);
% LL=getLines([ScanP 'SpiGradVec.cpp']);
LL=getLines(['SpiGradVec.cpp']);
% LL=LL(4:end-2);
% RR=LL{11}(2:end-3);
% XX=regexp(RR,'f,','split');
% DD=str2double(XX);
% MaxAmp=DD(end);
% DD(end)=0;
% DD=DD*MaxAmp;
% GTraja=squeeze(PartitionDim(DD,2,2));
% GTrajaC=GTraja(:,1)+1i*GTraja(:,2);
%% QuickAnalyzeTrajJumps;
% figure;subplot(1,2,1);
% plot(kK)
% zIdxs=find((real(kK(1:end-1)).*real(kK(2:end)))<0);
% hold on;
% plot(kK(zIdxs),'ro');
% SI=sort(imag(kK(zIdxs)));
% dSI=diff(SI);
% subplot(1,2,2);
% plot(dSI)
%% Do something with the noise!
noise=AData{1}.noise();
%%
for s=1:nSlices
    try
        SlbLoc(1,s)=asSlice(s).sPosition.dSag;
    catch
        SlbLoc(1,s)=0;
    end
    try
        SlbLoc(2,s)=asSlice(s).sPosition.dCor;
    catch
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
dFOV=FOVx/1000;

FOV_mm=FOVx;
disp('ok');
%%
MaxReps=2;
% ROPerADC Channels 1 1 Slices 3 1 1 Reps 1 ADCs
% ADataIsL=ADatax.image(:,:,:,:,1:nSlicesNoMB,3,:,:,1:MaxReps,:,:,:,:);

ADCsToRead=1:10;
% ADCsToRead=1;
SlicesToRead=[32];
RepsToRead=1:nRepsHdr;
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,ADCsToRead,:,:,:,:,:,:);

ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);
%ADataIsL=AData.image(:,:,:,:,1:nSlicesNoMB,3,:,:,1:5,:,:,:,:);

%ADataIsL=AData.image(:,:,:,:,1:nSlicesNoMB,1,:,:,1,:,:,:,:);
% Oredr: [1024 Channels 1 1 Slices 1 1 1 Reps 1 ADCs];
% ADataIsL=squeeze(ADataIsL);
disp('Read data');
%%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_ms=10e-3;

g=GTrajaC;
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

kK=k*FOV_mm/1000/2/pi;
figure;
subplot(2,3,1);
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
CurSens=SensBx(:,:,:,18);
%%
SensCC=permute(sum(CurSens.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
disp('ok SensCC');
%%
TrgSz=[220 220];
OnesSens=ones(TrgSz);
%%
Traj=interp1(1:5001,kK(1:5001),1:0.125:5001-0.01);

Extra=4;
kKA=[zeros(4,1); kK];
Reduce=21;
Delay=-6; % clockwise
Traj=[interp1(0:2500+Extra,kKA(1:2501+Extra),Extra-Delay/10+linspace(0,2500-Reduce,20000)) interp1(0:2500+Extra,kKA(2500+(1:2501+Extra)),Extra-Delay/10+linspace(0,2500-Reduce,20000))]; 
%%
figure;plot(grmss(ADataIsL,2));hold on
plot(abs(Traj)/50000,'k');
xlim([0 70])
xlim([4900 5100])
xlim(5000+[4800 5200])
xlim(10000+[4800 5200])
xlim(15000+[4800 5200])
xlim(15000+[4900 5100])
%%
figure;
Ranges=[0 70; 4900 5100; 5000+[4800 5200]; 5000*2+[4800 5200]; 5000*3+[4800 5200]; 5000*4+[4800 5200]; 5000*5+[4800 5200]; 5000*6+[4800 5200]];
for i=1:8
    subplot(2,4,i);
    plot(grmss(ADataIsL,2));hold on
    plot(abs(Traj)/50000,'k');
    xlim(Ranges(i,:))
end
title(Reduce);
% xlim(15000+[4970 5030])
%%
for i=1:size(ADataIsL,3)
    tmp=grmss(ADataIsL(5000:5060,:,i),2).';
    [~,PeakIdx]=max(tmp);
    Peak(i)=PeakIdx;
end
figure;
for i=1:6
    subplot(2,3,i);
    tmp=grmss(ADataIsL(5000:5060,:,i+13),2).';
    plot(tmp);hold on;
%     CenterMass=sum(tmp.*(1:numel(tmp)))/sum(tmp);
%     plot([CenterMass CenterMass],[0 max(tmp)]);
    [~,PeakIdx]=max(tmp);
    plot([PeakIdx PeakIdx],[0 max(tmp)]);title(PeakIdx);
end
%%
nRepsToUse=30;
TrajM=Traj.*exp(-1i*2*pi*110/360*(0:nRepsToUse-1)');
% Grad delay
% GradDelayPoints=0; % positive clockwise
TrajMa=TrajM;
DataC=permute(ADataIsLCC,[1 3 2]);
% DataC=circshift(DataC,-GradDelayPoints,1);
% TrajMa=[zeros([size(TrajM,1) GradDelayPoints]) TrajM(:,1:end-GradDelayPoints)];
%
% DataC=permute(ADataIsL,[1 3 2]);
%
TrajM=TrajMa*exp(1i*2*pi*-10/360)/1.03; % positive counter clockwize

BARTTraj=cat(3,real(TrajM),imag(TrajM),imag(TrajM)*0);
BARTTrajP=permute(BARTTraj,[3 2 1]);

dx=RotatedLocs(2,1)/ADatax.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
dy=RotatedLocs(1,1)/ADatax.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;

kx=BARTTrajP(1,:,:)*2*pi;
ky=BARTTrajP(2,:,:)*2*pi;

% modx=conj(exp(1i.*(dx*kx+dy*ky)));
modx=exp(-1i.*(dx*kx+dy*ky));
%
DataPC=permute(DataC,[4 1 2 5 6 7 8 3]).*modx;
% DataP=permute(Data,[3 1 2]).*modx;

% figure;plot(TrajM(:,:).')
% figure;plot(TrajM(:,1:8000).')
%
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 ncc]);
% Rec1c=bart('pics -t ',BARTTrajP,(DataPC),OnesSensC);
% fgmontage(Rec1c,[0 1e-3])
% ShowAbsAngle(Rec1)
%
TrajPartToUse=30000+(1:2000);
RepsToUse=1:30;
Rec1cp=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
% fgmontage(Rec1cp)
fgmontage(grmss(Rec1cp,8));title(GradDelayPoints);
% fgmontage(Rec1cp,[0 1e-3])
%%
nccToUse=31;
RepsToUse=1:30;
% TrajPartToUse=[1:5000 20000:25000];
TrajPartToUse=30000+(1:2500);

DataCCP=permute(DataPC,[1:3 8 4:7]);
% Rec1ccp=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,:),SensCC);
% fgmontage(Rec1ccp)
% ShowAbsAngle(Rec1ccp)
Rec1ccp=bart('pics -m -R W:3:0:0.00010 -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,1:nccToUse),SensCC(:,:,:,1:nccToUse));
fgmontage(Rec1ccp)
%%
for i=1:4
    nccToUse=31;
    RepsToUse=1:30;
    TrajPartToUse=(i-1)*10000+(1:2500);
    TrajPartToUse=(i-1)*10000+(1:1000);
    
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    Rec1ccpM(:,:,i)=bart('pics -m -R W:3:0:0.00010 -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,1:nccToUse),SensCC(:,:,:,1:nccToUse));
%     Rec1ccpM(:,:,i)=bart('pics -m -R W:3:0:0.10 -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,1:nccToUse),SensCC(:,:,:,1:nccToUse));
%     Rec1ccpM(:,:,i)=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,1:nccToUse),SensCC(:,:,:,1:nccToUse));
end
%%
HankelTemporalLen=2;
WhichIdxs=1:4;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(numel(WhichIdxs),HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[TrgSz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[TrgSz, numel(HankelMat)]).*DeHankelizingMatP,3));
%%
[ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(Rec1ccpM(:,:,WhichIdxs)));

R1=V_LLR(:,:,2,1)./V_LLR(:,:,1,1); % R1 is simply the decay
% fgmontage(angle(R1));
% fgmontage(-log(abs(R1)),[0 0.2])
% fgmontage(-1./log(abs(R1)),[0 100]) % TE is ~1ms, so ignore

TotalTime_ms=30000*1.25/1000;
UpdatedB0Map=-angle(R1)*1e3*(4/TotalTime_ms)/(2*pi);

% R1s{Iter}=R1;
% 
% Updated_TSC_L2=exp(-1i*angle(R1).*permute(0:(nTS_L2-1),[1 3 2]));
% % Updated_TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);

%%
%%
ADataIsLCCx=ADataIsLCC;

for i=1:size(ADataIsL,3)
    ADataIsLCCx(:,:,i)=circshift(ADataIsLCCx(:,:,i),33-Peak(i),1);
end
 
DataC=permute(ADataIsLCCx,[1 3 2]);

DataPC=permute(DataC,[4 1 2 5 6 7 8 3]).*modx;

figure;
subplot(1,2,1);plot(grmss(ADataIsLCC(5000:5060,:,:),2))
subplot(1,2,2);plot(grmss(ADataIsLCCx(5000:5060,:,:),2))
%%
RepsToUse=1:30;
% TrajPartToUse=[1:5000 20000:25000];
nPointsToUse=400;
TrajPartToUse=0000+(1:nPointsToUse);

% DataPC=permute(circshift(DataC,2,1),[4 1 2 5 6 7 8 3]).*modx;

DataCCP=permute(DataPC,[1:3 8 4:7]);
Rec1ccp1=bart('pics -m -R W:3:0:0.00010 -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,:),SensCC);
TrajPartToUse=10000+(1:nPointsToUse);
DataCCP=permute(DataPC,[1:3 8 4:7]);
Rec1ccp2=bart('pics -m -R W:3:0:0.00010 -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),DataCCP(:,TrajPartToUse,RepsToUse,:),SensCC);
fgmontage(Rec1ccp1);MaximizeFig
fgmontage(Rec1ccp2);MaximizeFig
%%
Rec1=bart('pics -t ',BARTTrajP,conj(DataP),OnesSens);
ShowAbsAngle(Rec1)
%%
Rec1=bart('pics -t ',BARTTrajP(:,1:18000,:),DataP(:,1:18000,:),OnesSens);
ShowAbsAngle(Rec1)
%%
Rec2=bart('pics -u -R T:3:3:1 -t ',BARTTrajP,DataP,OnesSens);
ShowAbsAngle(Rec2)
%%

%%
% file 0 is Traj - [3 #Traj]
% file 1 is Sens - [x y 1 Coils 1 1 MB]
% file 2 is TS image - [X Y 1 1 1 1 MB #TS]
% file 3 is TS signal - [1 #Traj 1 Coils 1 1 MB #TS]
% Data is [1 #Traj 1 Coils]
%%
paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{20};
spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{6};
paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{5};
paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};

CAIPISep_mm=AData.hdr.MeasYaps.sWiPMemBlock.adFree{7};
CAIPIPeriod_us=AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
CAIPIDelay_us=AData.hdr.MeasYaps.sWiPMemBlock.adFree{9};

% if(MB>1)
if(isempty(spBW))
    spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{14};
    paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
    paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};
    paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
end
disp('Read data');
%%
OutInOut=paramLongROSamples>9000;

% ADataIsP=double(permute(ADataIsL,[1 5 2 3 4]));
ADataIsPy=permute(ADataIsL(:,:,1:nSlicesNoMB,:,:),[1 5 2 3 4]);

nReps=size( ADataIsPy,5);
nChannels=size( ADataIsPy,3);

ADataIsPy=reshape(ADataIsPy,[paramLongROSamples,nChannels nSlicesNoMB nReps]);
if(OutInOut)
    ADataIsPy=PartitionDim(ADataIsPy,1,2);
    ADataIsPy=ADataIsPy(:,:,:,:,WhichE);
end
% ADataIsPy=ADataIsPy(1:size(BARTTrajMS,2),:,:,:,:);
disp('ok b')
%%
clear ADataIsL     
%%
SummaryStr.MB=MB;
SummaryStr.paramLongSpGradAmp=paramLongSpGradAmp;
SummaryStr.paramLongSpSlewRate=paramLongSpSlewRate;
SummaryStr.paramLongROSamples =paramLongROSamples;
SummaryStr.spBW=spBW;
SummaryStr.AccR=AccR;
SummaryStr.paramLongInterleaves=paramLongInterleaves;
SummaryStr.VD=VD;

if(MB>1)
    SummaryStr.CAIPISep_mm=CAIPISep_mm;
    SummaryStr.CAIPIPeriod_us=CAIPIPeriod_us;
    SummaryStr.CAIPIDelay_us=CAIPIDelay_us;
end
%%
% save([mainP filesep 'Data.mat'],'ADataIsL','AData');
% disp(['Saved ' [mainP filesep 'Data.mat']]);
%%
% save('StatusForTrajDesign1Band_2.mat');
%% given everything but VD and Acc
% OutInOut=false;
if(OutInOut)
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);
else
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
        paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
end

CAIPIVec=CAIPIBlips([paramLongSpGradAmp, paramLongSpSlewRate,CAIPISep_mm,CAIPIDelay_us,CAIPIPeriod_us,...
    2560*paramLongROSamples/1024]);

gamma=42.5774806;
cCAIPIVec=cumsum(CAIPIVec)*gamma*10*2*pi/1e6;

if(OutInOut)
    BaseLen=size(kTraj,1)/3;
    TrajIdxs=(1:BaseLen)+BaseLen*(WhichE-1);
    kTraj=kTraj(TrajIdxs,:);
    GradBuf=GradBuf(TrajIdxs,:);
    
    CAIPIVec=CAIPIVec(TrajIdxs,:);
    cCAIPIVec=cCAIPIVec(TrajIdxs,:);
end

EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));


cCAIPIVec=cumsum(CAIPIVec)*gamma*10*2*pi/1e6;
cCAIPIVecX=interp1(1:numel(cCAIPIVec),cCAIPIVec,1:1e5/spBW:numel(cCAIPIVec));
% 
% %% Perparing for comparison to measured trajectory
% OrigTraj=[kTraj cCAIPIVec*1000];
% for i=1:3
%     OrigTrajQ(:,i)=interp1(1:size(OrigTraj,1),OrigTraj(:,i),1:1e5/spBW:(size(OrigTraj,1)-0.0001));
% end
% 
% %% Here we look for starting point
% 
% for StartPoint=180
%     
% RFPoints=10*(0:1e5/spBW:size(OrigTraj,1)-1.001);
% clear IData
% for i=1:3
%     IData(:,i)=interp1(1:size(RData3,1),RData3(:,i),StartPoint+(RFPoints));
% end
% ChannelsToFlip=[1 3];
% 
% 
% IData(:,ChannelsToFlip)=-IData(:,ChannelsToFlip);
% % Showing with exaggerated CAIPI blips
% MidKz=max(abs(OrigTrajQ(:,3)))/2;
% OrigTrajQX=OrigTrajQ;
% OrigTrajQX(:,3)=(OrigTrajQX(:,3)-MidKz)*20;
% 
% IDataX=IData;
% IDataX(:,3)=(IDataX(:,3)-MidKz)*20;
% XsToShow=1:100;
% XsToShow=1:size(OrigTrajQ,1);
% figure;
% plot(OrigTrajQX(XsToShow,:));hold on;
% plot(IDataX(XsToShow,:),'--');
% setXaxis([1 XsToShow(end)]);
% title('Planned (solid) and measured (dashed) tajectory (after rotation optimization)');
% end
% %% Now we replace the designed with the measured, for the rest of the recon
% %kTrajQ=IData(:,1:2);
% kTrajQ=IData(:,[2 1]);
% cCAIPIVecX=IData(:,3)/1000;
% 
% figure; plot(kTrajQ(:,1)); hold on;
% plot(kTrajQ(:,2)); hold on; plot(cCAIPIVecX);

deltaT=1e6/spBW;
ADCAdvantagePoints=round(5/deltaT);
StartPoint=ADCAdvantagePoints+1;
cCAIPIVecY=cCAIPIVecX(StartPoint:end);
%%
BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;

BARTTrajMS=BARTTrajx;
BARTTrajxC=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
for i=2:paramLongInterleaves
    CurRotTrajC=BARTTrajxC*exp(1i*2*pi*(i-1)/paramLongInterleaves);
    BARTTrajMS(1,:,i)=real(CurRotTrajC);
    BARTTrajMS(2,:,i)=imag(CurRotTrajC);
end
BARTTrajMS(3,end)=0;

nTrajOrig=size(ADataIsPy,1);
nTraj=nTrajOrig-ADCAdvantagePoints;
nTraj=5116*2;
TimeInMs2=(0:nTraj-1)*deltaT/1e3;
T2SCompStr='';

BARTTrajAct=BARTTrajMS(:,1:(nTraj));
BARTTrajCorr=BARTTrajAct./[1.01; 1.01; 1.0];

disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
% nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
XResmm=dFOV*1000/(2*MaxK);
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));

setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc) ]);
xlabel(['Res ' num2str(XResmm,'%.2f') 'mm']);
ylabel(['Nominal acc: ' num2str(AccR) ' VD=' num2str(VD)]);

subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['                Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
hold on;
plot(CAIPIVec,'k');
subplot(2,2,4);
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
subplot(2,2,3);
plot(cumsum(CAIPIVec));
title('CAIPI cumsum');
%% Fancy plotting 3D
% Traj3=BARTTrajMS;
% 
% CAIPIkI = interp1(1:size(kTraj,1),cumsum(CAIPIVec),1:1e5/spBW:(size(kTraj,1)-0.01));
% 
% Traj3(3,:)=CAIPIkI;
% figure;plot3(Traj3(1,:),Traj3(2,:),Traj3(3,:))
% B=CAIPIkI<(max(CAIPIkI)+min(CAIPIkI)/2);
% Traj3X=Traj3;
% Traj3X(:,B)=NaN;
% figure;plot(Traj3X(1,:),Traj3X(2,:),'LineWidth',2);hold on
% Traj3X=Traj3;
% Traj3X(:,~B)=NaN;
% plot(Traj3X(1,:),Traj3X(2,:),'LineWidth',2);
% setXaxis([-1.1 1.1]*ceil(MaxK));
% setYaxis([-1.1 1.1]*ceil(MaxK));
%%
gprint(get(gcf,'Number'),[mainP filesep 'Traj'],[]) 
close(gcf);
save([mainP filesep 'Traj.mat'],'BARTTrajMS');
disp('Saved traj fig');
%% Summary P2
SummaryStr.MaxK=MaxK;
SummaryStr.nTraj=nTraj;
SummaryStr.Acc=Acc;
SummaryStr.XResmm=XResmm;
SummaryStr.nSlices=nSlices;


fileName = ['Summary_' BaseFN(6:11) '.txt'];
gStruct2txt(SummaryStr,[mainP filesep fileName]);
disp('Saved summary txt');
%%
Trajm2=BARTTrajMS(1:2,1:end-2);

Sz128=[96 96];

[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
% save('ForTFNUFT.mat','SN','Kd','P','A','NUbyFS3');
save([mainP filesep 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');
disp('Saved TrajForNUFT');
%%
% Traj=BARTTrajMS(1:2,1:end-2);
% save([ScanP BaseFN filesep 'Traj.mat'],'Traj');
%%
RLocs=load([RefFldMapP 'Locs.mat']);
%%
SensB=load([RefFldMapP 'Sens.mat']);
SensB=SensB.SensB;

SensBx=padarray(SensB(:,:,:,:,1),[0 2 0 0],'both');
SensBx=rot90(imresizeBySlices(SensBx,TrgSz));


FirstEcho=load([RefFldMapP 'FirstEcho.mat']);
FirstEcho=FirstEcho.FirstEcho;

FirstEcho=gflip(FirstEcho(:,:,:,SliOffset+(1:nSlices)),[]);
Mg=grmss(FirstEcho,3);

Mgx=padarray(Mg,[0 2 0],'both');
Mgx=rot90(imresizeBySlices(Mgx,TrgSz));
disp('ok');
%%
B0S=load([RefFldMapP 'B0S.mat']);
B0S=B0S.B0S;

T2S=load([RefFldMapP 'T2S.mat']);
T2S=T2S.T2SAll;
disp('ok');
% load('CurStatus_Ben4MinASL_x.mat','SensB','B0Q','Mg'); ;
%%
% SensX=SensB(:,:,:,6+(1:12),:);
SensX=permute(SensB(:,:,:,SliOffset+(1:nSlices),:),[1 2 3 5 4]);
SensX=gflip(SensX,1:2);
Sz2=gsize(SensX,1:2);
% B0Q2=B0Q(:,:,6+(1:12));
B0Q2=B0S(:,:,SliOffset+(1:nSlices));
B0Q2=gflip(B0Q2,1:2);

T2S=imresizeBySlices(T2S,gsize(SensX,1:2));
T2S=gflip(T2S,1:2);

MgcS=imresizeBySlices(Mg,Sz2);
MskcS=MgcS>7e-5;

MskcSE=imdilate(imfillholesBySlices(MskcS),strel('disk',5,8));
B0M2S=B0Q2;
B0M2S(~MskcSE)=0;

disp('ok');            
%%
try
    dx=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
catch
    disp('No x shift!');
    dx=0;
end

dx=RotatedLocs(2,1)/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
dy=RotatedLocs(1,1)/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;

kx=BARTTrajMS(1,:)*2*pi;
ky=BARTTrajMS(2,:)*2*pi;

modx=double(exp(1i*(dx*kx+dy*ky))');


%Sz2=gsize(SensX,1:2);
disp('ok k');
%% SensCCS for all slices
nBands=MB;
WhichSlices=nSlicesNoMB:-1:1;
for SliI=WhichSlices
    SliIs=[SliI SliI+nSlicesNoMB];

    nukData=ADataIsPy(:,:,SliI,min(25,size(ADataIsPy,4)-1)).';
    nukData=nukData(:,StartPoint:end);
    nukData=double(nukData);

    [U,S,sccmtx] = svd(nukData.','econ');
    
    sccmtxS(:,:,SliI)=sccmtx;
    
    ncc=16;

    for b=1:nBands
        CurSens=SensX(:,:,:,1,SliIs(b));
        SensCC=permute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2 4])),[2 3 1 4]);
        SensCCS(:,:,:,:,SliIs(b))=SensCC;
    end
end
disp('SensCCS ok');
%%
save([mainP filesep 'sccmtxS.mat'],'sccmtxS');
disp('sccmtxS saved');
%% All B0 effects across time
MgF=imresizeBySlices(Mg,Sz2);
MskcF=MgF>7e-5;
MskcEF=imdilate(imfillholesBySlices(MskcF),strel('disk',5,8));
B0MF=B0Q2;
B0MF(~MskcEF)=0;
% fgmontage(B0MF,[-500 500])
disp('ok Loaded B0 from field map, all slices');
%% TS for all slices
clear TSBC TSCC
MaxCond=1e4;
for WhichB=1:nSlices
    [TSBC{WhichB},TSCC{WhichB},nTSXC(WhichB),CondGHG]=GetTSCoeffsByMinMaxHist(B0MF(:,:,WhichB),TimeInMs2,500,MaxCond);
    
    TSBC{WhichB}=GetTSCoeffsByLinear(nTraj,nTSXC(WhichB));
    
    FesTimePoints=linspace(TimeInMs2(1),TimeInMs2(end),nTSXC(WhichB));
    FesTimePointsP=permute(FesTimePoints,[1 3 2]);
    CurT2S=T2S(:,:,WhichB);
    CurT2S(CurT2S<4)=20;
    TSCC{WhichB}=TSCC{WhichB}.*exp(-FesTimePointsP./CurT2S);
end
nTSX=max(nTSXC);
TSC=ones([Sz2 nTSX 2]);
TSB=zeros([nTraj nTSX 2]);
clear TSCS TSBS
for WhichB=1:nSlices
    TSCS(:,:,1:nTSXC(WhichB),WhichB)=TSCC{WhichB};
    TSBS(:,1:nTSXC(WhichB),WhichB)=TSBC{WhichB};
end
nTS=nTSX;
disp('ran GetTSCoeffsByMinMaxHist on all')
%% GPU TS
% cCAIPIVecX=interp1(1:numel(cCAIPIVec),cCAIPIVec,1:1/4:numel(cCAIPIVec));
% cCAIPIVecY=cCAIPIVecX(3:end-1);
% cCAIPIVecZ=exp(1i*cCAIPIVecY.*(RotatedLocs(3,SliIs)'));

DLocs=RotatedLocs(3,[1 1+nSlicesNoMB])';
DLocs=DLocs-DLocs(1);
cCAIPIVecZ=exp(1i*cCAIPIVecY'.*DLocs);

osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

NTrg1=Sz2;
nBands=MB;


nShots=paramLongInterleaves;
disp('ok')
% TSBF=permute(TSB(IB_TimeInMs2,:,:),[2 1 3]);
% Sens=squeeze(SensX(:,:,:,1,SliIs));
%%
Delay=(WhichE-1)*5;

nTrajS=nTraj-Delay;
nTrajS=5114;
for SliI=1:nSlicesNoMB
    SliIs=[SliI SliI+nSlicesNoMB];
    
    TSBF=permute(TSBS(1:nTrajS,:,SliIs),[2 1 3]);
    Sens=squeeze(SensCCS(:,:,:,:,SliIs));
    TSC=TSCS(:,:,:,SliIs);
    
    GOP_MCSMBMSC{SliI} = ggpuNUFT_TS_MC_MB_MS_CAIPI(BARTTrajCorr(:,1:nTrajS),NTrg1,osf,wg,sw,TSBF,TSC,Sens,cCAIPIVecZ(:,1:nTrajS));
end
% GOP_MC=GOP_MCSMBMS;
nTrajP2=nTrajS;
disp('ok');
% %%
% SliIP=[mainP filesep YLbl filesep];
% mkdir(SliIP);
% save([SliIP 'TSBF_TSC_Sens.mat'],'TSBF','TSC','Sens');
% disp(['Saved ' SliIP 'TSBF_TSC_Sens.mat']);
% %%
GOP_MC=GOP_MCSMBMSC{6};
x = randn([Sz2 MB]) + 1j*randn([Sz2 MB]);
y = randn([size(Sens,3) nTrajS]) + 1j*randn([size(Sens,3) nTrajS]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
if(abs(Out)>1e-3)
    disp('Operator error?');
end
%%
tic
% if(nShots==1)
    TVOP=TVOP_MSlice;
% else
%     TVOP=TVOP_MTC_W([1 1 0 1e1]);
% end
WhichReps=min(25,nReps);
WhichRep=WhichReps;
WhichReps=1:nReps;

% WhichSliIs=1:nSlicesNoMB;
WhichSliIs=1:12;
for WhichRep=WhichReps
    for SliI=WhichSliIs

        disp([SliI WhichRep]);
        
        SliIs=[SliI SliI+nSlicesNoMB];
        
        AOdd = GOP_MCSMBMSC{SliI};
        
        nukData=ADataIsPy(:,:,SliI,WhichRep).';
        nukData=double(nukData);
        kx=BARTTrajCorr(1,:)*2*pi;
        ky=BARTTrajCorr(2,:)*2*pi;
        
        modx=double(exp(1i*(dx*kx+dy*ky))');
        modx(size(nukData,2),1)=0;
        
        nukData=circshift(nukData,-(Delay+2),2);
        nukData=nukData.*(modx.');
        
        % CAIPI Slice base
        nukData(:,1:nTraj)=nukData(:,1:nTraj).*exp(-1i*cCAIPIVecY.*(RotatedLocs(3,SliI)'));
        
        nukData=nukData(:,1:end-Delay-2);
        
        %         nukData=nukData(:,(Delay+3):end);
        
        sccmtxBoth = sccmtxS(:,:,SliI);
        nScc=16;
        %     SensFCCBoth=squeeze(SensCCS(:,:,:,:,SliIs));
        
        DataCC=CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc));
        nukDataP=permute(DataCC,[2 1]);
        DataP=nukDataP;
        
        TVW=1e-6;
        
        % filterType:   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies','Symmlet', 'Vaidyanathan','Battle'
        % Use suffix _TI for translation invariance, for example, 'Daubechies_TI'
        % filterSize: related to the support and vanishing moments of the particular wavelet (See MakeONFilter in wavelab)
        % wavScale: 	scallest scale of wavelet decomposition
        
        XFMStr='Daubechies';
        %     XFMStr='Daubechies_TI';
        XFMStr='None';
        
        filterSize=4;
        wavScale=4;
        
        if(strcmp(XFMStr,'None'))
            XFM=1;
            xfmWeight=0;
        else
            XFM = Wavelet(XFMStr,filterSize,wavScale);
            xfmWeight = 1e-6;	% Weight for Transform L1 penalty
        end
        
        
        param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);
        
        param.data =     DataP;
        
        param.ShowFig=false;
        
        nfnlCgIters=40;
        RunFnlViewAmp=1;
        res=zeros([Sz2 MB]);
        if(nShots>1)
            res=repmat(resA,[1 1 1 nShots]);
            res=res+randn(size(res))*max(abs(resA(:)))/20;
        end
        
        FigH=4000;
        if(param.ShowFig)
            figure(FigH);close(FigH);
        end
        
        if(~isfield(param,'ShowFig'))
            param.ShowFig=true;
        end
        StartTime_fnl=now;
        param.Verbose=false;
        clear ObjConv Score
        for n=1:nfnlCgIters
            [res, CurObj] = fnlCg(res,param);
            ObjConv(n)=CurObj;
            im_res = param.XFM'*res;
            if(param.ShowFig)
                figure(FigH); subplot(1,3,1);
                gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
                cx=caxis;
                caxis(cx/RunFnlViewAmp);
                subplot(1,3,2);
                gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
                subplot(1,3,3);
                plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
                
            end
            if(n>1)
                dObjP=(ObjConv(n-1)-CurObj)/ObjConv(n-1);
                disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
                if(dObjP<2e-3)
                    disp('Not advancing. Stopping.');
                    break;
                end
            else
                disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
            end
            
            if(nShots==1)
                resA=res;
            end
        end
        if(param.ShowFig)
            close(FigH);
        end
        disp('ok im_res');
        
        ImResR(:,:,SliIs,WhichRep)=im_res;
    end
end
DefRange=[0 7e-3];
if(WhichE==2)
    DefRange=[0 5e-3];
end
save([mainP filesep 'Sparse_recon.mat'],'ImResR');
toc
fgmontage(im_res,DefRange);

XFMStrFull=['[' XFMStr ',' num2str(filterSize) ',' num2str(wavScale) ',' num2str(xfmWeight) ']'];
XLbl=['L' num2str(param.pNorm) ',TVW=' num2str(param.TVWeight) ',' XFMStrFull ];
xlabel(XLbl)
YLbl=['Sli' num2str(SliI,'%02d')];
ylabel(YLbl);
%%
gprint(get(gcf,'Number'),[mainP filesep YLbl '_' XLbl T2SCompStr],[]) 
close(gcf);
save([mainP filesep YLbl '_' XLbl T2SCompStr '.mat'],'im_res');
disp(['Saved' mainP filesep YLbl '_' XLbl T2SCompStr '.mat']);
%%
save([mainP filesep 'ImResR.mat'],'ImResR');
%%
% pause(60*60*3+200);
%cCAIPIVecZ=cCAIPIVecZ(:,1:end-5);
nTraj=5114;
for SliI=1:12
   

    SliIs=[SliI SliI+nSlicesNoMB];

nukData=ADataIsPy(:,:,SliI,WhichRep).';
kx=BARTTrajCorr(1,:)*2*pi;
ky=BARTTrajCorr(2,:)*2*pi;

modx=double(exp(1i*(dx*kx+dy*ky))');
modx(size(nukData,2),1)=0;

nukData=circshift(nukData,-(Delay+2),2);
nukData=nukData.*(modx.');

% CAIPI Slice base
nukData2(:,1:nTraj)=nukData(:,1:nTraj).*exp(-1i*cCAIPIVecY'.*(RotatedLocs(3,SliI)'));

%nukData=nukData2(:,1:end-Delay-2);
nukData=nukData2;

%         nukData=nukData(:,(Delay+3):end);

sccmtxBoth = sccmtxS(:,:,SliI);
nScc=16;
%     SensFCCBoth=squeeze(SensCCS(:,:,:,:,SliIs));

DataCC=CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc));
nukDataP=permute(DataCC,[2 1]);
DataPCC=double(permute(DataCC,[3 2 1]));
disp('Applied coil compression to this slice(s)'); 
%Trajm2=BARTTrajCorr; % BARTTrajMS(1:2,1:end-2);   %before correction
% Sz128=[128 128];
%Trajm2=Trajm2(:,1:nTrajS); %after correction



[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
% save('ForTFNUFT.mat','SN','Kd','P','A','NUbyFS3');
save([mainP filesep 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');
disp('Saved TrajForNUFT');
% Clean save 
TSBF=permute(TSBS(1:nTrajS,:,SliIs),[2 1 3]);
%TSBF=permute(TSBS(:,:,SliIs),[2 1 3]);
TSC=TSCS(:,:,:,SliIs);
SensCC=squeeze(SensCCS(:,:,:,:,SliIs));

SensCCA=SensCC(:,:,:,1);
SensCCB=SensCC(:,:,:,2);
SensMskA=single(grmss(SensCCA(:,:,:),3)>0.01);
SensMskB=single(grmss(SensCCB(:,:,:),3)>0.01);



TSBFA=TSBF(:,:,1).*cCAIPIVecZ(1,:);
TSBFB=TSBF(:,:,2).*cCAIPIVecZ(2,:);
TSCA=TSC(:,:,:,1);
TSCB=TSC(:,:,:,2);

CurOurP=[mainP filesep];
mkdir(CurOurP);
CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];
mkdir(CurOurPSli);
% save([CurOurPSli 'Sens.mat'],'CurSens');
% save([CurOurPSli 'SensCC.mat'],'SensCC','sccmtx','SensMsk');
save([CurOurPSli 'SensCC1.mat'],'SensCCA','SensCCB','sccmtx','SensMskA','SensMskB');
save([CurOurPSli 'B0TS.mat'],'TSBFA','TSBFB','TSB','TSCA','TSCB','osf','wg','sw','TimeInMs2');
% save([CurOurPSli 'TrajAndRealData.mat'],'BARTTrajAct','DataP2','DataPCC');
save([CurOurPSli 'TrajAndRealData.mat'],'BARTTrajAct','DataPCC');
disp('Saved');
%
% tmp=squeeze(DataPCC);
tmp=squeeze(DataCC);
% RealDataFac=grmss(AllData(5656:34:6767,:,:))/grmss(tmp)
RealDataFac=100;
CurIDataV=Row(tmp)*RealDataFac;
CurIDataVR=[real(CurIDataV) imag(CurIDataV)];

BatchSizeInNN=12;
Data=repmat(single(CurIDataVR),[BatchSizeInNN 1]);
if(WhichE==2)
    Data=Data*2.5;
end

RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
save(RealDataFN,'Data');
disp('Saved real data');


end
%% At this point ready for TF
% Call TF
% pause(60*60*11)
for SliI=[7];
CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];

%ParamsSDefaults=getParamsStructFromFN('/media/a/H2/home/a/TF/srez/RegridTry3C2_7TS_XXX_Sli06__2018-09-20_13-34-00_train/');
% ParamsSDefaults=getParamsStructFromFN('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/');
ParamsSDefaults=getParamsStructFromFN('/home/deni/TF/');
ParamsS=ParamsSDefaults;

ParamsS.LabelsH=Sz2(1);
ParamsS.LabelsW=Sz2(1);
ParamsS.aDataH=Sz2(1);
ParamsS.aDataW=Sz2(1);
ParamsS.kMax=48;

ParamsS.SessionNameBase=['RegridTry3C2_7TS_' MIDStr '_Sli' num2str(SliI,'%02d') ];

ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
ParamsS.BaseTSDataP=CurOurPSli;
ParamsS.BaseNUFTDataP=[mainP filesep];


%if(strhas(UserName,'deni'))
    ParamsS.DatasetMatFN='/home/deni/HCPData_256x256_int16.mat';
%else
 %   ParamsS.DatasetMatFN='/media/a/H1/HCPData_256x256_int16.mat';
%end
ParamsS.nToLoad=10000;

% Parames for MB:p
ParamsS.nccInData=ncc;

%ParamsS.InputMode='RegridTry3FMB';
%ParamsS.NetMode='RegridTry3C2FT_TS_MB';
%ParamsS.aNetMode='RegridTry3C2_TS_MB';
%ParamsS.NetMode='RegridTry3C2_TS_MB';


%ParamsS.InputMode='RegridTry3F';
%ParamsS.NetMode='RegridTry3C2FT_TS';
%ParamsS.kSideNet='64,5,62,7,14,5';

% To change and decide
ParamsS.nccToUse=16;    %13
ParamsS.nNeighbors=12;  %12
ParamsS.nTimeSegments=5; %11

ParamsS.batch_size=BatchSizeInNN;
ParamsS.WL2_Lambda=0;

ParamsS.RandomPhaseLinearFac=3;
ParamsS.RandomPhaseQuadraticFac=0.005;
ParamsS.RandomPhaseScaleFac=2;

ParamsS.checkpoint_period=20;

ParamsS.InitForRFN='None';
ParamsS.InitForLFN='None';

%ParamsS=getParamsStructFromFN([LastP filesep]);
ParamsS.learning_rate_start=0.003;
ParamsS.learning_rate_half_life=30;
ParamsS.train_time=120;

%noise
ParamsS.noise_level=0.007;

 %ParamsS.InitForRFN='/home/deni/TF/srez/RegridTry3C2_7TS_MID543_Sli06__2018-09-25_13-32-53_train/TrainSummary_023584.mat';
%ParamsS.batch_size=12;

ParamsS.BankSize=0; % Not use a signal bank
% ParamsS.BankSize=1440;
% ParamsS.BankK=7;

ParamsS.ShowRealData=1;
ParamsS.DataH=size(Data,2);
%ParamsS.ISide2Net=0;
% ParamsS.DataH=size(Data,2);

%Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    Txt=gStruct2txt(ParamsS,'~/TF/Params.txt');

CImages_b=ParamsS.nToLoad*256*256*2;
Reg_batch_b=128*128*ParamsS.nNeighbors*ParamsS.nccToUse*4*2*ParamsS.batch_size;
Kside_b=ParamsS.batch_size*128*128*ParamsS.nNeighbors*ParamsS.nccToUse*ParamsS.nTimeSegments*2*4*2;
SBank_b=ParamsS.DataH*4*ParamsS.BankSize*2;
IBank_b=128*128*4*2*ParamsS.BankSize*2;
Ttl_b=CImages_b+Reg_batch_b+Kside_b+SBank_b+IBank_b;
disp(['Constant images: ' num2str(CImages_b/1e9,'%.2f') 'GB']);
disp(['Regridded batch data ' num2str(Reg_batch_b/1e6,'%.2f') 'MB']);
disp(['Kside ' num2str(Kside_b/1e9,'%.2f') 'GB']);
disp(['Signal Bank ' num2str(SBank_b/1e9,'%.2f') 'GB']);
disp(['Image Bank ' num2str(IBank_b/1e6,'%.2f') 'MB']);
disp(['Total ' num2str(Ttl_b/1e9,'%.2f') 'GB']);
disp(['Bank batches: ' num2str(ParamsS.BankSize/ParamsS.batch_size)]);
%
%system('sudo -H -u a /media/a/H2/home/a/RunTFForMatlabx.sh');
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
system('/home/deni/RunTFForMatlabx.sh');
end
% system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
% system('sudo -H -u a /media/a/H2/home/a/RunTFForMatlabx.sh');

% /home/deni/TF/srez/RegridTry3C2_7TS_g18_Sli06__2018-08-16_12-35-57_train
% gSendMail('Subject',Txt,Image FileName,EmailAddr)
%% Look at the convergence graph
LastTrainedNet=[TFFolder FindNewestFolderWithPrefix(TFFolder,'Re','_train')];
GraphOptFromFolderf(LastTrainedNet,true);
%% Look at the image-side segments combination coefficient maps
W=getWeightsOfNetC(LastTrainedNet);
WIm=reshape(W.gene_GEN_L010_PixelwiseMultC_weightC,128,128,2,[]);
%%
TFFolder='/home/deni/TF/srez/';
TrainedNetP=[TFFolder FindNewestFolderWithPrefix(TFFolder,'Re','_check')];
TrainedNetPT=[TrainedNetP(1:end-10) 'train' filesep];
UsedParams=getParamsStructFromFN(TrainedNetPT);
%%

LastP=[TFFolder FindNewestFolderWithPrefix(TFFolder,'Re','_train')];
[~, LastFN]=FindNewestFileWithPrefix(LastP,'Tra','mat');

[~, LastPNG]=FindNewestFileWithPrefix(LastP,'b','png');
NNres=double(imread(LastPNG));
NNres=reshape(NNres(1:128,256*5+(1:256),1),[128 128 MB]);

BartResN=A.RecTSAllX(:,:,[6 18],11);
BartResN=BartResN*grmss(NNres)/grmss(BartResN);


%% Now Multislice, mostly ok
WhichSlices=nSlicesNoMB:-1:1
%%
for SliI=WhichSlices
    CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];
    mkdir(CurOurPSli);
    
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,StartPoint:end);
    
    DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
    
    DataPCC=double(permute(DataCC,[3 2 1]));
    
    tmp=squeeze(DataPCC);
%     RealDataFacS(SliI)=grmss(AllData(5656:34:6767,:,:))/grmss(tmp);
    CurIDataV=Row(tmp)*100; %*RealDataFacS(SliI); % Change by Echo! WhichE
    CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
    
    Data=repmat(single(CurIDataVR),[BatchSizeInNN 1]);
    RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    save(RealDataFN,'Data');
end
disp('Save slices real data for NN');
% save([BaseOutLoc BaseFN filesep 'RealDataFacS.mat'],'RealDataFacS');
%% All reps
CurRealDataP=[mainP filesep 'RealData' filesep];
mkdir(CurRealDataP);

kx=BARTTrajCorr(1,:)*2*pi;
ky=BARTTrajCorr(2,:)*2*pi;
        
for SliI=1:12 %WhichSlices
    for r=1:nReps
        disp([SliI r]);
        nukData=ADataIsPy(:,:,SliI,r).';
        
        modx=double(exp(1i*(dx*kx+dy*ky))');
        modx(size(nukData,2),1)=0;
        
        nukData=circshift(nukData,-(Delay+2),2);
        nukData=nukData.*(modx.');
        
        % CAIPI Slice base
        nukData2(:,1:nTraj)=nukData(:,1:nTraj).*exp(-1i*cCAIPIVecY'.*(RotatedLocs(3,SliI)'));
        
        %nukData=nukData2(:,1:end-Delay-2);
        
        DataCC=CC(permute(nukData2,[3 2 1]),sccmtxS(:,1:nScc,SliI));
        
        CurIDataV=Row(DataCC)*200; % *RealDataFac;
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
        Data=repmat(single(CurIDataVR),[BatchSizeInNN 1]);
        RealDataFN=[CurRealDataP 'Sli' num2str(SliI) '_r' num2str(r,'%02d') '.mat'];
        %     RealDataFN=['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RealData/b_Ben14May_Sli5_r' num2str(r,'%02d') '.mat'];
        save(RealDataFN,'Data');
    end
end
disp('Saved slices real data for NN all reps');
%%
CurRealDataOutP=[mainP filesep 'RealDataOut' filesep];      
mkdir(CurRealDataOutP);
for SliI=WhichSlices
    mkdir([CurRealDataOutP 'Sli' num2str(SliI,'%02d')]);
end


%% Run on all repetitions (save all repetitions' data before)

% The network before trained and is on the following folder:
% TrainedNetP='/home/deni/TF/srez/RegridTry3C2_7TS_MID448_Sli06__2018-09-20_09-17-39_checkpoint';
TrainedNetP='//media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/RegridTry3C2_7TS_MID156_Sli02__2018-10-09_02-40-10_checkpoint';
TrainedNetP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/RegridTry3C2_7TS_MID156_Sli05__2018-10-07_21-06-36_checkpoint';
TrainedNetP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/Noise_recon/noise_0.005/RegridTry3C2_7TS_MID30__Sli06__2018-11-30_12-13-33_checkpoint';

TrainedNetPT=[TrainedNetP(1:end-10) 'train' filesep];
St=getParamsStructFromFN(TrainedNetPT);
St.BaseTSDataP=ParamsS.BaseTSDataP;
St.DatasetMatFN=ParamsS.DatasetMatFN;
St.BaseNUFTDataP=ParamsS.BaseNUFTDataP;
St.LoadAndRunOnData=1;
St.LoadAndRunOnData_checkpointP=TrainedNetP;
St.LoadAndRunOnData_Prefix=[mainP filesep 'RealData' filesep 'Sli' num2str(SliI) '_r'];
St.LoadAndRunOnData_OutP=[mainP filesep 'RealDataOut' filesep 'Sli' num2str(SliI,'%02d') filesep];

St.HowManyToRun=nReps;
Txt=gStruct2txt(St,'~/TF/Params.txt');

system('/home/deni/RunTFForMatlabx.sh');

    % disp('Prepared Params');                  
    %
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
%     system('sudo /home/deni/RunTFForMatlabx.sh');
% end
%% TF for all repetitions ALL SLICES 

%re-organize file paths 
ScanP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S01/RawData/';   
BaseFN='meas_MID32_gBP_ep2d_bold_multiecho_ASL_SMS_Spic_4min_FID33413';
MIDStr=BaseFN(6:11);
FN=[ScanP BaseFN '.dat'];

WhichE=1;
mainP=[ScanP BaseFN '_E' num2str(WhichE)];
%%
%baseP1= folder of the MLN recon (train/checkpoint)
%BaseP1='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S02/Recon/S02_Spi_r02_MLN_E2';
%BaseP1='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S02/Recon/S02_Spi_r01_MLN_E1';
BaseP1='/media/deni/bigdrive/SPIRAL_ASL/S04/Recon/S04_SPI_2mm_run1_skopeTraj_Sli1-10';
Prefix='RegridTry3C2_7TS_MID363_Sli';
%Prefix='RegridTry3C2_7TS_MID426_Sli';
%BaseP1='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S02/Recon_Noise/S02_spi_r2e1';
%Prefix='RegridTry3C2_7TS_MID433_Sli';

for SliI=1:12
    CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];
D=dir([BaseP1 filesep Prefix num2str(SliI,'%02d') '*']);
DNames={D.name}';
D=D(strhas(DNames,'check'));
TrainedNetP=[BaseP1 filesep D.name];

TrainedNetPT=[TrainedNetP(1:end-10) 'train' filesep];
St=getParamsStructFromFN(TrainedNetPT);
% St.BaseTSDataP=CurOurPSli;
St.BaseTSDataP='/media/deni/bigdrive/SPIRAL_ASL/S04/RawData/meas_MID363_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_FID46893_E1/Sli04/';
St.DatasetMatFN=ParamsS.DatasetMatFN;
St.BaseNUFTDataP=ParamsS.BaseNUFTDataP;

%
% TrainedNetP='/home/deni/TF/srez/RegridTry3C2_7TS_MID156_Sli12__2018-10-09_05-41-30_checkpoint';
% TrainedNetP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/RegridTry3C2_7TS_MID156_Sli06__2018-10-08_01-12-23_checkpoint';
% St=ParamsS;

St.LoadAndRunOnData=1;
St.LoadAndRunOnData_checkpointP=TrainedNetP;
St.LoadAndRunOnData_Prefix=[mainP filesep 'RealData' filesep 'Sli' num2str(SliI) '_r'];
St.LoadAndRunOnData_OutP=[mainP filesep 'RealDataOut' filesep 'Sli' num2str(SliI,'%02d') filesep];

St.HowManyToRun=nReps;
Txt=gStruct2txt(St,'~/TF/Params.txt');

system('/home/deni/RunTFForMatlabx.sh');
end 
%% Read NN output all reps
% NNOutP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnBP_9Aug18/meas_MID488_gBP_ep2d_bold_multiecho_ASL_SMS_Spic_4min_FID26056/RealDataOut/Sli06/';
NNOutP=St.LoadAndRunOnData_OutP;    
D=dir([NNOutP '*.mat']);
for r=1:numel(D)
    tmp=load([NNOutP D(r).name]);
    NNRec(:,:,r)=double(squeeze(tmp.x(1,:,:,1))+1i*squeeze(tmp.x(1,:,:,2)));
end
disp('ok')
fgmontage(NNRec(:,:,1))
%%
for SliI=1:12
    NNOutP=[mainP filesep 'RealDataOut' filesep 'Sli' num2str(SliI,'%02d') filesep];
    D=dir([NNOutP '*.mat']);
    for r=1:numel(D)
%         tmp=load([NNOutP D(r).name]);
        tmp=load([NNOutP 'OnRealData' num2str(r,'%02d') '.mat']);
        NNRecS(:,:,SliI,r)=double(squeeze(tmp.x(1,:,:,1))+1i*squeeze(tmp.x(1,:,:,2)));
    end
end
disp('ok')
NNRecS(1,1,nSlicesNoMB,1)=0;
NNRecSx=PartitionDim(NNRecS,2,2);
NNRecSx=CombineDims(NNRecSx,[5 3]);
fgmontage(NNRecSx(:,:,:,2))
%%
save([BaseP1 '/S04_SPI_2mm_r1e1skopeTraj_recon.mat'],'NNRecSx');
save('/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/Noise_recon/noise_0.01/S01_r1e1_no01_recon/S01_r1e1_no01_recon_2.mat','NNRecSx')
niftiN='/media/deni/bigdrive/SPIRAL_ASL/S04/Nifti_Untouched/S04_EPI_2mm_run1.nii';
Raw2Nii(abs(NNRecSx)*1000,[BaseP1 '/S04_SPI_2mm_r1e1_skopeTraj.nii'],'float32',niftiN);
Raw2Nii(fp_spi,[BaseP1 '/S04_SPI_2mm_r1e1_007_flipped.nii'],'float32',niftiN);

%% 
WhichE=1;
runN=1;
pathN='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S01/CleanNifti/';
niftiN='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S01/181012_DK_2DEPI_NIIs/S33_BP_FAIR_2x2_z3_RS_M0_E00_M_Echo_0_e11.nii';
%niftiN='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S02/Nifti/S02_ASL_epi_r01.nii';
Raw2Nii(abs(NNRecSx),[pathN  'S01_spi_r01_e1_noScale.nii'],'float32',niftiN);

Raw2Nii(echo2,[pathN  'S02_spi_sparse_r' num2str(runN,'%02d') '_e' num2str(WhichE) '.nii'],'float32',niftiN);
disp('ok nifti')


%%

Raw2Nii(abs(NNRecSx)*1000,'/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnDK_12Oct/MID30_E1_MLN.nii','float32');

SurSub=SurroundSubtraction(abs(NNRecSx)*1000);
Raw2Nii(SurSub,'/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/Echo1MLN_SurSub.nii','float32');
%%
NNRecSy=NNRecSx(:,:,[5 6 17 18],:);
Raw2Nii(abs(NNRecSy)*1000,'/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/Echo1MLN_4Sli_R2.nii','float32');

SurSuby=SurroundSubtraction(abs(NNRecSy)*1000);
Raw2Nii(SurSuby,'/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/Echo1MLN_SurSub_4Sli.nii','float32','/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/4Oct18_DK/Echo1.nii');
%%
MLNS=SmoothBySlices(abs(NNRecSx(:,:,:,:)),[20 20],2);
%%
MEven=ImResR(:,:,:,2:2:end);
MOdd=ImResR(:,:,:,1:2:end);

MEven=NNRec(:,:,:,2:2:end);
MOdd=NNRec(:,:,:,1:2:end);
Perf=abs(MEven)-abs(MOdd);
MPerf=mean(Perf,4);
SPerf=std(Perf,0,4);
tSNR_Perf=MPerf./SPerf;

TMsk1=tSNR_Perf>0.3;
TMsk1(:,[1:10 50:138 188:256])=0;
TwoD=Reshape4d22d(permute(Perf,[1 2 4 3]),TMsk1);
%%
Design=[zeros(1,18) repmat([ones(1,8) zeros(1,16)],[1 10])];
%%
RecTSAll=NaN([gsize(SensCC,1:2) MB nSlicesNoMB nReps]);
%% BART TS try: To do after having coil compression matrices, B0 correction matrices for all slices
% setenv('TOOLBOX_PATH','~/bart-0.4.03TS')
setenv('TOOLBOX_PATH','~/gUM/bart-0.4.03TS')
setenv('TOOLBOX_PATH','~/bart')
WhichRepToRun=10;
% WhichRepToRun=1:nReps;

WhichSliIsToRun=6;
% WhichSliIsToRun=1:nSlicesNoMB;
Delay=(WhichE-1)*6;
for RepToRun=WhichRepToRun    
    for SliI=WhichSliIsToRun
        SliIs=[SliI SliI+nSlicesNoMB];
        disp([RepToRun SliIs]);
%         nukData=ADataIsPy(:,:,SliI,RepToRun).';
%         nukData=nukData(:,3:end);
%         nukData=nukData(:,1:end-2);

%         BARTTrajCorr=BARTTrajAct./[1.01; 1.01; 1.0];
        
        nukData=ADataIsPy(:,:,SliI,RepToRun).';
        kx=BARTTrajCorr(1,:)*2*pi;
        ky=BARTTrajCorr(2,:)*2*pi;

        modx=double(exp(1i*(dx*kx+dy*ky))');
        modx(size(nukData,2),1)=0;

        nukData=circshift(nukData,-(Delay+2),2);
        nukData=nukData.*(modx.');
        nukData(:,1:nTraj)=nukData(:,1:nTraj).*exp(-1i*cCAIPIVecY.*(RotatedLocs(3,SliI)'));

        nukData=nukData(:,1:end-Delay-2);
%         nukData=nukData(:,(Delay+3):end);

        sccmtxBoth = sccmtxS(:,:,SliI);
        nScc=16;
        SensFCCBoth=squeeze(SensCCS(:,:,:,:,SliIs));
        
        DataCC=CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc));
        nukDataP=permute(DataCC,[3 1 4 2]);
        
        nBands=MB;
        TSBFA=permute(TSBS(:,:,SliIs),[4 1 6 5 3 2]);
        TSBFAm=repmat(TSBFA,[1 1 1 nScc 1 1 1 1 1]);
%         cCAIPIVecZ_CurSlis=exp(1i*cCAIPIVecY.*(RotatedLocs(3,SliIs)'));
        TSBFAm=TSBFAm.*permute(cCAIPIVecZ,[3 2 4 5 1 6]);
        
        TSBFm=TSBFAm;
        TSBFm=TSBFm/nTS;
        
        TSBFm=TSBFm(:,1:end-Delay,:,:,:,:);
        %writecfl('/tmp/TSB',TSBFm);
        %writecfl('/tmp/',TSB);

        
        TSCToUse=TSCS(:,:,:,SliIs);        
        SensP=permute(SensFCCBoth,[1 2 5 3 4]);
        SensW=SensP.*permute(TSCToUse(:,:,:,:),[1 2 5 6 4 3])*nTS;
        if(WhichE==2)
            SensW=SensW.*permute(exp(1i*2*pi*B0MF(:,:,SliIs)*TimeInMs2(end)/1000),[1 2 4 5 3 6]);
        end
        
        TrajW=repmat(BARTTrajCorr,[1 1 1 nScc 2 nTS]);
%         StartIdx=14392*2-45;
%         MTrajx=MTraj(1:2,StartIdx:5:(StartIdx+5*(nTraj-1)));
%         MTrajx(3,1)=0;
%         TrajW=repmat(MTrajx,[1 1 1 nScc 2 nTS]);
        
        disp(datestr(now));
        
        TrajW=TrajW(:,1:end-Delay,:,:,:,:);
        RecTS=bart(['pics -S -m -R W:3:0:' num2str(1e-4) ' -t'],TrajW, nukDataP, SensW);
       
        
         %    RecTS=bart(['pics -S -m -R T:3:3:' num2str(1e-5) ' -t'],TrajW, nukDataP, SensW);
%         RecTS=bart(['pics -S -m -R T:3:3:' num2str(1e-9) ' -R W:3:0:' num2str(1e-9) ' -t'],TrajW, nukDataP, SensW);
        disp(datestr(now));
        RecTSAll(:,:,:,SliI,RepToRun)=RecTS;
    end
end
% save([mainP filesep 'RecTSAll.mat'],'RecTSAll');
fgmontage(RecTS);title('BART TS')
%%
RecTSAllX=CombineDims(RecTSAll,[3 4]);
save([mainP filesep 'RecTSAllX.mat'],'RecTSAllX');
%%
RecTSAllXCurRep=CombineDims(RecTSAll(:,:,:,:,RepToRun),[3 4]);
save([mainP filesep 'RecTSAllXCurRep.mat'],'RecTSAllXCurRep');
%% End BART TS Try














%%
% pause(60*60*11);
for SliI=WhichSlices
    CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];

    ParamsSDefaults=getParamsStructFromFN('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/');
    ParamsS=ParamsSDefaults;
    ParamsS.SessionNameBase=['RegridTry3C2_7TS_' ScanP(end-3:end-1) '_Sli' num2str(SliI,'%02d')];
    ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    ParamsS.BaseTSDataP=CurOurPSli;
    ParamsS.BaseNUFTDataP=[mainP filesep];
    Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    
    system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
end

%% ESPIRIT
GOP2=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,double(TSC));
GOP_TS=ggpuNUFT_TS(GOP2,double(TSBF));
GOP_TS_ESP=ggpuNUFT_TS_MCnosum(GOP_TS,[1 1 32]);

CurSens=SensX(:,:,:,:,SliI);
weights=ones(gsize(CurSens,[1 2 4]));
ESP = ESPIRiT(CurSens,weights);

DataP2=double(permute(DataP,[3 2 1]));
disp('ok');
%%
XOP = Wavelet('Daubechies_TI',4,6);

splitWeight=0.5;
lam=1e-6;
% nIterSplit=15;
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, im_res*0, GOP_TS_ESP, ESP, 5,XFM,lam,splitWeight,nIterSplit);
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, weights*0, GOP_TS_ESP, ESP, 5,XOP,lam,splitWeight,nIterSplit);
% This works but we skiup and go directly to CC
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, weights*0, GOP_TS_ESP, ESP, 5,XOP,lam,splitWeight,2e-3);
disp('ok cgL1ESPIRiT');
%%
% save('CurStatusForL1ESPIRIT_MM_Ben4Min.mat')
%%
% resL1ESPIRiT1=resL1ESPIRiT(:,:,1);
% fgmontage(resL1ESPIRiT1,[0 7e-3])
% title('resL1ESPIRiT 2 maps, B0');
% ylabel(YLbl);
% xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);
% 
% gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr ],[]) 
% close(gcf);
% save([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '.mat'],'resL1ESPIRiT1');
% disp(['Saved ' ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '.mat']);
%%
weights=ones(gsize(CurSens,[1 2 4]));
ESP = ESPIRiT(CurSens,weights);

weightsCC=ones(gsize(CurSens,[1 2 4]));
ESPCC = ESPIRiT(SensCC,weightsCC);

GOP_TS_ESPCC=ggpuNUFT_TS_MCnosum(GOP_TS,[1 1 ncc]);

% GOP_MC_CC = ggpuNUFT_TS_MCx(BARTTrajAct,Sz2,osf,wg,sw,TSB.',TSC,SensCC(:,:,:,1));


    
DataCC=(sccmtx(:,1:ncc).'*nukData);

DataPCC=double(permute(DataCC,[3 2 1]));

XOP = Wavelet('Daubechies_TI',4,6);

splitWeight=0.5;
lam=1e-7;
nIterSplit=2e-3;
% nIterSplit=15;
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, im_res*0, GOP_TS_ESP, ESP, 5,XFM,lam,splitWeight,nIterSplit);
[resL1ESPIRiTCC] = cgL1ESPIRiT(DataPCC, weightsCC*0, GOP_TS_ESPCC, ESPCC, 5,XOP,lam,splitWeight,nIterSplit);

disp('ok cgL1ESPIRiTCC');
%%
resL1ESPIRiTCC1=resL1ESPIRiTCC(:,:,1);
% fgmontage(cat(3,resL1ESPIRiT1, resL1ESPIRiTCC1),[0 7e-3])
% title(['resL1ESPIRiT 2 maps, B0. Right - with CC -> ' num2str(ncc)]);
fgmontage(resL1ESPIRiTCC1,[0 7e-3])
title(['resL1ESPIRiT ' num2str(nMapsToUse) ' maps, B0  with CC -> ' num2str(ncc)]);
ylabel(YLbl);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);
%%
gprint(get(gcf,'Number'),[mainP filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC' num2str(ncc)],[]) 
close(gcf);
save([mainP filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat'],'resL1ESPIRiTCC1');
disp(['Saved ' mainP filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat']);
%%
CurBartTraj=BARTTrajAct;
% CurBartTraj=BARTTrajMS;
nTraj=size(CurBartTraj,2);
% figure;plot(CurBartTraj(1,:),CurBartTraj(2,:),'.');
kMax=ceil(max(max(abs(CurBartTraj),[],2)));
% xlabel(nTraj);
Acc=(kMax*2)^2/nTraj;
% title(['kMax: ' num2str(kMax) ' Acc (for kMax): ' num2str(Acc)]);
%%
% % TBaseP='~/HomeA/TF/srez/';
% TBaseP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';
% 
% DataH=size(NMap,1);
% DataW=size(NMap,2);
% DataCh=size(NMap,3)*nChToUseInNN*2;
% LabelsH=Sz2(1);
% LabelsW=Sz2(2);
% LabelsCh=2;
% 
% ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
%   'dataset',TFDataP,'learning_rate_start',0.002,...
%   'learning_rate_half_life',30,... % in minutes if <1000
%   'summary_period',0.5,'checkpoint_period',20,...
%   'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
%   'WL1_Lambda',0,'WL2_Lambda',0,...
%   'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'DiscStartMinute',500,...
%   'ShowRealData',1,'CmplxBias',0,...
%   'InputMode','RegridTry1',...
%   'NetMode','RegridTry1C2_TS',...
%   'SessionNameBase','RegridTry1C2_TS',...
%   'ImgMode','Cmplx',...
%   'nTimeSegments',7,...
%   'UseSharedWightesInRelaxedFT',1,...
%   'WPhaseOnly',0.001,...
%   'NMAP_FN',NMapFN,...
%   'RealDataFN',RealDataFN);
% 
% ParamsS=ParamsSDefaults;
% Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
% %%
% ParamsS=ParamsSDefaults;
% ParamsS.WPhaseOnly=0.01;
% ParamsS.train_time=120;
% Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
% system('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
% 
% system('sudo bash -c ''echo "I am $USER, with uid $UID"''')
% system('sudo -H -u a bash -c ''echo "I am $USER, with uid $UID"''')
% 
% sendTFMail(TBaseP,ParamsS,Txt);
% %%
% ParamsS=ParamsSDefaults;
% ParamsS.MapSize=5;
% ParamsS.train_time=120;
% ParamsS.dataset='dataFaceP4C';
% ParamsS.learning_rate_start=0.001;
% ParamsS.learning_rate_half_life=60;
% ParamsS.WL1_Lambda=0;
% ParamsS.WL2_Lambda=0;
% Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
% system('/home/a/RunTFForMatlab.sh');
% sendTFMail(TBaseP,ParamsS,Txt);









%% Now Multislice, mostly ok
WhichSlices=nSlicesNoMB:-1:1
%%
for SliI=WhichSlices
    CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];
    mkdir(CurOurPSli);
    
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,StartPoint:end);
    
    DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
    
    DataPCC=double(permute(DataCC,[3 2 1]));
    
    tmp=squeeze(DataPCC);
%     RealDataFacS(SliI)=grmss(AllData(5656:34:6767,:,:))/grmss(tmp);
    CurIDataV=Row(tmp)*100; %*RealDataFacS(SliI); % Change by Echo! WhichE
    CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
    
    Data=repmat(single(CurIDataVR),[12 1]);
    RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    save(RealDataFN,'Data');
end
disp('Save slices real data for NN');
% save([BaseOutLoc BaseFN filesep 'RealDataFacS.mat'],'RealDataFacS');
%% All reps
CurRealDataP=[mainP filesep 'RealData' filesep];
mkdir(CurRealDataP);

kx=BARTTrajCorr(1,:)*2*pi;
ky=BARTTrajCorr(2,:)*2*pi;
        
for SliI=1:12 %WhichSlices
    for r=261:nReps
        disp([SliI r]);
        nukData=ADataIsPy(:,:,SliI,r).';
        
        modx=double(exp(1i*(dx*kx+dy*ky))');
        modx(size(nukData,2),1)=0;
        
        nukData=circshift(nukData,-(Delay+2),2);
        nukData=nukData.*(modx.');
        
        % CAIPI Slice base
        nukData(:,1:nTraj)=nukData(:,1:nTraj).*exp(-1i*cCAIPIVecY.*(RotatedLocs(3,SliI)'));
        
        nukData=nukData(:,1:end-Delay-2);
        
        DataCC=CC(permute(nukData,[3 2 1]),sccmtxS(:,1:nScc,SliI));
        
        CurIDataV=Row(DataCC)*200; % *RealDataFac;
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
        Data=repmat(single(CurIDataVR),[BatchSizeInNN 1]);
        RealDataFN=[CurRealDataP 'Sli' num2str(SliI) '_r' num2str(r,'%02d') '.mat'];
        %     RealDataFN=['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RealData/b_Ben14May_Sli5_r' num2str(r,'%02d') '.mat'];
        save(RealDataFN,'Data');
    end
end
disp('Saved slices real data for NN all reps');
%%
CurRealDataOutP=[mainP filesep 'RealDataOut' filesep];      
mkdir(CurRealDataOutP);
for SliI=WhichSlices
    mkdir([CurRealDataOutP 'Sli' num2str(SliI,'%02d')]);
end
%%
WhichSlices=nSlices:-1:1;
%%
clear TSBS TSCS
%%
nTS=11;
for SliI=WhichSlices
    B0M2=B0M2S(:,:,SliI);
    Mgc=MgcS(:,:,SliI);
    
    Mskc=Mgc>7e-5;
    MskcE=imdilate(imfill(Mskc,'holes'),strel('disk',5,8));
    B0M2(~MskcE)=0;
    
    AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);

   % AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
    MgcN=Mgc./grmss(Mgc);
    WE=Col(MgcN);
    % WE=Col(MgcN)*0+1;
    WeightedE=WE.*E;
    
    % Fessler time segmentation
    clear ErrTS
    TS_Thresh=1e-5;
        FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
        TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
        TSC2=reshape(TSC,prod(Sz2),nTS);
        WTSC2=WE.*TSC2;
        TSB=(WeightedE.')/(WTSC2.');% W in both sides

        CurErrTS=grmss(WeightedE-WTSC2*(TSB.')); %include W
        ErrTSS(nTS,SliI)=CurErrTS;

        disp([datestr(now) ' Sli #' num2str(SliI) ' nTS ' num2str(nTS) ' err=' num2str(CurErrTS)]);
    TSBS(:,:,SliI)=TSB;
    TSCS(:,:,:,SliI)=TSC;
end
disp('ok TSBS TSCS');
%%
% load([ScanP BaseFN filesep 'TSBS_TSCS.mat'])
save([mainP filesep 'TSBS_TSCS.mat'],'TSBS','TSCS','ErrTSS');
disp('Saved TSBS TSCS');
%%





%% save TSBF TSC Per slice
for SliI=WhichSlices
    YLbl=['Sli' num2str(SliI,'%02d')];
    
    SliIP=[mainP filesep YLbl filesep];
    mkdir(SliIP);
    
    TSB=squeeze(TSBS(:,:,SliI));
    TSC=squeeze(TSCS(:,:,:,SliI));
    TSBF=TSB(IB_TimeInMs2,:).';
    Sens=squeeze(SensX(:,:,:,1,SliI));
    
    save([SliIP 'TSBF_TSC_Sens.mat'],'TSBF','TSC','Sens');
    disp(['Saved ' SliIP 'TSBF_TSC_Sens.mat']);
end
%% Clean save per slice
% BaseOutLoc='/media/a/DATA/';
for SliI=WhichSlices
    SensCC=squeeze(SensCCS(:,:,:,:,SliI));
    SensMsk=single(grmss(SensCC(:,:,:),3)>0.01);
    CurOurP=[mainP filesep];
    CurOurPSli=[mainP filesep 'Sli' num2str(SliI,'%02d') filesep];
%     save([CurOurPSli 'Sens.mat'],'CurSens');
    save([CurOurPSli 'SensCC.mat'],'SensCC','sccmtx','SensMsk');
    SensCC=squeeze(SensCC(:,:,:,1));
    save([CurOurPSli 'SensCC1.mat'],'SensCC','sccmtx','SensMsk');
    
    TSB=squeeze(TSBS(:,:,SliI));
    TSC=squeeze(TSCS(:,:,:,SliI));
    TSBF=TSB(IB_TimeInMs2,:).';
    Mgc=MgcS(:,:,SliI);

    save([CurOurPSli 'B0TS.mat'],'TSBF','TSB','TSC','osf','wg','sw','Mgc','TimeInMs2');
    disp(['Saved ' num2str(SliI)]);
end
%% Call TF
% pause(60*60*11);
for SliI=WhichSlices
    CurOurPSli=[ScanP BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];

%     ParamsSDefaults=getParamsStructFromFN('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/');
    ParamsSDefaults=getParamsStructFromFN('/home/deni/TF/srez/RegridTry3C2_7TS_S02_Sli01__2018-07-23_09-40-59_train/');
    ParamsS=ParamsSDefaults;
    ParamsS.SessionNameBase=['RegridTry3C2_7TS_' ScanP(end-3:end-1) '_Sli' num2str(SliI,'%02d')];
    ParamsS.SessionNameBase=['RegridTry3C2_7TS_' MIDStr '_Sli' num2str(SliI,'%02d')];
    
    ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    ParamsS.BaseTSDataP=CurOurPSli;
    ParamsS.BaseNUFTDataP=[ScanP BaseFN filesep];
%     Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    Txt=gStruct2txt(ParamsS,'~/TF/Params.txt');
    
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
    system('sudo /home/deni/RunTFForMatlabx.sh');
end
% /home/deni/TF/srez/RegridTry3C2_7TS_g18_Sli06__2018-08-16_12-35-57_train
%% GPU TS
for SliI=WhichSlices
%     TSBF=squeeze(TSBS(IB_TimeInMs2,:,SliI)).';
%     Sens=squeeze(SensX(:,:,:,1,SliI));
    TSBF=squeeze(TSBS(:,:,SliI)).';
    Sens=squeeze(SensCCS(:,:,:,:,SliI));
    
    GOP_MCSMBMSS{SliI} = ggpuNUFT_TS_MC_MB_MS(BARTTrajAct,NTrg1,osf,wg,sw,TSBF,TSCS(:,:,:,SliI),Sens);
    GOP_MCS{SliI}=GOP_MCSMBMSS{SliI};
end
disp('ok');
%%
nfnlCgIters=40;
RunFnlViewAmp=1;

TVW=1e-6;

XFMStr='Daubechies';
filterSize=4;
wavScale=4;

FigH=4000;
figure(FigH);clf;

Delay=(WhichE-1)*6;

for SliI=WhichSlices
    disp(SliI);
    
    BARTTrajCorr=BARTTrajAct./[1.01; 1.01; 1.0];
    
    nukData=ADataIsPy1(:,:,SliI,RepToRun).';
    kx=BARTTrajCorr(1,:)*2*pi;
    ky=BARTTrajCorr(2,:)*2*pi;
    
    modx=double(exp(1i*(dx*kx+dy*ky))');
    modx(size(nukData,2),1)=0;
    
    nukData=circshift(nukData,-(Delay+2),2);
    nukData=nukData.*(modx.');
    nukData=nukData(:,1:end-Delay-2);
    %         nukData=nukData(:,(Delay+3):end);
    
    sccmtxBoth = sccmtxS(:,:,SliI);
    nScc=16;
%     SensFCCBoth=squeeze(SensCCS(:,:,:,:,SliIs));
    
    DataCC=CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc));
    nukDataP=permute(DataCC,[2 1]);
    DataP=nukDataP;
        
%     nukData=ADataIsPy(:,:,SliI,1).';
%     nukData=nukData(:,3:end);
%     % nukData=nukData.*T2SEstComp;
%     
%     DataP=nukData;

    AOdd = GOP_MCSMBMSS{SliI};
    % AOdd = GOP_MC;


    if(strcmp(XFMStr,'None'))
        XFM=1;
        xfmWeight=0;
    else
        XFM = Wavelet(XFMStr,filterSize,wavScale);
        xfmWeight = 1e-6;	% Weight for Transform L1 penalty
    end


    param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);

    param.data =     DataP;

    res=zeros(Sz2);
    if(nShots>1)
        res=repmat(resA,[1 1 1 nShots]);
        res=res+randn(size(res))*max(abs(resA(:)))/20;
    end

    if(~isfield(param,'ShowFig'))
        param.ShowFig=true;
    end
    StartTime_fnl=now;
    param.Verbose=false;
    clear ObjConv Score
    
    clf;
    
    for n=1:nfnlCgIters
        [res, CurObj] = fnlCg(res,param);
        ObjConv(n)=CurObj;
        im_res = param.XFM'*res;
        if(param.ShowFig)
%             figure(FigH); 
            subplot(1,3,1);
            gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
            cx=caxis;
            caxis(cx/RunFnlViewAmp);
            subplot(1,3,2);
            gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
            subplot(1,3,3);
            plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
        end
        if(n>1)
            dObjP=(ObjConv(n-1)-CurObj)/ObjConv(n-1);
            disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
            if(dObjP<1e-3)
                disp('Not advancing. Stopping.');
                break;
            end
        else
            disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
        end

        if(nShots==1)
            resA=res;
        end
    end
    
    im_resS(:,:,SliI)=im_res;
end
disp('Finished im_resS');
%%
fgmontage(im_resS)
title(['SparseMRI With TS 1map no CC']);

XFMStrFull=['[' XFMStr ',' num2str(filterSize) ',' num2str(wavScale) ',' num2str(xfmWeight) ']'];
XLbl=['L' num2str(param.pNorm) ',TVW=' num2str(param.TVWeight) ',' XFMStrFull ];
xlabel(XLbl)

gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'SparseMRI_TS_1map_noCC_' XFMStrFull],[]) 
close(gcf);
%%
CurOurP=[ScanP BaseFN filesep];
mkdir(CurOurP);
save([CurOurP 'im_resS.mat'],'im_resS');
%% L1cgESPIRIT_CC on all:
% nMapsToUse=size(SensCC,4);
nMapsToUse=1;
for SliI=WhichSlices
    SensCC=squeeze(SensCCS(:,:,:,1:nMapsToUse,SliI));    
    TSBF=squeeze(TSBS(IB_TimeInMs2,:,SliI)).';
    GOP2S{SliI}=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(size(BARTTrajAct,2),1),osf,wg,sw,Sz2,double(TSCS(:,:,:,SliI)));
    GOP_TSS{SliI}=ggpuNUFT_TS(GOP2S{SliI},double(TSBF));
    GOP_TS_ESPCCS{SliI}=ggpuNUFT_TS_MCnosum(GOP_TSS{SliI},[1 1 ncc]);
    
    weightsCC=ones(gsize(SensCC,[1 2 4]));
    ESPCCS{SliI} = ESPIRiT(SensCC,weightsCC);
end
disp('ok');
%%
clear resL1ESPIRiTCCS
%%
for SliI=WhichSlices
    disp(['Sli #' num2str(SliI) ' ' datestr(now)]);
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,StartPoint:end);
    DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
    DataPCC=double(permute(DataCC,[3 2 1]));

    XOP = Wavelet('Daubechies_TI',4,6);
    splitWeight=0.5;
    lam=1e-7;
    nIterSplit=5e-4;
    % nIterSplit=15;
    [resL1ESPIRiTCCS(:,:,:,SliI)] = cgL1ESPIRiT(DataPCC, weightsCC*0, GOP_TS_ESPCCS{SliI}, ESPCCS{SliI}, 5,XOP,lam,splitWeight,nIterSplit);
end
disp('ok resL1ESPIRiTCCS');
%%
resL1ESPIRiTCCS1=squeeze(resL1ESPIRiTCCS(:,:,1,:));
fgmontage(resL1ESPIRiTCCS1,[0 7e-3])
title(['resL1ESPIRiT ' num2str(nMapsToUse) ' maps, B0, with CC -> ' num2str(ncc)]);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC' num2str(ncc)],[]) 
close(gcf);
save([ScanP BaseFN filesep 'L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat'],'resL1ESPIRiTCCS1');
disp(['Saved ' ScanP BaseFN filesep 'L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat']);
%% All reps
resL1ESPIRiTCCSA=zeros([size(weightsCC) nSlices nReps]);
%%
for SliI=WhichSlices
    for r=1:nReps
        disp(['Sli #' num2str(SliI) ' rep #' num2str(r) ' ' datestr(now)]);
        nukData=ADataIsPy(:,:,SliI,r).';
        nukData=nukData(:,StartPoint:end);
        DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
        DataPCC=double(permute(DataCC,[3 2 1]));
        
        XOP = Wavelet('Daubechies_TI',4,6);
        splitWeight=0.5;
        lam=1e-6;
        nIterSplit=5e-4;
        % nIterSplit=15;
        [resL1ESPIRiTCCSA(:,:,:,SliI,r)] = cgL1ESPIRiT(DataPCC, weightsCC*0, GOP_TS_ESPCCS{SliI}, ESPCCS{SliI}, 5,XOP,lam,splitWeight,nIterSplit);
    end
end
disp('ok resL1ESPIRiTCCSA');
%%
resL1ESPIRiTCCS1A=squeeze(resL1ESPIRiTCCSA(:,:,1,:,:));
save([ScanP BaseFN filesep 'L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC_Allreps.mat'],'resL1ESPIRiTCCS1A');
disp(['Saved ' ScanP BaseFN filesep 'L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC_Allreps.mat']);
%%








%% Call TF to get results
DBase=dir([ScanP BaseFN filesep 'RegridTry3C2_7TS*_train']);
for SliI=WhichSlices
    St=getParamsStructFromFN([ScanP BaseFN filesep DBase(1).name filesep]);
    St.LoadAndRunOnData=1;
    D=dir([ScanP BaseFN filesep 'RegridTry3C2_7TS_' ScanP(end-3:end-1) '_Sli' num2str(SliI,'%02d')  '*checkpoint*']);
    St.LoadAndRunOnData_checkpointP=[ScanP BaseFN filesep D.name];
    St.LoadAndRunOnData_Prefix=[ScanP BaseFN filesep 'RealData' filesep 'Sli' num2str(SliI) '_r'];
    St.LoadAndRunOnData_OutP=[ScanP BaseFN filesep 'RealDataOut' filesep 'Sli' num2str(SliI,'%02d') filesep];
    Txt=gStruct2txt(St,'~/HomeA/TF/Params.txt');
    % disp('Prepared Params');
    %
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
    system('sudo /home/deni/RunTFForMatlabx.sh');
end


%% Run already trained network on real data all reps
TrainNNRunP='/home/deni/TF/srez/RegridTry3C2_7TS_g18_Sli06__2018-08-16_12-35-57_train/';

St=getParamsStructFromFN(TrainNNRunP);
St.LoadAndRunOnData=1;
St.LoadAndRunOnData_checkpointP=[TrainNNRunP(1:end-6) 'checkpoint'];
St.LoadAndRunOnData_Prefix=[ScanP BaseFN filesep 'RealData' filesep 'Sli' num2str(SliI) '_r'];
St.LoadAndRunOnData_OutP=[ScanP BaseFN filesep 'RealDataOut' filesep 'Sli' num2str(SliI,'%02d') filesep];
% Txt=gStruct2txt(St,'~/HomeA/TF/Params.txt');
Txt=gStruct2txt(St,'~/TF/Params.txt');
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
system('sudo /home/deni/RunTFForMatlabx.sh');
%% Read NN output all reps
NNOutP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnBP_9Aug18/meas_MID488_gBP_ep2d_bold_multiecho_ASL_SMS_Spic_4min_FID26056/RealDataOut/Sli06/';
D=dir([NNOutP '*.mat']);
for r=1:numel(D)
    tmp=load([NNOutP D(r).name]);
    NNRec(:,:,r)=double(squeeze(tmp.x(1,:,:,1))+1i*squeeze(tmp.x(1,:,:,2)));
end
%%
MEven=NNRec(:,:,2:2:end);
MOdd=NNRec(:,:,1:2:end);
Perf=abs(MEven)-abs(MOdd);
MPerf=mean(Perf,3);
SPerf=std(Perf,0,3);
tSNR_Perf=MPerf./SPerf;
%%
RecP='/media/a/DATA/Sep19_OnDeni/RegridTry3C2_7TS_MID543_Sli06__2018-09-29_11-22-06_train/';
[ScrN,BatchN,MinN,LastFN]=GraphOptFromFolderf(RecP,false);
X=double(rgb2gray(imread([RecP LastFN])));
X=X(1:128,256*5+(1:256));
fgmontage(X)
%%