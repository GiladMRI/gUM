ScanP='/autofs/cluster/kawin/Gilad/gep2d_Try1_Phantom_Traj1/';
BaseFN='meas_MID00394_FID40730_gep2d_gese_sms_mgh_EPTI';
RefFldMapP=[ScanP 'meas_MID04675_FID21637_gre_te4_9' filesep];

ScanP='/autofs/cluster/kawin/Gilad/gep2d_Try1_Phantom_Traj1/';
BaseFN='meas_MID01330_FID41644_SEPTI';

ScanP='/autofs/cluster/kawin/Gilad/gep_Phantom/';
BaseFN='meas_MID00924_FID42864_gepBase_T2_50rep';
BaseFN='meas_MID00912_FID42852_gepBase_12_50rep';
% BaseFN='meas_MID00904_FID42844_gepBase_8cor_50rep';
BaseFN='meas_MID00914_FID42854_gepBase_13_50rep';

ScanP='/autofs/cluster/kawin/Gilad/gep11mm/';
BaseFN='meas_MID01551_FID43495_gepBase_13_50rep';

ScanP='/autofs/space/daisy_001/users/data/Gilad/gep_CL/';
RefFldMapP=[ScanP 'meas_MID01928_FID43869_gre_4echo_24_26_G2' filesep];
BaseFN='meas_MID01934_FID43875_gepBase_T1_50rep';
BaseFN='meas_MID01936_FID43877_gepBase_T2_50rep';
BaseFN='meas_MID01938_FID43879_gepBase_T3_50rep';
BaseFN='meas_MID01940_FID43881_gepBase_T4_50rep';
% BaseFN='meas_MID01942_FID43883_gepBase_T5_50rep';
% BaseFN='meas_MID01960_FID43901_gepBase_T6cor_50rep';
% BaseFN='meas_MID01962_FID43903_gepBase_T7cor_50rep';
% BaseFN='meas_MID01950_FID43891_gepBase_12_50rep';
% BaseFN='meas_MID01952_FID43893_gepBase_13_50rep';
% BaseFN='meas_MID01970_FID43911_gepBase_T2_50rep_ContinuousSlices';


ScanP='/autofs/space/daisy_001/users/data/Gilad/gep_GL/';
BaseFN='meas_MID03807_FID45716_gepBase_T1_50rep';

ScanP='/autofs/space/daisy_001/users/data/Gilad/gEPTI_Try1/';
BaseFN='meas_MID00491_FID49056_gepBase_T1_50rep_continuousSlices';

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
disp('ok 2');
%% Read raw
ADataFN=[mainP filesep 'AData.mat'];
if(exist(ADataFN,'file'))
    load(ADataFN);
else
    AData = mapVBVD(FN);
    save(ADataFN,'AData');
end
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
ParamA=WipMemBlock.alFree{41};

TrajType=mod(ParamA,100);

TE0_ms=2.38;

G2mm=load('Grads2mmb.mat');
% QQ=load('GAll1p9mmVD1PAT3Pause.mat');
QQ=load('GAll1p9mmVD1PAT3Pauseb.mat');
AllGrads=cat(1,G2mm.Grads{:}).';
AllGrads=cat(2,AllGrads,QQ.GAll(:,7:8));
G11mm=load('Grads11mmb.mat');
AllGrads11=cat(1,G11mm.Grads{:}).';
AllGrads=cat(2,AllGrads,AllGrads11);

nInnerShotsAll=[8 6 4 3 2 8 6 4 3 2 3 3 6 4 3 2];
nInnerShots=nInnerShotsAll(TrajType);
% GTrajaCBase=G2mm.Grads{TrajType};
GTrajaCBase=AllGrads(:,TrajType).';
% load('GTraj16.mat');
% TrajType=WipMemBlock.adFree{12};
% GTrajaCBase=GTraj16;
% ResType=floor((TrajType-10)/2)+1; % 1.9,1.3
% TimingType=mod(TrajType-10,2)+1; % 6ms, 8ms
% load('GAll68.mat');
% GNav=load('GNav1ms.mat');
% GNav=GNav.GNav;
% GTrajaCBase=GAll(:,TimingType,ResType);
% 
% if(TimingType==1)
%     nInnerShots=8;
% else
%     nInnerShots=6;
% end
%%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;
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

% FOVx=220;

dFOV=FOVx/1000;

FOV_mm=FOVx;

disp('ok');

GradReduceFac=1; % WipMemBlock.adFree{4};
PhiRotPerRep=110; % WipMemBlock.adFree{8};
AcqDwellTime_us=1; %WipMemBlock.adFree{13}/1000;

dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;
%%
RLocs=load([RefFldMapP 'Locs.mat']);
sTwix=load([RefFldMapP 'sTwixX.mat']);
RefLocs=RLocs.RotatedLocs;
CurzLocs=RotatedLocs(3,:);
RefzLocs=RefLocs(3,:);

RefMaps=load([RefFldMapP 'B0T2S.mat']);
RefSens=load([RefFldMapP 'Sens.mat']);

RefSV1=grmss(RefMaps.s_vals,3);

% RefSV1x=padarray(RefSV1(:,:,Ord(1)),[0 3],'Both');
% RefSV1x=padarray(RefSV1(:,:,Ord(1)),[0 6],'Post');
% RefSV1x=padarray(padarray(RefSV1(:,:,Ord(1)),[0 4],'Post'),[0 2],'Pre');
% RefSV1x=padarray(padarray(RefSV1(:,:,Ord(1)),[0 5],'Post'),[0 1],'Pre');
RefSV1x=padarray(padarray(RefSV1(:,:,:),[0 5],'Post'),[0 1],'Pre');
RefSV1x=rot90(RefSV1x);
RefSV1x=imresize(RefSV1x,Sz);

% figure;imshowpair(grmss(Rec1p,8)/grmss(Rec1p),RefSV1x/grmss(RefSV1x));
% figure;imshowpair(grmss(Rec1p,8)/gmax(grmss(Rec1p,8)),RefSV1x/gmax(RefSV1x),'checkerboard');

RefSensx=padarray(padarray(RefSens.SensB(:,:,:,:,1),[0 5],'Post'),[0 1],'Pre');
RefSensx=rot90(RefSensx);
RefSensx=imresize(RefSensx,Sz); % X Y Ch SliI

% RefB0x=padarray(padarray(RefMaps.B0M_Hz,[0 5],'Post'),[0 1],'Pre');
RefB0x=padarray(padarray(RefMaps.UpdatedB0Map_Hz,[0 5],'Post'),[0 1],'Pre');
RefB0x=rot90(RefB0x);
RefB0x=imresize(RefB0x,Sz); % X Y Ch SliI
%%
GTrajaC=GTrajaCBase/GradReduceFac;
% GNavC=GNav/GradReduceFac;

for j=1:5
%     GTrajaCM(:,j)=[zeros((j-1)*50,1); GTrajaC; zeros((5-j)*50,1); GNavC; -GNavC];
    GTrajaCM(:,j)=GTrajaC;
end

for j=1:5
    g=GTrajaCM(:,j);
    k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
    s=diff(g)/GradDwellTime_ms;
    
    kK=k*FOV_mm/1000/2/pi;
    kkM(:,j)=kK;
end

Kmax=ceil(max(abs(kK)));
NTrg=Kmax*2;

HighRes=false;
if(NTrg<150)
    disp('High res mode');
    disp(['Kmax is ' num2str(Kmax) ' aimed for ' num2str(NTrg) ' but setting to 120']);
    NTrg=120;
else
    HighRes=true;
    disp('High res mode');
    disp(['Kmax is ' num2str(Kmax) ' aimed for ' num2str(NTrg) ' but setting to 216']);
    NTrg=216;
end

res_mm=FOVx/(max(abs(kK))*2);
% figure;plot(kK);axis square;axis equal;ylabel('k');
%
TrgSz=[NTrg NTrg];
OnesSens=ones(TrgSz);
%%
Extra=4;
kKA=[zeros(4,5); kkM; repmat(kkM(end,:),[10 1])]; % in us

% Delay_us=-2.5; % positive data rightward %-6; % clockwise
% Delays_us=-2:0.1:4;

Delays_us=1.6;
dlyI=1;
% for dlyI=1:numel(Delays_us)
    disp('-------');
    disp(dlyI);
    Delay_us=Delays_us(dlyI);
    AcqTimePoints_us=Extra*GradDwellTime_us+Delay_us+(0:AcqDwellTime_us:48000);

    nRepsToUse=nRepsHdr;
    
    clear TrajM
    for i=1:nRepsToUse
        CurJitterIdx=mod(i-1,5)+1;
        Traj=interp1((0:(size(kKA,1)-1))*GradDwellTime_us,kKA(:,CurJitterIdx),AcqTimePoints_us);
        TrajM(i,:)=Traj.*exp(-1i*2*pi*PhiRotPerRep/360*(i-1));
    end
    
    AcqLen_us=AcqTimePoints_us(end)-AcqTimePoints_us(1)+AcqDwellTime_us;
    
    BARTTraj=cat(3,real(TrajM),imag(TrajM),imag(TrajM)*0);
    BARTTrajP=permute(BARTTraj,[3 2 1]);

    kx=BARTTrajP(1,:,:)*2*pi;
    ky=BARTTrajP(2,:,:)*2*pi;

    modx=exp(-1i.*(dx*kx+dy*ky));
    disp('ok mod');
    
    nAcqPoints=numel(Traj);

Sz=TrgSz;
Sz16=FillOnesTo16(Sz);
%%
% image_datan = AData{1}.noise.unsorted(); % no slicing supported atm
% image_data = ADatax.image.unsorted(); % no slicing supported atm
% XX=PartitionDim(PartitionDim(image_data,3,nRepsHdr*nSlices),4,nRepsHdr);
% XX=CombineDims(XX,[3 1]);
%%
SlicesToRead=1;
clear KernsP_TSTHLRBase KernsPMMedBase
for SlicesToRead=4:nSlices
    disp(['Slice ' num2str(SlicesToRead)]);
RepsToRead=1:(nRepsHdr);
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepsToRead,:,:,:,:,:,:,:,:);
ADataIsL=squeeze(ADataIsL);
ADataIsL=CombineDims(ADataIsL,[4 1]);
disp('Read data');
%%
% % ChRMS=grmss(ADataIsL,[1 3]);
% % [SChRMS ChOrd]=sort(ChRMS,'descend');
% % 
% % Ch2D=CombineDims(ADataIsL,[3 1]);
% % [~,S,sccmtx] = svd(Ch2D(1:10:end,:),'econ');
% % clear Ch2D 
% % 
% % ncc=31;
% % ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
% % for i=1:ncc
% %     ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
% % end
% % DataC=permute(ADataIsLCC,[1 3 2]);
% % disp('ok cc');
% % 
% % sccmtxS(:,:,SlicesToRead)=sccmtx;
%% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
% TrajPartToUse=0+(1:9000);
TrajPartToUse=0+(1:2000);
RepsToUse=1:(nRepsHdr-1);

RepsToUse=5:38;

% % DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 3 5 6 7 8 2]).*modx;
% DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx;
% OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 size(DataPC,8)]);
% 
% Rec1p=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
%%
% Rots=mod(360-mod(Ord*50-50,360),360);
Rots=mod(360-mod((Ord-1)*50,360),360);
CurRot=Rots(SlicesToRead);
RotsX=mod(110*(1:36),360);
ShiftIdx=find(CurRot==RotsX);
RepsI=mod( (1:36) +ShiftIdx -1,36)+1;
% RepsI=mod( (1:36) +8 -1,36)+1;
% RepsI=mod( (1:36) +22 -1,36)+1;
% DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsI,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,RepsI);
DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsI,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,1:36);
Rec1p=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsI),(DataPC(:,TrajPartToUse,RepsI,1,1,1,1,:)),OnesSensC);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
%%

Rec1pRMS=grmss(Rec1p,8);
Rec1pRMSS(:,:,SlicesToRead)=Rec1pRMS;
Rec1pS(:,:,:,SlicesToRead)=squeeze(Rec1p);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
end

%%
% Raw2Nii('/autofs/space/daisy_002/users/Gilad/gUM/AA.nii',Rec1pRMSS,');
Raw2Nii(Rec1pRMSS(:,:,ROrd),'/autofs/space/daisy_002/users/Gilad/gUM/AA.nii','float32');
% /autofs/space/daisy_002/users/Gilad/mricron_lx
% system('/autofs/space/daisy_002/users/Gilad/mricron_lx /autofs/space/daisy_002/users/Gilad/gUM/AA.nii');
% title([num2str(numel(RepsToUse)) ' shots data x ' num2str(numel(TrajPartToUse)) 'points, pre channel recon, RMS']);
%%

[Out B BN]=CalcSlicesSNR(Rec1pRMS,false,7);
Msk=imfillholesBySlices(~BN);
% Msk=getLargestComponent(Msk);
% Msk=imfillholesBySlices(Rec1pRMS>0.001);
se = strel('disk', 3);
DMsk=imdilate(Msk,se,'same');
DMsk=imfillholesBySlices(DMsk);
try
    SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p),0,TrgSz);
catch
    SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p),8,TrgSz);
end
SelfSens1=SelfSens(:,:,:,1);
% SelfSens1=SelfSens1.*DMsk;
disp('ok SelfSens');

SelfSens1S(:,:,:,SlicesToRead)=SelfSens1;
%% Recon each channel separately, CC
% DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;

disp('ok rec pre channel cc');
%
SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
disp('ok SensCC');

DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
DataCCP=permute(DataPC,[1:3 8 4:7]);
disp('ok 0');

Rec1x=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),perm84(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),SensCC);

% RSensCC=permute(sum(RefSensx(:,:,:,Ord(SlicesToRead)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
% RSensCC=permute43(RSensCC);
% Rec1xR=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),perm84(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),RSensCC);
%% Try to get all TSC, multi-shot, for B0,T2*
CurReps=1:(nRepsHdr-1);

% nTS_THLR=15;
nTS_THLR=nInnerShots*2+1;

TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

nTSMed=nTS_THLR;

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC=FillOnesTo16([TrgSz 1 1 1 1 nTS_THLR]);
%
TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

dTS_planned_ms=2.5;

% nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);
TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
%
nTraj=numel(TrajPartMed);
TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

nPointsNoNav=nTraj;
% floor(50000/AcqDwellTime_us);
NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
NoNavB=zeros(1,nTraj);
NoNavB(1:nPointsNoNav)=1;

% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nTraj,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
% TSBMed(nPointsMed,1)=0;
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
if(~exist('KernsPMMedBase','var'))
    KernsPMMedBase=getKernsFromTrajM(TrajM(1:36,TrajPartMed),Sz,TSBMed);
end
%
KernsPMMed=KernsPMMedBase(:,:,mod(CurReps-1,36)+1,:,:,:,:,:,:,:,:);
% KernsPMMed=getKernsFromTrajM(TrajM(1:nRepsHdr,TrajPartMed),Sz,TSBMed);

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC=FillOnesTo16([Sz 1 1 1 1 nTSMed]);

disp('Prepared TSB,Kerns for Med');
%%
% TSB_THLR=GetTSCoeffsByLinear(nPointsMed,nTS_THLR);
[TSB_THLR, dT_THLR, TimePointsR_THLR]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTS_THLR);
dT_THLR_ms=dT_THLR*NoNavTime_ms;
FirstT_THLR_ms=TimePointsR_THLR(1)*NoNavTime_ms;
TSB_THLR(nPointsMed,1)=0;
TSB_THLRP=permute(TSB_THLR,[3 1 4 5 6 7 2]);
Sz16_THLR=FillOnesTo16([Sz 1 1 1 1 nTS_THLR]);
if(~exist('KernsP_TSTHLRBase','var'))
    KernsP_TSTHLRBase=getKernsFromTrajM(TrajM(1:36,TrajPartMed),Sz,TSB_THLR);
end
KernsP_TSTHLR=KernsP_TSTHLRBase(:,:,mod(CurReps-1,36)+1,:,:,:,:,:,:,:,:);
disp('Prepared TSB,Kerns for THLR');
%%
clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
%     disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
end
%%
THLR_lambda=10;
RhoStr=[' -u ' num2str(1e-3) ' '];
%%
clear THLRMultiShot_RS
% RepsSets={1:(nRepsHdr-1), 1:3};
RepsSets={6:45, [1 13 25]+7};
if(HighRes)
    RepsSets={6:45, [36    21    15    11    40    14    19    25    37]};
end

nRepsSets=numel(RepsSets);
for rs=1 % :nRepsSets
    CurReps=RepsSets{rs};
    disp(['RepSet ' num2str(rs) ' ! ' datestr(now)]);
    tmp=bart(['picsS -w 1 -d 5 -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16_THLR,NoNavB.*DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSB_THLRP,1,...
        sum(KernsP_TSTHLR(:,:,CurReps,:,:,:,:),3));
    THLRMultiShot_RS(:,:,:,rs)=squeeze(tmp);
end
%%
WhichTSToUs=2:12;
WhichTSToUs=2:(nTS_THLR-1);

clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
for rs=1 %:nRepsSets
    [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=...
        FitToModel_MPBD1CSf(THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end
%%
% RepsSets{2}=randi(36,1,9)+5;
%%
T2svalues_ms=linspace(5,300,200);
Decays=exp(-TimePoints_Med_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

nComponentsToUse=3;

clear TSCMed_RS TSCxPMed_RS CompsP_RS
rsRef=1;
for rs=2:nRepsSets
    TSCMed_RS(:,:,:,rs)=exp(-TimePoints_Med_ms3./UpdatedT2SMap_ms_RS(:,:,rsRef)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,rsRef).*TimePoints_Med_ms3/1e3);
    TSCxPMed_RS(:,:,rs,1,1,1,:)=permute(TSCMed_RS(:,:,:,rs),[1 2 7 6 5 4 3]);
    
%     TempRMS=grmss(THLRMultiShot_RS(:,:,:,rs),3:30);
%     CurT2SMap=UpdatedT2SMap_ms_RS(:,:,rs);
%     [HT2s,edges,T2sbin] = histcounts(CurT2SMap(:),T2Sbins);
%     WHT2s=HT2s*0;
%     for i=1:numel(HT2s)
%         WHT2s(i)=sum(TempRMS(T2sbin==i));
%     end
%     SWHT2s=max(WHT2s,sum(WHT2s)*0.03/numel(WHT2s));
% 
%     clear WDecays WDecays4
%     for i=2:numel(HT2s)-1
%         WDecays(i,:)=SWHT2s(i)*exp(-TimePointsMed_ms./T2sCenters(i-1));
%     end
%     [WUd,WSd,WVd]=svd(WDecays,'econ');
%     WVd_RS(:,:,rs)=WVd;
%     
%     CompsP_RS(rs,1,1,1,1,:,:)=permute(WVd_RS(:,1:nComponentsToUse,rs),[7:-1:3 2 1]);
    CompsP_RS(rs,1,1,1,1,:,:)=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
end
disp('Prepared for LLR RepSets');

LLR_lambda=0.1;
RhoStr=[' -u ' num2str(1e-3) ' '];
% RhoStr=[' -u ' num2str(1e-1) ' '];
BlkSz=4;
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);

TSCxPMedOnlyB0_RS=exp(1i.*angle(TSCxPMed_RS));
%%
% [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=FitToModel_MPBD1CSf(THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,TotalAcqTime_ms,TE0_ms);
% save(['For_NU_MPBD2_S' num2str(SliI) '_.mat'],'THLRMultiShot_RS','PDBase_RS','UpdatedB0Map_RS','UpdatedT2SMap_ms_RS','Fitted0_RS');
%%
clear Rec_CompgB0_C_RS
for rs=2:nRepsSets
    CurRep=RepsSets{rs};
%     Rec_CompgB0_C_RS{rs}=bart(['picsS -m -w 1 ' RhoStr ' -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
%     Rec_CompgB0_C_RS{rs}=bart(['picsS -m -S -d 5 -R W:3:0:0.1 ' ScriptFN_CompgBo],...
    Rec_CompgB0_C_RS{rs}=bart(['picsS -m -S -d 5 -R W:3:0:0.001 ' ScriptFN_CompgBo],...
        Sz16CompgB0,NoNavB.*DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed,TSCxPMedOnlyB0_RS(:,:,rs,1,1,1,:),...
        sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP_RS(rs,1,1,1,1,:,:));
end
%%
clear PDBase_RS_LLR UpdatedB0Map_RS_LLR UpdatedT2SMap_ms_RS_LLR s_vals_RS_LLR Fitted0_RS_LLR PDBase0_RS_LLR
Rec_CompgB0_RS_M=cat(4,Rec_CompgB0_C_RS{:});
Rec_CompgB0_RS_MX=perm43(squeeze(sum(Rec_CompgB0_RS_M.*gpermute(CompsP_RS,[4 1]),6)));
Rec_CompgB0_RS_MX=Rec_CompgB0_RS_MX.*perm43(squeeze(TSCxPMedOnlyB0_RS));
Rec_CompgB0_RS_MX=Rec_CompgB0_RS_MX(:,:,:,2);
WhichTSToUse_LLR=2:(nTSMed-1);

for rs=1 % :nRepsSets
    [PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs), Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
        FitToModel_MPBD1CSf(Rec_CompgB0_RS_MX(:,:,:,:,rs),WhichTSToUse_LLR,dT_Med_ms,TE0_ms+FirstT_Med_ms(1));
end
% ShowAbsAngle(Rec_CompgB0_RS_MX(:,:,2:5:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('LLR');
% ShowAbsAngle(Fitted0_RS_LLR(:,:,2:5:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('1CS Fitted LLR');
% ShowAbsAngle(Rec_CompgB0_RS_MX(:,:,2:5:end,2));ylabel('3 shots','FontSize',16);title('LLR');
% ShowAbsAngle(Fitted0_RS_LLR(:,:,2:5:end,2));ylabel('3 shots','FontSize',16);title('1CS Fitted LLR');
%%
CurReps=RepsSets{2};
% CurSig=DataCCP(:,TrajPartMed,CurReps,1:nccToUse);
% CurKernsA=KernsPMMed(:,:,CurReps,:,:,:,:);
CurSens=SensCC(:,:,:,1:nccToUse);
% CurTraj=STraj3MMed(:,:,CurReps);

SliI=SlicesToRead;

sccmtxS(:,:,SliI)=sccmtx;

save([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat'],'CurSens','Fitted0_RS_LLR',...
    'UpdatedB0Map_RS_LLR','UpdatedT2SMap_ms_RS_LLR','TotalAcqTime_ms','TE0_ms','THLRMultiShot_RS','Fitted0_RS','Rec_CompgB0_RS_MX',...
    'sccmtx');
%%
% DataForSlice{SlicesToRead}={UpdatedB0Map,UpdatedT2SMap_ms,THLRMultiShot,SelfSens1,sccmtx};

close all;

end % Slices loop
%%
SmallVarT=50e6;
getSmallVars;
save([mainP filesep 'SmallVars.mat'],SmallVars{:});

save([mainP filesep 'KernsPMMedBase' '.mat'],'KernsPMMedBase');
save([mainP filesep 'KernsP_TSTHLRBase' '.mat'],'KernsP_TSTHLRBase');
save([mainP filesep 'TrajM' '.mat'],'TrajM');
%%
% save([mainP filesep 'DataForSlice.mat'],'DataForSlice');
% % load([mainP filesep 'DataForSlice.mat']);
% %%
% clear KernsByRepMed              	1513758480    		1513.76MB
% clear KernsPMMed                 	1513754112    		1513.75MB
% clear ADataIsL                   	499200000     		499.20MB
% clear DataForSlice               	492372000     		492.37MB
% clear ADataIsLCC                 	483600000     		483.60MB
% clear DataC                      	483600000     		483.60MB
% clear DataCCP                    	439640760     		439.64MB
% clear DataPC                     	439640760     		439.64MB
% clear InnerShotData              	193440000     		193.44MB
% clear InnerShotDataIn            	193440000     		193.44MB
% save([mainP filesep 'CurStatus.mat']);
%%
clear Fitted0_RS_LLRS Fitted0_RSS UpdatedB0Map_RS_LLRS UpdatedT2SMap_ms_RS_LLRS PDBase0_RS_LLRS UpdatedT2SMap_ms_RS_LLRSb
clear PDBase_RS_LLR UpdatedB0Map_RS_LLR UpdatedT2SMap_ms_RS_LLR s_vals_RS_LLR Fitted0_RS_LLR PDBase0_RS_LLR
clear UpdatedB0Map_RSS UpdatedT2SMap_ms_RSS THLRMultiShot_RSS
SliI=7;
for SliI=1:nSlices
    disp(SliI);
    QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
    WhichTSToUse_LLR=2:(nTSMed-1);
    Fitted0_RS_LLRS(:,:,:,SliI)=QQ.Fitted0_RS_LLR(:,:,:);
    Fitted0_RSS(:,:,:,SliI)=QQ.Fitted0_RS;
    UpdatedB0Map_RS_LLRS(:,:,:,SliI)=QQ.UpdatedB0Map_RS_LLR;
    UpdatedT2SMap_ms_RS_LLRS(:,:,SliI)=QQ.UpdatedT2SMap_ms_RS_LLR;
% end
% %%
    rs=1;
    [PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs), Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
            FitToModel_MPBD1CSf(QQ.Rec_CompgB0_RS_MX(:,:,:,rs),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    PDBase0_RS_LLRS(:,:,SliI)=PDBase0_RS_LLR(:,:,rs);
    UpdatedT2SMap_ms_RS_LLRSb(:,:,SliI)=UpdatedT2SMap_ms_RS_LLR(:,:,rs);
    
%     WhichTSToUs=2:12;
    WhichTSToUs=2:(nTS_THLR-1);

    clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
    for rs=1 % :nRepsSets
        [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RSS(:,:,SliI,rs)]=...
            FitToModel_MPBD1CSf(QQ.THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    end
    UpdatedB0Map_RSS(:,:,:,SliI)=UpdatedB0Map_RS;
    UpdatedT2SMap_ms_RSS(:,:,:,SliI)=UpdatedT2SMap_ms_RS;
% end
% %%
% for SliI=1:nSlices
%     disp(SliI);
    QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
    THLRMultiShot_RSS(:,:,:,SliI)=QQ.THLRMultiShot_RS;
end
%%
MidEchoI=(nTSMed+1)/2;

MskS=(grmss(THLRMultiShot_RSS(:,:,:,:),3)>1e-4);
fgmontagex(Fitted0_RS_LLRS(:,:,12,ROrd))
fgmontagex(Fitted0_UpdatedT2SMap_ms_RS_LLRSRS_LLRS(:,:,12,:))
fgmontagex(Fitted0_RSS(:,:,9,:))
fgmontagex(Fitted0_RS_LLRS(:,:,MidEchoI,:));xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);
fgmontagex(Fitted0_RS_LLRS(:,:,MidEchoI,ROrd),'Size',[4 nSlices/4]);title('Fitted0_RS_LLRS');xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(gflip(Fitted0_RS_LLRS(:,:,MidEchoI,ROrd),1),'Size',[4 nSlices/4]);title('Fitted0_RS_LLRS');xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(Fitted0_RS_LLRS(:,:,MidEchoI+[-2 2],ROrd([6 11 16])),[0 5e-4]);xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(UpdatedT2SMap_ms_RS_LLRSb(:,:,12,:),[0 200]);colormap hot

fgmontagex(gflip(Fitted0_RS_LLRS(:,:,MidEchoI+[-2 2],ROrd([6 11 16])),1),[0 9e-4]);xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

Fitted0_RS_LLRSm=Fitted0_RS_LLRS.*perm43(MskS);

% ReFitted=perm43(PDBase0_RSS.*MskS).*exp(-perm32(0:1:49)./UpdatedT2SMap_ms_RSS);
ReFitted=perm43(PDBase0_RS_LLRS.*MskS.*exp(-perm42(0:1:49)./UpdatedT2SMap_ms_RS_LLRS));
ReFitted=ReFitted(:,:,:,ROrd);
ReFitted=ReFittedb;
%%
TE1=20;
TE2=50;
ForFig=cat(4,ReFitted(:,:,TE1,(32:-1:17)),perm43(ReFitted(:,:,TE1:2:TE2,(17))),ReFitted(:,:,TE2,(16:-1:5)),perm43(ReFitted(:,:,TE2:-2:TE1,(5))),ReFitted(:,:,TE1,(5:32)));
Mx=6e-4;
QQ=min(max(abs(ForFig),0),Mx)/Mx*256;
filename = 'testAnimated2.gif'; % Specify the output file name
imwrite(QQ,gray(256),filename,'gif','LoopCount',Inf,'DelayTime',.07);
%%
UpdatedT2SMap_ms_RS_LLRSbm=UpdatedT2SMap_ms_RS_LLRSb.*MskS;
fgmontagex(UpdatedT2SMap_ms_RSS,[0 200]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);
fgmontagex(UpdatedT2SMap_ms_RS_LLRSbm(:,:,ROrd),[0 200],'Size',[4 8]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(UpdatedT2SMap_ms_RS_LLRSbm(:,:,ROrd([7 17 25])),[0 200],'Size',[1 3]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(UpdatedT2SMap_ms_RS_LLRSbm(:,:,1:16),[0 200],'Size',[1 3]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(gflip(UpdatedT2SMap_ms_RS_LLRSbm(:,:,1:16),1),[0 200],'Size',[1 3]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

fgmontagex(gflip(UpdatedT2SMap_ms_RS_LLRSbm(:,:,ROrd),1),[0 200],'Size',[4 6]);colormap hot;xlabel(mainP,'Interpreter','None');set(get(gcf,'Children'),'Position',[0.05 0.05 0.9 0.9]);

%%





%%
for SliI=1:nSlices
    disp(SliI);
    Rec1S(SliI)=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
end
%%
for SliI=1:nSlices
    disp(SliI);
    THLRMultiShot_RSS(:,:,:,SliI)=Rec1S(SliI).THLRMultiShot_RS;
    Fitted0_RS_LLRS(:,:,:,SliI)=Rec1S(SliI).Fitted0_RS_LLR;
    UpdatedT2SMap_ms_RS_LLRS(:,:,SliI)=Rec1S(SliI).UpdatedT2SMap_ms_RS_LLR;
end
%%
for SliI=1:nSlices
    SensCCS(:,:,:,:,SliI)=Rec1S(SliI).CurSens;
end
%% get Data and sccmtxS
RepsToRead=1:(nRepsHdr);
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    disp(SliI);
    
    ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepsToRead,:,:,:,:,:,:,:,:);
    ADataIsL=squeeze(ADataIsL);
    ADataIsL=CombineDims(ADataIsL,[4 1]);

%     ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,:,:,:,:,:,:,:);
%     ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
%     ADataIsL=CombineDims(ADataIsL,[4 1]);
    disp([datestr(now) ' Read data']);
    
    ChRMS=grmss(ADataIsL,[1 3]);
    [SChRMS ChOrd]=sort(ChRMS,'descend');
    
    Ch2D=CombineDims(ADataIsL,[3 1]);
    [~,S,sccmtx] = svd(Ch2D(1:10:end,:),'econ');
    clear Ch2D
    
    sccmtxS(:,:,SliI)=sccmtx;
    
    ncc=31;
    
    ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
    for i=1:ncc
        ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
    end
    DataC=permute(ADataIsLCC,[1 3 2]);
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    
    DataCCPS(SliI,:,:,:,:)=DataCCP;
    if(SliI==1)
        DataCCPS(nSlices,1,1,1,1)=0;
    end
end
%%
for i=1:ncc
    disp(['Compressed channel ' num2str(i)]);
    CurChData=DataCCPS(:,1,:,:,i);
    save([mainP filesep 'DataCCPS_Ch' num2str(i) '.mat'],'CurChData');
end
%%
for SliI=1:nSlices
    sccmtx=sccmtxS(:,:,SliI);
    RSensCC=permute(sum(RefSensx(:,:,:,Ord(SliI)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    RSensCC=permute43(RSensCC);
    RSensCCS(:,:,:,SliI)=RSensCC;
end
SensCCS=permute(RSensCCS,[1 2 5 3 4]);
%%
% CurReps=1:(nRepsHdr-1);
% 
% % CurReps=1:39;
% 
% % nTS_THLR=15;
% 
% TrajPartMed=1:nTrajToUse;
% 
% nPointsMed=numel(TrajPartMed);
% 
% nTSMed=nTS_THLR;
% 
% ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
% Sz16AllTSC=FillOnesTo16([TrgSz 1 1 1 1 nTS_THLR]);
% %
% TrajPartMed=1:nTrajToUse;
% 
% nPointsMed=numel(TrajPartMed);
% 
% dTS_planned_ms=2.5;
% % dTS_planned_ms=0.5;
% dTS_planned_ms=8;
% 
% % nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);
% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
% TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
% 
% TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
% TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
% %
% nTraj=numel(TrajPartMed);
% TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;
% 
% % nPointsNoNav=floor(50000/AcqDwellTime_us);
% % NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
% % NoNavB=zeros(1,nTraj);
% % NoNavB(1:nPointsNoNav)=1;
% 
% nPointsNoNav=nTraj;
% NoNavTime_ms=TotalAcqTime_ms;
% NoNavB=ones(1,nTraj);
% 
% % TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
% % [TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTSMed);
% [TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTSMed);
% dT_Med_ms=dT_Med*NoNavTime_ms;
% FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
% TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
% TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
% TSBMed(nPointsMed,1)=0;
% TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
% 
% RepsSets={1:(nRepsHdr-1), [1 13 25]};
% RepsSets={1:(nRepsHdr-1), [6 18 30]};
% 
% RepsToCompute=RepsSets{2};
% if(~exist('KernsPMMedx','var'))
% %     KernsPMMed=getKernsFromTrajM(TrajM(1:(nRepsHdr-1),TrajPartMed),Sz,TSBMed);
% %     KernsPMMedx=getKernsFromTrajM(TrajM(RepsToCompute,TrajPartMed),Sz,TSBMed);
% end
% KernsPMMedy=getKernsFromTrajM(TrajM(1:(nRepsHdr-1),TrajPartMed),Sz,TSBMed);
% 
% ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
% Sz16AllTSC=FillOnesTo16([Sz 1 1 1 1 nTSMed]);
% 
% disp('Prepared TSB,Kerns for Med');
%%
KernsPMMedy=KernsPMMedBase(:,:,mod((1:(nRepsHdr-1))-1,36)+1,:,:,:,:,:,:,:,:);
%%
WhichTSToUs=3:nTS_THLR-3;
% WhichTSToUs=3:7
clear PDBase_RSS UpdatedB0Map_RSS UpdatedT2SMap_ms_RSS s_vals_RSS Fitted0_RSS PDBase0_RSS
for SliI=1:nSlices % :nRepsSets
    disp(SliI);
    [PDBase_RSS(:,:,SliI), UpdatedB0Map_RSS(:,:,SliI), UpdatedT2SMap_ms_RSS(:,:,SliI), s_vals_RSS(:,:,:,SliI), Fitted0_RSS(:,:,:,SliI), PDBase0_RSS(:,:,SliI)]=...
        FitToModel_MPBD1CSf(THLRMultiShot_RSS(:,:,:,SliI),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end

for SliI=1:nSlices % :nRepsSets
    disp(SliI);
    UpdatedB0Map_RSxS(:,:,SliI)=SmoothByW(UpdatedB0Map_RSS(:,:,SliI),grmss(s_vals_RSS(:,:,:,SliI),3));
end
disp('Prepared sB0');
%%
nComponentsToUse=3;

T2svalues_ms=linspace(5,300,200);
Decays=exp(-TimePoints_Med_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

clear TSCMed TSCxPMed CompsP

CompsP(1,1,1,1,1,:,:)=permute(Vd(:,1:nComponentsToUse,1),[7:-1:3 2 1]);
%%
SliIs=ROrd(8:11);
SliIs=ROrd;
nSliIs=numel(SliIs);
TSCMedS=exp(-1i*2*pi*perm93(UpdatedB0Map_RSxS(:,:,SliIs)).*TimePoints_Med_ms3/1e3);
TSCxPMedS=perm73(TSCMedS);

disp('Prepared for LLR RepSets');

LLR_lambda=0.1;
% RhoStr=[' -u ' num2str(1e-3) ' '];
RhoStr=[' -u ' num2str(1e-1) ' '];
BlkSz=4;
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
Sz16CompgB0S=FillOnesTo16([Sz 1 1 1 nComponentsToUse 1 1 numel(SliIs)]);

TSCxPMedOnlyB0S=exp(1i.*angle(TSCxPMedS));
%%
nccToUse=15;
RepSetX=[1 13 25; 6 18 30; 4 16 28; 9 21 33];
RepSetX=repmat(RepSetX,[8 1]);

% 3-shot
BaseRepSetX=[1 13 25];
BIdxs= [randi(3,nSliIs,1), randi(2,nSliIs,1)*2-3];
BIdxs=mod([BIdxs(:,1) sum(BIdxs,2)],3)+1;
BIdxs(:,3)=6-sum(BIdxs,2);
RepSetX1=BaseRepSetX(BIdxs);
% RepSetX=mod([1 13 25]-5+(1:numel(SliIs)).'*5,36)+6;
RepSetX=mod(RepSetX1+((1:nSliIs)-1).'*7,36)+8;

% % 2-shot
% BaseRepSetX=[1 10];
% BIdxs=randi(2,nSliIs,1);
% BIdxs=[BIdxs 3-BIdxs];
% RepSetX1=BaseRepSetX(BIdxs);
% RepSetX=mod(RepSetX1+((1:nSliIs)-1).'*1,36)+8;
% 
% % 1-shot
% RepSetX=(8+(1:nSliIs)).';

% 9-shot
RepSetX=randi(36,nSliIs,9)+8;

CurSens=perm95(SensCCS(:,:,:,1:nccToUse,SliIs));
clear CurData CurKerns CurSTraj
disp('Preparing data...');
for i=1:nSliIs
%     disp(i);
    CurSTraj(:,:,:,:,:,:,:,:,i)=STraj3MMed(:,:,RepSetX(i,:));
    CurKerns(:,:,:,:,:,:,:,1,i)=sum(KernsPMMedy(:,:, RepSetX(i,:) ,:,:,:,:),3);
    CurData(1,:,:,:,1,1,1,1,i)=DataCCPS(SliIs(i),:,:,RepSetX(i,:),1:nccToUse);
end
if(size(RepSetX,2)==1)
    CurData=perm43(CurData);
end
disp('Prepared data');
%% High res: 3-shot 1735 sec 9shot 1988.327092
%     RecStr=[' -m -S -d 5 -R T:259:0:0.1 '];
RecStr=[' -m -S -d 5 -R W:3:256:0.1 '];
RecStr=[' -m -S -d 5 -R W:259:0:0.0001 '];
RecStr=[' -m -S -d 5 -R T:256:0:0.1 '];
Rec_CompgB0_C_RSS=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData,...
    CurSens,CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);
%%
Rec_CompgB0_RSS_M=cat(4,Rec_CompgB0_C_RSS);
Rec_CompgB0_RSS_MX=perm43(squeeze(sum(Rec_CompgB0_RSS_M.*gpermute(CompsP,[4 1]),6)));
Rec_CompgB0_RSS_MX=Rec_CompgB0_RSS_MX.*perm43(squeeze(TSCxPMedOnlyB0S));
%% That's it
Rec_CompgB0_RSS_MX_WithTVOverZ=Rec_CompgB0_RSS_MX;
save([mainP filesep 'Rec_CompgB0_RSS_MX_WithTVOverZ.mat'],'Rec_CompgB0_RSS_MX_WithTVOverZ');
%%
[PDBase_RS_A, UpdatedB0Map_RS_A, UpdatedT2SMap_ms_RS_A, s_vals_RS_A, Fitted0_RS_A, PDBase0_RS_A]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);

[PDBase_RS_LLRSX, UpdatedB0Map_RS_LLRSX, UpdatedT2SMap_ms_RS_LLRSX, s_vals_RS_LLRSX, Fitted0_RS_LLRSX, PDBase0_RS_LLRSX]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX_WithTVOverZ),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);

[PDBase_RS_LLRSX0, UpdatedB0Map_RS_LLRSX0, UpdatedT2SMap_ms_RS_LLRSX0, s_vals_RS_LLRSX0, Fitted0_RS_LLRSX0, PDBase0_RS_LLRSX0]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX0),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    
[PDBase_RS_LLRSX3, UpdatedB0Map_RS_LLRSX3, UpdatedT2SMap_ms_RS_LLRSX3, s_vals_RS_LLRSX3, Fitted0_RS_LLRSX3, PDBase0_RS_LLRSX3]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX3),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    
[PDBase_RS_LLRSX2, UpdatedB0Map_RS_LLRSX2, UpdatedT2SMap_ms_RS_LLRSX2, s_vals_RS_LLRSX2, Fitted0_RS_LLRSX2, PDBase0_RS_LLRSX2]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX2),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    
[PDBase_RS_LLRSX1, UpdatedB0Map_RS_LLRSX1, UpdatedT2SMap_ms_RS_LLRSX1, s_vals_RS_LLRSX1, Fitted0_RS_LLRSX1, PDBase0_RS_LLRSX1]=...
        FitToModel_MPBD1CSf(perm43(Rec_CompgB0_RSS_MX1),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
%%
ReFittedb=perm43((PDBase0_RS_A.*MskS(:,:,SliIs)).*exp(-perm42(0:1:49)./UpdatedT2SMap_ms_RS_A));

%%
BynShots=cat(5,Fitted0_RS_LLRSX1,Fitted0_RS_LLRSX2,Fitted0_RS_LLRSX3);

fgmontagex(perm43(squeeze(BynShots(10:end-10,10:end-10,8,[6 11 16 25],:))))
%%
SmallVarT=100e6;
getSmallVars;
save([mainP filesep 'SmallVars2.mat'],SmallVars{:});

    %%
MskS=(grmss(s_vals_RS_A(:,:,:,:),3)>1e-4);
UpdatedT2SMap_ms_RS_A=UpdatedT2SMap_ms_RS_A.*MskS;
fgmontagex(UpdatedT2SMap_ms_RS_A(:,:,:),[0 200]);colormap hot
fgmontagex(gflip(perm31(UpdatedT2SMap_ms_RS_A([60],:,:)),1),[0 200]);colormap hot;daspect([1 2 1])
%%
QQ=min(max(UpdatedT2SMap_ms_RS_LLRSX0,0),200)/200*256;
filename = 'testAnimated.gif'; % Specify the output file name
imwrite(perm43(QQ),hot(256),filename,'gif','LoopCount',Inf,'DelayTime',.25);
%%
for SliI=1:nSlices
    disp(SliI);
    [PDBase_RS_LLRSX(:,:,SliI), UpdatedB0Map_RS_LLRSX(:,:,SliI), UpdatedT2SMap_ms_RS_LLRSX(:,:,SliI), s_vals_RS_LLRSX(:,:,:,SliI), Fitted0_RS_LLRSX(:,:,:,SliI), PDBase0_RS_LLRSX(:,:,SliI)]=...
        FitToModel_MPBD1CSf(squeeze(Rec_CompgB0_RSS_MX_WithTVOverZ(:,:,SliI,:)),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end
%%
fgmontagex(UpdatedT2SMap_ms_RS_LLRSX,[0 200]);colormap hot
fgmontagex(UpdatedT2SMap_ms_RS_LLRSX(:,:,1:6),[0 200]);colormap hot
fgmontagex(UpdatedT2SMap_ms_RS_LLRSX(:,:,21:26),[0 200]);colormap hot

fgmontagex(rot90(perm32(UpdatedT2SMap_ms_RS_LLRSX(:,[40 70],:))),[0 200]);colormap hot;daspect([1 2 1]);
fgmontagex(gflip(perm31(UpdatedT2SMap_ms_RS_LLRSX([40 60],:,:)),1),[0 200]);colormap hot;daspect([1 2 1]);
%% LLR
RecStr=[' -m -S -d 5 -b 4 -R L:259:259:10.1 '];
Rec_CompgB0_C_RSSL=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData,...
    CurSens,CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

Rec_CompgB0_RSS_ML=cat(4,Rec_CompgB0_C_RSSL);
Rec_CompgB0_RSS_MXL=perm43(squeeze(sum(Rec_CompgB0_RSS_ML.*gpermute(CompsP,[4 1]),6)));
Rec_CompgB0_RSS_MXL=Rec_CompgB0_RSS_MXL.*perm43(squeeze(TSCxPMedOnlyB0S));

%% Multi lambda,rho test: doesn't change much. 1e-1 and lower good for both
Rhos=10.^(-3:1);
Lambdas=10.^(-3:1);

for i=1:numel(Rhos)
    for j=1:numel(Lambdas)
        disp([i j]);
        disp(datestr(now));
        RecStrS=[' -m -S -u ' num2str(Rhos(i)) ' -R W:259:0:' num2str(Lambdas(j)) ' '];
        Rec_CompgB0_C_RSSC{i,j}=bart(['picsS ' RecStrS ScriptFN_CompgBo],...
            Sz16CompgB0S,CurData,...
            CurSens,CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
            CurKerns,CompsP);
    end
end
%%
clear Rec_CompgB0_C_RSSCM
for i=1:numel(Rhos)
    Rec_CompgB0_C_RSSCM{i}=cat(4,Rec_CompgB0_C_RSSC{i,:});
end
Rec_CompgB0_C_RSSCM=cat(5,Rec_CompgB0_C_RSSCM{:});
Rec_CompgB0_C_RSSCMX=sum(Rec_CompgB0_C_RSSCM.*CompsP,6);
%%
    WhichTSToUse_LLR=2:20;
    Fitted0_RS_LLRS(:,:,:,SliI)=QQ.Fitted0_RS_LLR(:,:,:,2);
    Fitted0_RSS(:,:,:,:,SliI)=QQ.Fitted0_RS;
    UpdatedB0Map_RS_LLRS(:,:,:,SliI)=QQ.UpdatedB0Map_RS_LLR;
    UpdatedT2SMap_ms_RS_LLRS(:,:,:,SliI)=QQ.UpdatedT2SMap_ms_RS_LLR;
    rs=2;
    [PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs), Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
            FitToModel_MPBD1CSf(QQ.Rec_CompgB0_RS_MX(:,:,:,rs),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    PDBase0_RS_LLRS(:,:,SliI)=PDBase0_RS_LLR(:,:,rs);
    
    WhichTSToUs=2:12;
    clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
    for rs=1:nRepsSets
        [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=...
            FitToModel_MPBD1CSf(QQ.THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    end
    UpdatedB0Map_RSS(:,:,:,SliI)=UpdatedB0Map_RS;
    UpdatedT2SMap_ms_RSS(:,:,:,SliI)=UpdatedT2SMap_ms_RS;
% end