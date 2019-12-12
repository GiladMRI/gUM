ScanP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68xms/';
% BaseFN='meas_MID00663_FID30522_gSpi2d_T10_Dw11_d110_VD1';
RefFldMapP=[ScanP 'meas_MID00654_FID30513_gre_4echo_32_22' filesep];
BaseFN='meas_MID00722_FID30571_gSpi2d_T12_Dw11_d110_VD1_Cor';

ScanP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
RefFldMapP=[ScanP 'meas_MID00856_FID32091_gre_4echo_32_22' filesep];
BaseFN='meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1';
BaseFN='meas_MID00864_FID32099_gSpi2d_T10_Dw11_d120_VD1';
% BaseFN='meas_MID00884_FID32119_gSpi2d_T12_Dw11_d110_VD1_Cor';
% BaseFN='meas_MID00918_FID32147_gSpi2d_T10_Dw11_d110_VD1_Sag';

ScanP='/autofs/cluster/kawin/Gilad/VD1PAT3_CL_CylPhantom/';
BaseFN='meas_MID00605_FID37709_gSpi2d_T10_Dw11_d110_VD1';
% BaseFN='meas_MID00608_FID37712_gSpi2d_T14_Dw11_d110_VD102';
% BaseFN='meas_MID00610_FID37714_gSpi2d_T15_Dw11_d110_VD102';

ScanP='/autofs/cluster/kawin/Gilad/gep2d_Try1_Phantom_Traj1/';
BaseFN='meas_MID01323_FID41637_gSpi2d_T17_Dw11_d110_VD102';

MIDStr=BaseFN(6:13);
FN=[ScanP BaseFN '.dat'];
disp('ok');
%%
mainP=[ScanP BaseFN];
mkdir(mainP);

system(['chmod +777 -R ' mainP]);
disp([mainP ' Created']);

setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/bart-0.4.04b')

BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';
ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];

nccToUse=31;

Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];

[~, HostName]=system('hostname');
HostName=HostName(1:(find(HostName=='.',1)-1));

disp('ok 2');

ADataxFN=[mainP filesep 'ADatax.mat'];
%% Read raw
if(exist(ADataxFN,'file'))
    load(ADataxFN);
else
    AData = mapVBVD(FN);
    if(iscell(AData))
        ADatax=AData{end};
    else
        ADatax=AData;
    end
    save(ADataxFN,'ADatax');
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
TE0_ms=2.38;

TrajType=WipMemBlock.adFree{12};
% ResType=floor((TrajType-10)/2)+1; % 1.9,1.3
ResType=1;
TimingType=mod(TrajType-10,8)+1; % 6ms, 8ms
% G1=load('GAll68x.mat');
% G2=load('GAll68cor.mat');
% GAll=cat(2,G1.GAll, G2.GAll);

% load('GAll1p9mmVD1PAT3.mat');
load('GAll1p9mmVD1PAT3Pause.mat');

GNav=load('GNav1ms.mat');
GNav=GNav.GNav;
GTrajaCBase=GAll(:,TimingType,ResType);

if(TimingType==1 || TimingType==3)
    nInnerShots=8;
else
    nInnerShots=6;
end
disp('Read traj base');
%
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

FOVx=220;

dFOV=FOVx/1000;

FOV_mm=FOVx;

disp('ok');

GradReduceFac=WipMemBlock.adFree{4};
PhiRotPerRep=WipMemBlock.adFree{8};
AcqDwellTime_us=WipMemBlock.adFree{13}/1000;

dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;
%%
RLocs=load([RefFldMapP 'Locs.mat']);
RefLocs=RLocs.RotatedLocs;
CurzLocs=RotatedLocs(3,:);
RefzLocs=RefLocs(3,:);

RefMaps=load([RefFldMapP 'B0T2S.mat']);
RefSens=load([RefFldMapP 'Sens.mat']);
disp('Loaded ref');
%%
GTrajaC=GTrajaCBase/GradReduceFac;
GNavC=GNav/GradReduceFac;

for j=1:5
    GTrajaCM(:,j)=[zeros((j-1)*50,1); GTrajaC; zeros((5-j)*50,1); GNavC; -GNavC];
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
res_mm=FOVx/(max(abs(kK))*2);
%
figure;
plot(kK);
axis square;
axis equal;
ylabel('k');
%%
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
    AcqTimePoints_us=Extra*GradDwellTime_us+Delay_us+(0:AcqDwellTime_us:52000);

    nRepsToUse=nRepsHdr-1;
    
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
SlicesToRead=[5];
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
RepsToRead=1:(nRepsHdr-1);
% ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,ADCsToRead,:,:,:,:,:,:);
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,:,:,:,:,:,:,:);

ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);

disp(['Read data ' num2str(SliI)]);
%%
ChRMS=grmss(ADataIsL,[1 3]);
[SChRMS ChOrd]=sort(ChRMS,'descend');
% Data=squeeze(ADataIsL(:,ChOrd(1),1:nRepsToUse));

Ch2D=CombineDims(ADataIsL,[3 1]);
[~,S,sccmtx] = svd(Ch2D(1:end,:),'econ');
clear Ch2D 
sccmtxS(:,:,SliI)=sccmtx;


ncc=31;
% ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
% for i=1:ncc
%     ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
% end
% DataC=permute(ADataIsLCC,[1 3 2]);
DataC=perm43(sum(perm32(ADataIsL).*permute(sccmtx(:,1:ncc),[3 4 1 2]),3));
disp('ok cc');

%     DataC=permute(ADataIsLCC,[1 3 2]);
    %% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
TrajPartToUse=0+(1:2000);
RepsToUse=1:(nRepsHdr-1);

DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 3 5 6 7 8 2]).*modx;
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 size(DataPC,8)]);

Rec1p=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
% fgmontage(grmss(Rec1p,8));removeTicks;
% title([num2str(numel(RepsToUse)) ' shots data x ' num2str(numel(TrajPartToUse)) 'points, pre channel recon, RMS']);
disp('ok per channel');
%%
Rec1pRMS=grmss(Rec1p,8);
[Out B BN]=CalcSlicesSNR(Rec1pRMS,false,7);
Msk=imfillholesBySlices(~BN);
% Msk=getLargestComponent(Msk);
% Msk=imfillholesBySlices(Rec1pRMS>0.001);
se = strel('disk', 3);
DMsk=imdilate(Msk,se,'same');
DMsk=imfillholesBySlices(DMsk);
SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p).*DMsk,0,TrgSz);
SelfSens1=SelfSens(:,:,:,1);
disp('ok SelfSens');
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

DataCCP1S(:,:,SliI,:)=DataCCP(:,:,1,:);

SensCCS(:,:,SliI,:)=SensCC;

SelfSens1S(:,:,SliI,:)=SelfSens1;
end
%%
save([mainP filesep 'sccmtxS.mat'],'sccmtxS');
save([mainP filesep 'SensCCS.mat'],'SensCCS');
save([mainP filesep 'SelfSens1S.mat'],'SelfSens1S');
save([mainP filesep 'DataCCP1S.mat'],'DataCCP1S');
%%
load([mainP filesep 'sccmtxS.mat'],'sccmtxS');
load([mainP filesep 'SensCCS.mat'],'SensCCS');
load([mainP filesep 'SelfSens1S.mat'],'SelfSens1S');
load([mainP filesep 'DataCCP1S.mat'],'DataCCP1S');
%%
TrajPartToUseX=1:14536; % 29082 % 43627
TrajPartToUseXC={1:14080 14500+(1:14080) 29100+(1:14080)};
nCCToUseX=1:25;

for SliI=1:nSlices
for i=1:3
%     RecC{i}=bart('pics -S -u -R T:3:0:0.0000000001 -t ',BARTTrajP(:,TrajPartToUseXC{i},1),DataCCP(:,TrajPartToUseXC{i},1,nCCToUseX),SensCC(:,:,:,nCCToUseX));
%     RecC{i}=bart('pics -S -u -R T:3:0:0.0000000001 -t ',BARTTrajP(:,TrajPartToUseXC{i},1),permute(DataCCP1S(TrajPartToUseXC{i},nCCToUseX,SliI),[4 1 3 2]),SensCCS(:,:,SliI,nCCToUseX));
    RecC{i}=bart('pics -S -u -R T:3:0:0.0001 -t ',BARTTrajP(:,TrajPartToUseXC{i},1),permute(DataCCP1S(TrajPartToUseXC{i},nCCToUseX,SliI),[4 1 3 2]),SensCCS(:,:,SliI,nCCToUseX));
end
RecCM=cat(3,RecC{:});
RecCMS(:,:,:,SliI)=RecCM;
end
%%
save([mainP filesep 'RecCMS.mat'],'RecCMS');
disp('saved RecCMS');
%%
load([mainP filesep 'RecCMS.mat']);
%%
RefSpiP='/autofs/cluster/kawin/Gilad/VD1PAT3_CL_CylPhantom/meas_MID00605_FID37709_gSpi2d_T10_Dw11_d110_VD1/';
for SliI=1:nSlices
    disp(SliI);
    QQ=load([RefSpiP 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
    THLRMultiShot_RSS(:,:,:,SliI)=QQ.THLRMultiShot_RS;
end

fgmontagex(THLRMultiShot_RSS(:,:,[2 5 11],ROrd));
%%

%%
Both=cat(5,RecCMS(:,:,:,ROrd),padarray(THLRMultiShot_RSS(:,:,[2 5 11],ROrd),[2 2],'Both'));
Both=Both./grms(Both,1:3);
%% Fit THLR ref
nTS_THLR=15;
nPointsNoNav=floor(50000/AcqDwellTime_us);
NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;

[TSB_THLR, dT_THLR, TimePointsR_THLR]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTS_THLR);
dT_THLR_ms=dT_THLR*NoNavTime_ms;
FirstT_THLR_ms=TimePointsR_THLR(1)*NoNavTime_ms;

WhichTSToUs=2:12;
for SliI=1:nSlices
    disp(SliI);
    [~, UpdatedB0MapTHLRS(:,:,SliI), UpdatedT2SMap_msTHLRS(:,:,SliI), s_valsTHLRS(:,:,:,SliI), ~, PDBase0THLRS(:,:,SliI)]=...
        FitToModel_MPBD1CSf(THLRMultiShot_RSS(:,:,:,SliI),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end
%% get good B0 from ref
SRange=[21 21];
SSigs=[0.0001 0.1:0.1:29];
MaxB0=400;
M0MedianFac=20;
for SliI=1:nSlices
    disp(SliI);
    W=s_valsTHLRS(:,:,1,SliI);
    W=min(W,median(W(:))*M0MedianFac);
    B0=UpdatedB0MapTHLRS(:,:,SliI);
    B0=min(max(B0,-MaxB0),MaxB0);
    
    SB0(:,:,SliI)=SmoothByW(B0,W,SRange,SSigs);
end
disp('got SB0');
SB0r=imresize(SB0,Sz);
%%
save([RefSpiP 'FitAnalysis.mat'],'UpdatedB0MapTHLRS','UpdatedT2SMap_msTHLRS','s_valsTHLRS','PDBase0THLRS','SB0');
disp(['Saved ' RefSpiP 'FitAnalysis.mat']);
%% Preparations
nTS_THLR=15;

TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

dTS_planned_ms=2.5;

nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);

% nTSMed=6;
% nTSMed=81;
% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
% TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
%
nTraj=numel(TrajPartMed);
TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

nPointsNoNav=floor(50000/AcqDwellTime_us);
NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
NoNavB=zeros(1,nTraj);
NoNavB(1:nPointsNoNav)=1;

% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
TSBMed(nPointsMed,1)=0;
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

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
disp('Prepared TSB,Kerns for THLR');
%%
clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
%     disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
    STraj3M(:,:,CurRep)=CTo3Rows(TrajM(CurRep,:));
end
%%
THLR_lambda=10;
RhoStr=[' -u ' num2str(1e-3) ' '];
%%
RepsForKerns=[1 13 25];
% RepsForKerns=1:3;
% RepsForKerns=1:(nRepsHdr-1);
% if(~exist('KernsPMMed','var'))
    KernsPMMed=getKernsFromTrajM(TrajM(RepsForKerns,TrajPartMed),Sz,TSBMed);
%     KernsP_TSTHLR=getKernsFromTrajM(TrajM(RepsForKerns,TrajPartMed),Sz,TSB_THLR);
% end
%% Fit CMS
WhichTSToUseA=1:3;
dTA=14500*AcqDwellTime_us/1000;
for SliI=1:nSlices
    disp(SliI);
    [~, UpdatedB0MapA(:,:,SliI), UpdatedT2SMap_msA(:,:,SliI), s_valsA(:,:,:,SliI), FittedA(:,:,:,SliI), PDBase0A(:,:,SliI)]=...
        FitToModel_MPBD1CSf(RecCMS(:,:,:,SliI),WhichTSToUseA,dTA,TE0_ms);
end
%% Now LLR

%% Recon parts of the traj
CurRep=1;
CurRepSet=RepsForKerns(CurRep);

RegStrC='-R W:3:0:10.1';
RhoStrC=[' -u ' num2str(1e+0) ' '];

RegStrC='-R T:3:0:0.1';
RhoStrC=[' -u ' num2str(1e-1) ' '];

TrajPartToUseXC={1:14080 14500+(1:14080) 29100+(1:14080)};
CurI=1;


for CurI=1:3
    disp(['------ ' num2str(CurI) ' --------']);
    CurTrajPart=TrajPartToUseXC{CurI};
%     CurData=DataCCP4(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
%     CurData=DataCCP8(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurData=DataCCP10(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurnTS=17;
    CurTimePoints_ms=(linspace(CurTrajPart(1),CurTrajPart(end),CurnTS)-1)*AcqDwellTime_us/1000;
    CurTSC=exp(-1i*2*pi.*SB0r(:,:,SliI).*perm72(CurTimePoints_ms)/1000);
    CurTSB=GetTSCoeffsByLinear(numel(CurTrajPart),CurnTS);
    CurTSBP=permute(CurTSB,[3 1 4 5 6 7 2]);
    CurKerns=getKernsFromTrajM(TrajM(RepsForKerns(CurRep),CurTrajPart),Sz,CurTSB);
    RecByPart(:,:,CurI)=bart(['picsS -m ' RhoStrC ' -d 5 -S ' RegStrC ' ' ScriptFN_CompgBo],...
        FillOnesTo16(Sz),CurData,...
        SensCCS(:,:,SliI,nCCToUseX),STraj3M(:,CurTrajPart,CurRepSet),CurTSBP,CurTSC,...
        sum(CurKerns(:,:,CurRep,:,:,:,:),3),ones(1,1,1,1,1,1,CurnTS));
end
%% Fit RecByPart
WhichTSToUs_byPart=1:3;
for SliI=1:nSlices
    disp(SliI);
    [~, UpdatedB0Map_byPart(:,:,SliI), UpdatedT2SMap_ms_byPart(:,:,SliI), s_vals_byPart(:,:,:,SliI), Fitted_byPart(:,:,:,SliI), PDBase0_byPart(:,:,SliI)]=...
        FitToModel_MPBD1CSf(RecByPart,WhichTSToUs_byPart,14.080*AcqDwellTime_us,TE0_ms);
end
%% 0.000008
    CurTrajPart=14500:43180;
%     CurData=DataCCP4(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
%     CurData=DataCCP8(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurData=DataCCP10(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurnTS=17*2-1;
    CurTimePoints_ms=(linspace(CurTrajPart(1),CurTrajPart(end),CurnTS)-1)*AcqDwellTime_us/1000;
    TP29100_ms=(29100-1)*AcqDwellTime_us/1000;
    CurTSC=exp(-1i*2*pi.*SB0r(:,:,SliI).*perm72(CurTimePoints_ms)/1000);
%     CurTSCb=exp(1i.*mean(angle(DD8),3).*perm72(CurTimePoints_ms-CurTimePoints_ms(1))/((29100-14500)*AcqDwellTime_us/1000));
%     CurTSCb=exp(1i.*(angle(DD10(:,:,2))).*perm72(CurTimePoints_ms-CurTimePoints_ms(1))/((29100-14500)*AcqDwellTime_us/1000));
    CurTSCb=exp(1i.*(angle(DD10(:,:,2))).*perm72(CurTimePoints_ms>=(TP29100_ms-0.001)));
    CurTSC=CurTSC.*CurTSCb;
    CurTSB=GetTSCoeffsByLinear(numel(CurTrajPart),CurnTS);
    CurTSBP=permute(CurTSB,[3 1 4 5 6 7 2]);
    CurKerns=getKernsFromTrajM(TrajM(RepsForKerns(CurRep),CurTrajPart),Sz,CurTSB);
    RecByPart23=bart(['picsS -m ' RhoStrC ' -d 5 -S ' RegStrC ' ' ScriptFN_CompgBo],...
        FillOnesTo16(Sz),CurData,...
        SensCCS(:,:,SliI,nCCToUseX),STraj3M(:,CurTrajPart,CurRepSet),CurTSBP,CurTSC,...
        sum(CurKerns(:,:,CurRep,:,:,:,:),3),ones(1,1,1,1,1,1,CurnTS));
%%
% FitLinearPhase;
LinPhase=angle(LinFieldFromParams(BestX));
%%
    CurTrajPart=14500:43180;
    CurTrajB=ismember(CurTrajPart,union(TrajPartToUseXC{2},TrajPartToUseXC{3}(1000:end)));
%     CurData=CurTrajB.*DataCCP4(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurData=CurTrajB.*DataCCP10(1,CurTrajPart,RepsForKerns(CurRep),nCCToUseX);
    CurnTS=17*2-1;
    CurTimePoints_ms=(linspace(CurTrajPart(1),CurTrajPart(end),CurnTS)-1)*AcqDwellTime_us/1000;
    TP29100_ms=(29100-1)*AcqDwellTime_us/1000;
    CurTSC=exp(-1i*2*pi.*SB0r(:,:,SliI)*1.3.*perm72(CurTimePoints_ms)/1000);
%     CurTSCb=exp(1i.*mean(angle(DDa),3).*perm72(CurTimePoints_ms-CurTimePoints_ms(1))/((29100-14500)*AcqDwellTime_us/1000));
%     CurTSCb=exp(1i.*(angle(DDa(:,:,2))).*perm72(CurTimePoints_ms-CurTimePoints_ms(1))/((29100-14500)*AcqDwellTime_us/1000));
%     CurTSCb=exp(-1i.*LinPhase.*perm72(CurTimePoints_ms-CurTimePoints_ms(1))/((29100-14500)*AcqDwellTime_us/1000));
%     CurTSCb=exp(1i.*(angle(DD10(:,:,2))).*perm72(CurTimePoints_ms>=(TP29100_ms-0.001)));
%     CurTSC=CurTSC.*CurTSCb;
    CurTSB=GetTSCoeffsByLinear(numel(CurTrajPart),CurnTS);
    CurTSBP=CurTrajB.*permute(CurTSB,[3 1 4 5 6 7 2]);
    CurKerns=getKernsFromTrajM(TrajM(RepsForKerns(CurRep),CurTrajPart),Sz,(CurTrajB.').*CurTSB);
    RecByPart23Wb=bart(['picsS -m ' RhoStrC ' -w 0.000008 -d 5 -S ' RegStrC ' ' ScriptFN_CompgBo],...
        FillOnesTo16(Sz),CurData,...
        SensCCS(:,:,SliI,nCCToUseX),STraj3M(:,CurTrajPart,CurRepSet),CurTSBP,CurTSC,...
        sum(CurKerns(:,:,CurRep,:,:,:,:),3),ones(1,1,1,1,1,1,CurnTS));
%%
figure;plot(abs(TrajM(1,:)));hold on;plot(14080*[1 1],[0 60],'r');plot(14500*[1 1],[0 60],'r');plot(29100*[1 1],[0 60],'r');plot(28580*[1 1],[0 60],'r');
%%
%% T2* histogram
%% Using components
T2svalues_ms=linspace(5,300,200);
Decays=exp(-TimePointsMed_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');
%
CurReps=1;
nComponentsToUse=2;

ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];

TSCxPMedOnlyB0=exp(-1i*2*pi.*SB0r.*perm72(TimePoints_Med_ms)/1000);

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);

LLR_lambda=0.1;
RhoStr=[' -u ' num2str(1e-3) ' '];
BlkSz=4;
disp('Prepared LLR');
%% just 3 echos, given B0
TSBSimple=GetTSCoeffsByLinear(nTSMed,4);
TSBSimple=TSBSimple(:,1:3);
TSBSimple((nTSMed*2/3+1):end,3)=1;

% TSPerOneEcho=nTSMed/3;
% TSBSimple=zeros(nTSMed,3);
% TSBSimple(1:TSPerOneEcho,1)=1;
% TSBSimple(TSPerOneEcho+(1:TSPerOneEcho),2)=1;
% TSBSimple(TSPerOneEcho*2+(1:TSPerOneEcho),3)=1;
CompsPSimple=permute(TSBSimple,[3 4 5 6 7 2 1]);
Sz16gB0S=FillOnesTo16([Sz 1 1 1 3]);

RegStrS='-R T:3:0:10.1';
% RegStrS='-R W:3:0:0.001';
RegStrS='-R W:3:0:10.1';
% RegStrS='-R W:3:32:10.1';

% RhoStrS=[' -u ' num2str(1e-3) ' '];
RhoStrS=[' -u ' num2str(1e-1) ' '];
RhoStrS=[' -u ' num2str(1e+0) ' '];

SliI=4;
CurRep=1:3;
CurRepSet=RepsForKerns(CurRep);

% CurData=permute(DataCCP1S(TrajPartMed,nCCToUseX,SliI),[4 1 3 2]);
CurData=DataCCP4(1,TrajPartMed,RepsForKerns,nCCToUseX);
CurTSCxPMedOnlyB0=TSCxPMedOnlyB0(:,:,SliI,:,:,:,:);
Rec_gB0=bart(['picsS -m ' RhoStrS ' -d 5 -S ' RegStrS ' ' ScriptFN_CompgBo],...
    Sz16gB0S,CurData,...
    SensCCS(:,:,SliI,nCCToUseX),STraj3M(:,:,CurRepSet),TSBPMed,CurTSCxPMedOnlyB0,...
    sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsPSimple);

Rec_gB0X=squeeze(sum(Rec_gB0.*CompsPSimple,6));
    
    
%%
SliI=4;
for SliI=1:nSlices
    disp(SliI);
%     CurData=NoNavB.*permute(DataCCP1S(TrajPartMed,nCCToUseX,SliI),[4 1 3 2]);
    CurData=permute(DataCCP1S(TrajPartMed,nCCToUseX,SliI),[4 1 3 2]);
    CurTSCxPMedOnlyB0=TSCxPMedOnlyB0(:,:,SliI,:,:,:,:);
%     RegStr='-R T:3:0:0.000001';
    RegStr='-R W:3:0:0.001';
    CurRep=1;
    Rec_CompgB0_C=bart(['picsS -m ' RhoStr ' -d 5 -S -b ' num2str(BlkSz) ' ' RegStr ' ' ScriptFN_CompgBo],...
        Sz16CompgB0,CurData,...
        SensCCS(:,:,SliI,nCCToUseX),STraj3M(:,:,CurRep),TSBPMed,CurTSCxPMedOnlyB0,...
        sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP);
%         SensCCS(:,:,SliI,nCCToUseX),STraj3MMed(:,:,CurRep),NoNavB.*TSBPMed,CurTSCxPMedOnlyB0,...
    
    Rec_CompgB0_MX=squeeze(sum(Rec_CompgB0_C.*CompsP,6).*CurTSCxPMedOnlyB0);
    Rec_CompgB0_MXS(:,:,:,SliI)=Rec_CompgB0_MX.*CurTSCxPMedOnlyB0;
end
%% 
save([mainP filesep 'Rec_CompgB0_MXS.mat'],'Rec_CompgB0_MXS');
%% Fit Subspace
WhichTSToUsSb=23:71;
for SliI=1:nSlices
    disp(SliI);
    [~, UpdatedB0MapSb(:,:,SliI), UpdatedT2SMap_msSb(:,:,SliI), s_valsSb(:,:,:,SliI), FittedSb(:,:,:,SliI), PDBase0Sb(:,:,SliI)]=...
        FitToModel_MPBD1CSf(Rec_CompgB0_MXS(:,:,:,SliI),WhichTSToUsSb,dT_Med_ms,TE0_ms);
end
%%

Both2=cat(5,RecCMS(:,:,:,ROrd),Rec_CompgB0_MXS4(:,:,[2 8 15],ROrd),FittedSb(:,:,[2 8 15],ROrd),imresize(THLRMultiShot_RSS(:,:,[2 5 10],ROrd),Sz));
% Both2=cat(5,RecCMS(:,:,:,ROrd),Rec_CompgB0_MXS4(:,:,[2 8 15],ROrd),imresize(THLRMultiShot_RSS(:,:,[2 5 10],ROrd),Sz));
% Both2=Both2./grms(Both2,1:3);

fgmontagex(perm43(CombineDims(Both2(:,:,:,[2 7],:),[3 2])));caxis(caxis/1.8)

fgmontagex(perm43(CombineDims(Both2(:,:,:,[3 8],:),[3 2])));caxis(caxis/1.8)

%% LLR?
% CurRep=1;
% Rec_CompgB0_C=bart(['picsS -m ' RhoStr ' -d 5 -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
%     Sz16CompgB0,NoNavB.*DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
%     SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),NoNavB.*TSBPMed,TSCxPMedOnlyB0,...
%     sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP);
%%
Rec_CompgB0_M=cat(4,Rec_CompgB0_C{:});
Rec_CompgB0_MX=squeeze(sum(Rec_CompgB0_M.*CompsP,6));
%% Split prox
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);

nCS=1;

nEchos=numel(TimePointsMed_ms);

PrefixBase='SKEPTIC_';
ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
LS_ScriptFN=ScriptFN_AllTS;

OneFN='/autofs/cluster/kawin/Gilad/One';
writecfl(OneFN,1);

ImSz16=FillOnesTo16(Sz);
ImSz16(TS_Dim)=nEchos;
ImSz16FN=[mainP filesep 'ImSz16'];
writecfl(ImSz16FN,ImSz16);

disp('Prepared folders');
%%
TEs_ms=TE0_ms+TimePointsMed_ms.';

NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms=TEs_ms.';
EchoTimes_ms3=permute32(EchoTimes_ms);
%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
ElemL={M,P,B,TT};
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};

disp('ok operators');
%%
nccToUse=23;
WhichRepsToUse=1;

T2SRange=[4 200];

RepSets={1};

clear BARTCmd
%%
ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
ElemsLambda=[1e-4,1e-4,1e-9,1e-10];
OutPrefix='Lm4c_';
OutPrefixBase='gB0xSli';

TSBFN=[mainP filesep 'TSB'];
writecfl(TSBFN,TSBPMed);

ElemsAlphas=[1e-4, 1e-0, 1e+4, 1e+6];
ElemsLambda=[1e-6,1e-3,1e-9,1e-10];
OutPrefix='Lm6c_';

ElemsAlphas=[1e-4, 1e-0, 1e+4, 1e+6];
ElemsLambda=[1e-5,1e-3,1e-9,1e-10];
OutPrefix='Lm5c_';

ElemsAlphas=[1e-4, 1e-1, 1e+4, 1e+6];
ElemsLambda=[1e-5,1e-3,1e-9,1e-9];
OutPrefix='Lm5d_';

ElemsAlphas=[1e-4, 1e-1, 1e+5, 1e+7];
ElemsLambda=[1e-5,1e-3,1e-9,1e-10];
OutPrefix='Lm5d_';

% ElemsAlphas=[1e-4, 1e-1, 1e+5, 1e+7];
% ElemsLambda=[1e-5,1e-3,1e-8,1e-10];
% OutPrefix='Lm5d_';

for i=1:numel(ElemsL)
    writecfl([mainP filesep 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[Sz 1 1 1 1 1]));
end

writecfl([mainP filesep 'ElemsAlpha'],ElemsAlphas.');
writecfl([mainP filesep 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([mainP filesep 'ElementTypes'],ElementTypes.');

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
end
% ninneriterBART(2)=0;

writecfl([mainP filesep 'ninneriter'],ninneriterBART);
disp('Wrote general splitProx files');
%%
SliI=6;
for SliI=4 % 1:nSlices
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(CurSPrefix);
    
    SensCC=SensCCS(:,:,SliI,1:nccToUse);
    DataCCP=permute(DataCCP1S(:,1:nccToUse,SliI),[3 1 4 2]);
    
    SensFN=[CurSPrefix 'Sens'];
    writecfl(SensFN,SensCC);
    
for rs=1 % 1:numel(RepSets)
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
        CurSRPrefix=[mainP filesep OutPrefixBase num2str(SliI) '_R' RepsStr '_'];
        
        CurRPrefix=[mainP filesep 'R' RepsStr '_'];       
        TrajFN=[CurRPrefix 'Traj'];
        writecfl(TrajFN,STraj3M(:,:,WhichRepsToUse));    
        sumKernsFN=[CurRPrefix 'sumKerns'];
        writecfl(sumKernsFN,sum(KernsPMMed(:,:,WhichRepsToUse,:,:,:,:),3));

% c0=PDBase0Sb(:,:,SliI);
c0=PDBase0_byPart(:,:,SliI);

% W=s_valsSb(:,:,1,SliI);
% T2S=UpdatedT2SMap_msSb(:,:,SliI);
% T2S=min(max(abs(T2S),T2SRange(1)),T2SRange(2));
% ST2SSb(:,:,SliI)=SmoothByW(T2S,W,SRange,SSigs);
% sST2SSb=SmoothBySlices(ST2SSb(:,:,SliI),[20 20],3);
% t0=ST2SSb(:,:,SliI);
% t0=sST2SSb;
t0=UpdatedT2SMap_ms_byPart(:,:,SliI);
% t0=UpdatedT2SMap_msSb(:,:,SliI)*0+50;
t0=min(max(abs(t0),T2SRange(1)),T2SRange(2));
b0=SB0r(:,:,SliI);
% b0=UpdatedB0Map_byPart(:,:,SliI);
m0=abs(c0);
m0=min(m0,median(m0(:))*M0MedianFac);
% p0=angle(c0);
p0=angle(PDBase0THLRS(:,:,SliI));
p0=angle(PDBase0_byPart(:,:,SliI));

Elem0={m0,p0,b0,t0};
for i=1:numel(Elem0)
    writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end

BARTS_Aopx.ImSz16=ImSz16FN;
BARTS_Aop=BARTS_Aopx;
BARTS_Aop.Others{1}=SensFN;
BARTS_Aop.Others{2}=TrajFN;
BARTS_Aop.Others{3}=TSBFN;
BARTS_Aop.Others{4}=OneFN;
BARTS_Aop.Others{5}=sumKernsFN;

SigToUse=DataCCP(:,:,WhichRepsToUse,:);
ksp_adjFN=[CurSRPrefix 'sig_adj'];
% if(~exist([ksp_adjFN '.cfl'],'file'))
    ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
    disp('got ksp_adj');
    writecfl(ksp_adjFN,ksp_adj);
% end

Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTt0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
Rec0X=squeeze(sum(Rec0,CS_Dim));
% AHA_Rec0=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Rec0,BARTS_Aop.Others{:});
        
disp(['saved all ' CurSRPrefix]);

for i=1:4
    delete([CurSRPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSRPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue?
ContinueRun=false;
if(ContinueRun)
    for i=1:numel(Elem0)
        writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end

CurSRPrefixOut=[CurSRPrefix OutPrefix];

% disp('Wrote maps');
BARTCmd{SliI,rs}=bartCmd(['splitProx -i 10000 -l 1.1 -I 1000 -s 60 -d 2 -g -f -S ' ksp_adjFN ' -W ' CurSRPrefix ' -F ' mainP filesep ' -O ' CurSRPrefixOut ' ' LS_ScriptFN],BARTS_Aop.ImSz16,OneFN,BARTS_Aop.Others{:});
% system(BARTCmd{SliI,rs});
end
end
disp('Prepared BARTCmd');
%%
% SliI=2;
% SliI=3;
% SliI=14;

% CurSRPrefix=[mainP filesep 'gB0Sli' num2str(SliI) '_R' RepsStr '_'];
% CurSRPrefix=[mainP filesep OutPrefixBase num2str(SliI) '_R' RepsStr '_'];

% CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];
% CurSRPrefixOut=[CurSRPrefix OutPrefix];


ErrVec=readcfl([CurSRPrefixOut 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
% figure;plot(ErrVec)

for i=1:4
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter1000']);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));

disp('Loaded maps');
%%
ElemRanges={[0 11e-4],[-pi pi],[-300 300],[0 200]};
mThresh=60e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); xlabel([num2str(numel(ErrVec)/sum(ninneriterBART)) ' : ' num2str(ErrVec(end),'%.7g')]); end
    if(i==3), xlabel(num2str(ElemsAlphas,' %.9g,')); end
    if(i==4), xlabel(num2str(ElemsLambda,' %.9g,')); end
end
%%
AllRecs=cat(4,RecByPart(:,:,[1 2 3 3]),RecMX(:,:,[2 7 14 19]),imresize(THLRMultiShot_RSS(:,:,[2 6 10 14],SliI),Sz));