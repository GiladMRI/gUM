ScanP='/autofs/space/daisy_001/users/data/Gilad/SEPTI_2Echo_Try1_1Dec19/';
BaseFN='meas_MID00364_FID52946_gep2d_327_T19_6Sli_9Rep_2Echo';
BaseFN='meas_MID00360_FID52942_gep2d_327_T18_20Sli_50Rep_2Echo';
% RefFldMapP=[ScanP 'meas_MID04675_FID21637_gre_te4_9' filesep];

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
G11mm36ms=load('Grads11mm36ms.mat');
AllGrads1136=cat(1,G11mm36ms.Grads{:}).';
AllGrads1136(4800,1)=0;
AllGrads=cat(2,AllGrads,AllGrads11,zeros(4800,1),AllGrads1136);

nInnerShotsAll=[8 6 4 3 2 8 6 4 3 2 3 3 6 4 3 2 0 5 4 3];
nInnerShots=nInnerShotsAll(TrajType);
% GTrajaCBase=G2mm.Grads{TrajType};
GTrajaCBase=AllGrads(:,TrajType).';
TrajGradnSamples=3600;
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
%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;

disp(['Traj type ' num2str(TrajType)]);
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

GradReduceFac=100/floor(ParamA/100);
if(~isfinite(GradReduceFac))
    GradReduceFac=1;
end
PhiRotPerRep=110; % WipMemBlock.adFree{8};
AcqDwellTime_us=2; %WipMemBlock.adFree{13}/1000;

dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;

CurzLocs=RotatedLocs(3,:);

Thickness_mm=ADatax.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
DistBetweenSlices_mm=CurzLocs(2)-CurzLocs(1);
%%
RLocs=load([RefFldMapP 'Locs.mat']);
sTwix=load([RefFldMapP 'sTwixX.mat']);
RefLocs=RLocs.RotatedLocs;
RefzLocs=RefLocs(3,:);

RefMaps=load([RefFldMapP 'B0T2S.mat']);
RefSens=load([RefFldMapP 'Sens.mat']);

RefSV1=grmss(RefMaps.s_vals,3);

% RefSV1x=padarray(RefSV1(:,:,Ord(1)),[0 3],'Both');
% RefSV1x=padarray(RefSV1(:,:,Ord(1)),[0 6],'Post');
% RefSV1x=padarray(padarray(RefSV1(:,:,Ord(1)),[0 4],'Post'),[0 2],'Pre');
% RefSV1x=padarray(padarray(RefSV1(:,:,Ord(1)),[0 5],'Post'),[0 1],'Pre');
% RefSV1x=padarray(padarray(RefSV1(:,:,:),[0 5],'Post'),[0 1],'Pre');

RefSV1x=padarray(RefSV1(:,:,:),[0 5],'Post');
RefSV1x=rot90(RefSV1x);
RefSV1x=circshift(RefSV1x,-2,1);
RefSV1x(:,:,[2 4 6 8])=(RefSV1x(:,:,[2 4 6 8]-1)+RefSV1x(:,:,[2 4 6 8]+1))/2;
RefSV1x=perm31(interp1(RefzLocs,perm31(RefSV1x),CurzLocs,'linear','extrap'));

% tmp=RefSV1x;
% RefSV1x(:,:,1:2:100)=tmp;
% RefSV1x(:,:,2:2:99)=(tmp(:,:,1:end-1)+tmp(:,:,2:end))/2;
% RefSV1x(:,:,100)=tmp(:,:,end);
% figure;imshowpair(grmss(Rec1p,8)/grmss(Rec1p),RefSV1x/grmss(RefSV1x));
% figure;imshowpair(grmss(Rec1p,8)/gmax(grmss(Rec1p,8)),RefSV1x/gmax(RefSV1x),'checkerboard');

RefSVsx=padarray(RefMaps.s_vals,[0 5],'Post');
RefSVsx=rot90(RefSVsx);
RefSVsx=circshift(RefSVsx,-2,1);
RefSVsx(:,:,:,[2 4 6 8])=(RefSVsx(:,:,:,[2 4 6 8]-1)+RefSVsx(:,:,:,[2 4 6 8]+1))/2;
RefSVsx=perm41(interp1(RefzLocs,perm41(RefSVsx),CurzLocs,'linear','extrap'));

% RefSensx=padarray(padarray(RefSens.SensB(:,:,:,:,1),[0 5],'Post'),[0 1],'Pre');
RefSensx=padarray(RefSens.SensB(:,:,:,:,1),[0 5],'Post');
RefSensx=rot90(RefSensx);
RefSensx=circshift(RefSensx,-2,1);
RefSensx(:,:,:,[2 4 6 8])=(RefSensx(:,:,:,[2 4 6 8]-1)+RefSensx(:,:,:,[2 4 6 8]+1))/2;
RefSensx=perm41(interp1(RefzLocs,perm41(RefSensx),CurzLocs,'linear','extrap'));

% RefB0x=padarray(padarray(RefMaps.B0M_Hz,[0 5],'Post'),[0 1],'Pre');
% RefB0x=padarray(padarray(RefMaps.UpdatedB0Map_Hz,[0 5],'Post'),[0 1],'Pre');
RefB0x=padarray(RefMaps.UpdatedB0Map_Hz,[0 5],'Post');
% RefB0x=padarray(RefMaps.B0M_Hz,[0 5],'Post');
RefB0x=rot90(RefB0x);
RefB0x=circshift(RefB0x,-2,1);
RefB0x(:,:,[2 4 6 8])=(RefB0x(:,:,[2 4 6 8]-1)+RefB0x(:,:,[2 4 6 8]+1))/2;
RefB0x=perm31(interp1(RefzLocs,perm31(RefB0x),CurzLocs,'linear','extrap'));

% RefB0x=imresize(RefB0x,Sz); % X Y Ch SliI
RefMsk=RefSV1x>1e-4;
SE = strel('disk',3,4);
RefMsk=imdilate(RefMsk,SE);
% RefMsk=SmoothBySlices(RefMsk,[21 21],7.6)>1e-9;
RefMsk=imfillholesBySlices(RefMsk);

% RefFatMaskx=squeeze(abs(RefSVsx(:,:,2,:)))>5e-5;
RefFatMaskx=squeeze(abs(RefSVsx(:,:,2,:)))>10e-5;

AllB0_EME=load([RefFldMapP 'AllB0_EME.mat']);
AllB0_EME=AllB0_EME.AllB0_EME;
RefB0EMEx=padarray(perm43(AllB0_EME),[0 5],'Post');
RefB0EMEx=rot90(RefB0EMEx);
RefB0EMEx=circshift(RefB0EMEx,-2,1);
RefB0EMEx(:,:,:,[2 4 6 8])=(RefB0EMEx(:,:,:,[2 4 6 8]-1)+RefB0EMEx(:,:,:,[2 4 6 8]+1))/2;
RefB0EMEx=perm41(interp1(RefzLocs,perm41(RefB0EMEx),CurzLocs,'linear','extrap'));
disp('Read ref');
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
    AcqTimePoints_us=Extra*GradDwellTime_us+Delay_us+(0:AcqDwellTime_us:(TrajGradnSamples*GradDwellTime_us-1));

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
RefB0x=imresize(RefB0x,Sz);
RefSensx=imresize(RefSensx,Sz); % X Y Ch SliI
RefSVsx=imresize(RefSVsx,Sz); % X Y Ch SliI
RefSV1x=imresize(RefSV1x,Sz);
RefSV1x=max(0,RefSV1x);
RefB0EMEx=imresize(RefB0EMEx,Sz); % X Y Ch SliI
RefB0EMEx=perm43(RefB0EMEx);

B0MapToUse_Hz=RefB0EMEx(:,:,:,4);

dB0dx=symD(B0MapToUse_Hz,1);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dz=symD(B0MapToUse_Hz,3)*Thickness_mm/DistBetweenSlices_mm;
%%
% image_datan = AData{1}.noise.unsorted(); % no slicing supported atm
% image_data = ADatax.image.unsorted(); % no slicing supported atm
% XX=PartitionDim(PartitionDim(image_data,3,nRepsHdr*nSlices),4,nRepsHdr);
% XX=CombineDims(XX,[3 1]);
%%
nCh=numel(ADatax.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList);
ncc=31;

SelfSens1S=zeros([Sz nCh nSlices]);
SensCCS=zeros([Sz nSlices ncc]);
Rec1xS=zeros([Sz nSlices]);
Rec1pRMSS=zeros([Sz nSlices]);
Rec1pS=zeros([Sz nCh nSlices]);

sccmtxS=zeros([nCh nCh nSlices]);
disp('Prepared mem for start');
%%
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 nCh]);

SlicesToRead=1;
clear KernsP_TSTHLRBase KernsPMMedBase
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    disp(['Slice ' num2str(SlicesToRead)]);
RepsToRead=1:(nRepsHdr);
% tmp=ADatax.image(); % SamplesInADC Ch 1 1 Slices 1 1 3 Reps 1 ADCs
% ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,1,RepsToRead,:,:,:,:,:,:,:,:);
WhichEcho=3;
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,WhichEcho,RepsToRead,:,:,:,:,:,:,:,:);
ADataIsL=squeeze(ADataIsL);
ADataIsL=CombineDims(ADataIsL,[4 1]);
disp('Read data');
%
% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
% TrajPartToUse=0+(1:9000);
TrajPartToUse=0+(1:2000);
RepsToUse=1:(nRepsHdr-1);

RepsToUse=8:49;

Rots=mod(360-(Ord-1)*50 + -(WhichEcho-1)*240,360);
CurRot=Rots(SlicesToRead);
RotsX=mod(PhiRotPerRep*(1:36),360);
ShiftIdx=find(CurRot==RotsX);
RepsI=mod( (1:36) +ShiftIdx -1,36)+1;
% RepsmIx=mod( RepsToUse-ShiftIdx -1,36)+1;
RepsmIx=mod( (1:nRepsHdr)-ShiftIdx -1,36)+1;
% RepsI=mod( (1:36) +8 -1,36)+1;
% RepsI=mod( (1:36) +22 -1,36)+1;
% DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsI,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,RepsI);
% DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsI,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,1:36);
DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,RepsmIx);
% Rec1p=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsI),(DataPC(:,TrajPartToUse,RepsI,1,1,1,1,:)),OnesSensC);
Rec1p=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsmIx(RepsToUse)),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
%

Rec1pRMS=grmss(Rec1p,8);
Rec1pRMSS(:,:,SlicesToRead)=Rec1pRMS;
Rec1pS(:,:,:,SlicesToRead)=squeeze(Rec1p);
fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
end
%%
%%
SlicesToRead=1;
clear KernsP_TSTHLRBase KernsPMMedBase
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    disp(['Slice ' num2str(SlicesToRead)]);
    RepsToRead=1:(nRepsHdr);
    ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepsToRead,:,:,:,:,:,:,:,:);
    ADataIsL=squeeze(ADataIsL);
    ADataIsL=CombineDims(ADataIsL,[4 1]);
    disp('Read data');
%
ChRMS=grmss(ADataIsL,[1 3]);
[SChRMS ChOrd]=sort(ChRMS,'descend');

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

sccmtxS(:,:,SlicesToRead)=sccmtx;

% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
TrajPartToUse=0+(1:2000);
RepsToUse=5:47;

Rots=mod(360-mod((Ord-1)*50,360),360);
CurRot=Rots(SlicesToRead);
RotsX=mod(110*(1:36),360);
ShiftIdx=find(CurRot==RotsX);
% RepsI=mod( (1:36) +ShiftIdx -1,36)+1;
% DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsI,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,1:36);
% Rec1p=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsI),(DataPC(:,TrajPartToUse,RepsI,1,1,1,1,:)),OnesSensC);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
% 
% RepsmI=mod( (1:36) -ShiftIdx -1,36)+1;
% DataPC=permute(ADataIsL(1:nTrajToUse,:,1:36,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,RepsmI);
% Rec1p=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsmI),(DataPC(:,TrajPartToUse,1:36,1,1,1,1,:)),OnesSensC);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);
RepsmI=mod( RepsToUse -ShiftIdx -1,36)+1;
RepsmIx=mod( (1:nRepsHdr)-ShiftIdx -1,36)+1;
DataPC=permute(ADataIsL(1:nTrajToUse,:,RepsToUse,:,:,:,:,:,:),[8 1 3 7 6 5 4 2]).*modx(:,:,RepsmI);
Rec1p=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsmI),(DataPC(:,TrajPartToUse,:,1,1,1,1,:)),OnesSensC);
% fgmontagex(grmss(Rec1p,8));title(['Slice #' num2str(SlicesToRead) ', ' num2str(ROrd(SlicesToRead))]);

Rec1pRMS=grmss(Rec1p,8);
Rec1pRMSS(:,:,SlicesToRead)=Rec1pRMS;
Rec1pS(:,:,:,SlicesToRead)=squeeze(Rec1p);


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
%     SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p),8,TrgSz);
    FF=squeeze(bart('fft 6',permute(squeeze(Rec1p),[4 1 2 3])));
    FF=permute(FF,[4 1 2 3]);
    SelfSens = squeeze(bart(['ecalib -r 24 -c 0.5 '], FF));
end
SelfSens1=SelfSens(:,:,:,1);
% SelfSens1=SelfSens1.*DMsk;
disp('ok SelfSens');

SelfSens1S(:,:,:,SlicesToRead)=SelfSens1;
% Recon each channel separately, CC
% DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;

disp('ok rec pre channel cc');
%
SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
disp('ok SensCC');

% DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,RepsmIx);
DataCCP=permute(DataPC,[1:3 8 4:7]);
disp('ok 0');

Rec1x=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsmIx(RepsToUse)),perm84(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),SensCC);

Rec1xS(:,:,SlicesToRead)=Rec1x;
SensCCS(:,:,SlicesToRead,:)=SensCC;

% using ref sens
RefSensCC=permute(sum(RefSensx(:,:,:,Ord(SliI)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
RefSensCC=permute43(RefSensCC);
% RefMsk(:,:,Ord(SliI)).*
Rec1xR=bart('pics -S -t ',BARTTrajP(:,TrajPartToUse,RepsmIx(RepsToUse)),perm84(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),RefSensCC);

Rec1xSR(:,:,SlicesToRead)=Rec1xR;
RefSensCCS(:,:,SlicesToRead,:)=RefSensCC;
end
%%
save([mainP filesep 'TrajM' '.mat'],'TrajM');
save([mainP filesep 'sccmtxS' '.mat'],'sccmtxS');
save([mainP filesep 'Rec1pRMSS' '.mat'],'Rec1pRMSS');
save([mainP filesep 'Rec1xS' '.mat'],'Rec1xS');
save([mainP filesep 'Rec1pS' '.mat'],'Rec1pS');
save([mainP filesep 'SelfSens1S' '.mat'],'SelfSens1S');
save([mainP filesep 'SensCCS' '.mat'],'SensCCS');

SmallVarT=50e6;
getSmallVars;
save([mainP filesep 'SmallVars.mat'],SmallVars{:});
%%
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    disp(['Slice ' num2str(SlicesToRead)]);
    sccmtx=sccmtxS(:,:,SlicesToRead);
    RefSensCC=permute(sum(RefSensx(:,:,:,Ord(SliI)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    RefSensCC=permute43(RefSensCC);
    RefSensCCS(:,:,SlicesToRead,:)=RefSensCC;
end
%%
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
%
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
%
clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
%     disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
end
%
THLR_lambda=10;
RhoStr=[' -u ' num2str(1e-3) ' '];
%%
save([mainP filesep 'KernsPMMedBase' '.mat'],'KernsPMMedBase');
save([mainP filesep 'KernsP_TSTHLRBase' '.mat'],'KernsP_TSTHLRBase');
%%
RepsSets={6:45, [1 13 25]+7};
if(HighRes)
    RepsSets={6:45, [36    21    15    11    40    14    19    25    37]};
end

nRepsSets=numel(RepsSets);
THLRMultiShot_RSS=zeros([Sz nTS_THLR nSlices]);
%%
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    disp([datestr(now) ' Slice ' num2str(SlicesToRead)]);
    
    sccmtx=sccmtxS(:,:,SlicesToRead);

    RepsToRead=1:(nRepsHdr);
    ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepsToRead,:,:,:,:,:,:,:,:);
    ADataIsL=squeeze(ADataIsL);
    ADataIsL=CombineDims(ADataIsL,[4 1]);
    
    ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
    for i=1:ncc
        ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
    end
    DataC=permute(ADataIsLCC,[1 3 2]);
    
    Rots=mod(360-mod((Ord-1)*50,360),360);
    CurRot=Rots(SlicesToRead);
    RotsX=mod(110*(1:36),360);
    ShiftIdx=find(CurRot==RotsX);
    ShiftIdxS(SliI)=ShiftIdx;
    RepsmIx=mod( (1:nRepsHdr)-ShiftIdx -1,36)+1;
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,RepsmIx);
%     DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    
    SensCC=SensCCS(:,:,SlicesToRead,:);
    SensCC=RefSensCCS(:,:,SlicesToRead,:);
    disp('Read data');

    clear THLRMultiShot_RS
    for rs=1 % :nRepsSets
        CurReps=RepsSets{rs};
        disp(['RepSet ' num2str(rs) ' ! ' datestr(now)]);
        tmp=bart(['picsS -w 1 -d 5 -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16_THLR,NoNavB.*DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
            SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,RepsmIx(CurReps)),TSB_THLRP,1,...
            sum(KernsP_TSTHLR(:,:,RepsmIx(CurReps),:,:,:,:),3));
%             sum(KernsP_TSTHLR(:,:,CurReps,:,:,:,:),3));
        THLRMultiShot_RS(:,:,:,rs)=squeeze(tmp);
    end
    THLRMultiShot_RSS(:,:,:,SliI)=THLRMultiShot_RS;
end % THLR on all slices
%%
save([mainP filesep 'THLRMultiShot_RSS' '.mat'],'THLRMultiShot_RSS');
%%
WhichTSToUs=2:(nTS_THLR-1);
[PDBase_THLRS, UpdatedB0Map_THLRS, UpdatedT2SMap_ms_THLRS, s_vals_THLRS, Fitted0_THLRS, PDBase0_THLRS]=...
        FitToModel_MPBD1CSf(THLRMultiShot_RSS,WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
%%
save([mainP filesep 'THLR_Fit' '.mat'],'PDBase_THLRS','UpdatedB0Map_THLRS','UpdatedT2SMap_ms_THLRS','s_vals_THLRS','PDBase0_THLRS');
disp(['Saved ' mainP filesep 'THLR_Fit' '.mat']);
%%
UpdatedB0Map_RSxS=SmoothByW(UpdatedB0Map_THLRS,grmss(s_vals_THLRS,3));
RefB0xS=SmoothByW(RefB0x(:,:,Ord),RefSV1x(:,:,Ord));
RefB0xS=min(300,max(-300,RefB0xS));
disp('Prepared sB0');
%%
for SlicesToRead=1:nSlices
    SliI=SlicesToRead;
    
    Rots=mod(360-mod((Ord-1)*50,360),360);
    CurRot=Rots(SlicesToRead);
    RotsX=mod(110*(1:36),360);
    ShiftIdx=find(CurRot==RotsX);
    ShiftIdxS(SliI)=ShiftIdx;
    
%     ShiftIdx=ShiftIdxS(SliI);
    RepsmIx=mod( (1:nRepsHdr)-ShiftIdx -1,36)+1;
    RepsmIxS(SliI,:)=RepsmIx;
end
disp('ok RepsmIxS');
%%
SliIs=ROrd;
SliIs=ROrd(30:35);
nSliIs=numel(SliIs);
%%
if(HighRes)
    RSz=Sz;
else
    RSz=Sz*2;
end
% RefB0xSr=imresize(RefB0xS,RSz);
RefB0xSr=imresize(min(500,max(-500,RefB0EMEx(:,:,Ord,4))),RSz);
dB0dzr=imresize(dB0dz(:,:,Ord),RSz);

nComponentsToUse=3;
nFatComponentsToUse=0;
nTtlComponentsToUse=nComponentsToUse+nFatComponentsToUse;

disp('Prepared for LLR RepSets');

LLR_lambda=0.1;
% RhoStr=[' -u ' num2str(1e-3) ' '];
RhoStr=[' -u ' num2str(1e-1) ' '];
BlkSz=4;
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
Sz16CompgB0S=FillOnesTo16([RSz 1 1 1 nTtlComponentsToUse 1 1 numel(SliIs)]);
%%
nTSMed=17;
%%
[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nTraj,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
% TSBMed(nPointsMed,1)=0;

T2svalues_ms=linspace(5,300,200);
Decays=exp(-TimePoints_Med_ms./(T2svalues_ms.'));
[Ud,Sd,Vd]=svd(Decays,'econ');

FatT2svalues_ms=linspace(3,20,20);
FatDecays=exp(-TimePoints_Med_ms./(FatT2svalues_ms.'));
[FatUd,FatSd,FatVd]=svd(FatDecays,'econ');
FatCS_kHz=-440/1000;
FatVdx=FatVd.*exp(1i*2*pi*(TimePoints_Med_ms.')*FatCS_kHz);
% BothCmps=cat(2,Vd(:,1:nComponentsToUse),FatVdx(:,1:nFatComponentsToUse));

clear CompsP CompsPW CompsPF CompsPWM CompsPFM
% CompsP(1,1,1,1,1,:,:)=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
CompsPW(1,1,1,1,1,:,:)=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
if(nFatComponentsToUse>0)
    CompsPF(1,1,1,1,1,:,:)=permute(FatVdx(:,1:nFatComponentsToUse),[7:-1:3 2 1]);
    CompsPWM=ones(RSz).*CompsPW;
    CurRefFatMsk=double(imresize(RefFatMaskx(:,:,Ord(SliIs)),RSz)>0);
    CurRefFatMsk=ones(RSz);
    % fgmontagex(CurRefFatMsk)
    % fgmontagex(Rec1pRMSS(:,:,SliIs))
    CompsPFM=CurRefFatMsk.*CompsPF;
    CompsP=cat(6,CompsPWM,CompsPFM);
else
    CompsP=CompsPW;
end
% CompsP(1,1,1,1,1,:,:)=permute(BothCmps,[7:-1:3 2 1]);

TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

RecalcKerns=false;
if(~exist('KernsPMMedBase','var'))
    RecalcKerns=true;
else
    RecalcKerns=any(gsize(KernsPMMedBase,1:2)~=RSz*2) || size(KernsPMMedBase,7)~=nTSMed;
end
if(RecalcKerns )
    KernsPMMedBase_FN=['/autofs/cluster/kawin/Gilad/KernsPMMedBase_' num2str(TrajType) '_' num2str(nTSMed) '_m' num2str(RSz(1)) 'x' num2str(RSz(2)) '.mat'];
    if(exist(KernsPMMedBase_FN,'file'))
        disp(['loading ' KernsPMMedBase_FN '...']);
        load(KernsPMMedBase_FN);
        disp(['loaded ' KernsPMMedBase_FN]);
    else
        KernsPMMedBase=getKernsFromTrajM(TrajM(1:36,TrajPartMed),RSz,TSBMed);
        KernsPMMedBase=single(KernsPMMedBase);
        disp(['saving ' KernsPMMedBase_FN]);
        save(KernsPMMedBase_FN,'KernsPMMedBase','-v7.3');
    end
end
%
% KernsPMMed=KernsPMMedBase(:,:,mod(CurReps-1,36)+1,:,:,:,:,:,:,:,:);

KernsPMMedy=KernsPMMedBase(:,:,mod((1:(nRepsHdr-1))-1,36)+1,:,:,:,:,:,:,:,:);
disp('ok KernsPMMedy');
%%
CurB0=RefB0xSr(:,:,SliIs);
% CurB0(CurB0>50)=0;
% CurB0=min(CurB0,50);
% TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(UpdatedB0Map_RSxS(:,:,SliIs)).*TimePoints_Med_ms3/1e3));
% tmp=SmoothByW(RefB0x(:,:,Ord(SliIs)),RefSV1x(:,:,Ord(SliIs)));
% TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(RefB0xSr(:,:,SliIs)).*TimePoints_Med_ms3/1e3));
TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(CurB0).*TimePoints_Med_ms3/1e3));
% TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(tmp).*TimePoints_Med_ms3/1e3));

ThroughPlanePhaseDiffS=2*pi.*TimePoints_Med_ms3.*perm43(dB0dzr(:,:,SliIs))/1000;
ThroughPlaneDecayS=sinc(ThroughPlanePhaseDiffS/(2*pi));

TSCxPMedOnlyB0S=TSCxPMedOnlyB0S.*perm73(ThroughPlaneDecayS);
disp('ok TSCxPMedOnlyB0S');
%%
nccToUse=15;
% nccToUse=32;
% RepSetX=[1 13 25; 6 18 30; 4 16 28; 9 21 33];
% RepSetX=repmat(RepSetX,[8 1]);

% 3-shot
% BaseRepSetX=[1 13 25];
% BIdxs= [randi(3,nSliIs,1), randi(2,nSliIs,1)*2-3];
% BIdxs=mod([BIdxs(:,1) sum(BIdxs,2)],3)+1;
% BIdxs(:,3)=6-sum(BIdxs,2);
% RepSetX1=BaseRepSetX(BIdxs);
% RepSetX=mod([1 13 25]-5+(1:numel(SliIs)).'*5,36)+6;
% RepSetX=mod(RepSetX1+((1:nSliIs)-1).'*7,36)+8;
% % RepSetX=mod(BaseRepSetX+((1:nSliIs)-1).'*7,36)+8;

if(HighRes)
    BaseRepSetX=7:12;
    BaseRepSetX=7:15;
else
    BaseRepSetX=7:9;
    BaseRepSetX=[1 13 25]+22;
    BaseRepSetX=11:13;
end
RepSetX=repmat(BaseRepSetX,[nSliIs 1]);
disp('ok RepSetX');
%%
if(nccToUse==nCh)
    CurSens=perm43(RefSensx(:,:,:,Ord(SliIs)));
else
    CurSens=perm93(RefSensCCS(:,:,SliIs,1:nccToUse));
end
CurSens=imresize(CurSens,RSz);

clear CurKerns CurSTraj
disp('Preparing Traj...');
CurKerns=single(zeros([RSz*2 1 1 1 1 nTSMed 1 nSliIs]));
CurSTraj=single(zeros([3 size(STraj3MMed,2) size(RepSetX,2) 1 1 1 1 1 nSliIs]));
for i=1:nSliIs
%     CurSTraj(:,:,:,:,:,:,:,:,i)=STraj3MMed(:,:,RepSetX(i,:));
    CurSTraj(:,:,:,:,:,:,:,:,i)=STraj3MMed(:,:,RepsmIxS(SliIs(i),RepSetX(i,:)));
%     CurKerns(:,:,:,:,:,:,:,1,i)=sum(KernsPMMedy(:,:, RepSetX(i,:) ,:,:,:,:),3);
    CurKerns(:,:,:,:,:,:,:,1,i)=sum(KernsPMMedy(:,:, RepsmIxS(SliIs(i),RepSetX(i,:)) ,:,:,:,:),3);
%     CurData(1,:,:,:,1,1,1,1,i)=DataCCPS(SliIs(i),:,:,RepSetX(i,:),1:nccToUse);
end
disp('Prepared Traj');
%%
disp('Preparing data...');
CurData=single(zeros([1,nTrajToUse,size(RepSetX,2),nccToUse,1,1,1,1,nSliIs]));
for i=1:nSliIs
    SlicesToRead=SliIs(i);
    SliI=SlicesToRead;
    dispstat(['Slice ' num2str(i)],'timestamp');
%     disp([datestr(now) ' Slice ' num2str(i)]);
    
    sccmtx=sccmtxS(:,:,SlicesToRead);

    ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepSetX(i,:),:,:,:,:,:,:,:,:);
%     ADataIsL=squeeze(ADataIsL);
    ADataIsL=permute(ADataIsL,[1 2 9 11 3:8 10 12:20]);
    ADataIsL=CombineDims(ADataIsL,[4 1]);
    
    if(nccToUse==nCh)
        DataC=perm32(ADataIsL);
    else
        ADataIsLCC=single(zeros([size(ADataIsL,1) nccToUse size(ADataIsL,3)]));
        for c=1:nccToUse
            ADataIsLCC(:,c,:)=sum(ADataIsL.*permute(sccmtx(:,c),[3 1 4 5 6 7 8 9 2]),2);
        end
        DataC=permute(ADataIsLCC,[1 3 2]);
    end
%     DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(1,:,RepSetX(i,:));
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(1,:,RepsmIxS(SliIs(i),RepSetX(i,:)));
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    CurData(1,:,:,:,1,1,1,1,i)=DataCCP;
end
CurData=single(CurData);

% if(size(RepSetX,2)==1)
%     CurData=perm43(CurData);
% end
disp('Prepared data');
%%
% DataCCPS=CurData(:,:,:,1:nccToUse,1,1,1,1,Ord);
% %%
% CurData=DataCCPS(:,:,:,:,:,:,:,:,SliIs);
%% High res: 3-shot 1735 sec 9shot 1988.327092
%     RecStr=[' -m -S -d 5 -R T:259:0:0.1 '];
% RecStr=[' -m -S -d 5 -R W:3:256:0.1 '];
% RecStr=[' -m -S -d 5 -R W:259:0:0.0001 '];
% RecStr=[' -m -S -d 5 -R T:256:0:0.1 '];
RecStr=[' -m -S -d 4 -R T:256:0:10 '];
% RecStr=[' -m -S -d 4 -R W:3:0:0.01 '];
nccy=15;
% RecStr=[' -m -S -d 5 -R T:256:0:10 -R T:3:0:.03 '];
RecStr=[' -m -S -d 4 -R T:256:0:10 -R T:3:0:.001 '];
RecStr=[' -m -S -d 4 -R T:256:0:1 -R T:3:0:.001 '];
% RecStr=[' -m -S -d 4 '];
RecStr=[' -m -S -d 4 -R W:3:0:0.02 -i 60 -C 10 '];
% RecStr=[' -m -S -d 4 -R W:3:0:0.01 '];
% -i iter      	max. number of iterations
% -C iter      	ADMM max. CG iterations
% -q cclambda      	(cclambda)
% -u rho      	ADMM rho
% -s step      	iteration stepsize
% maxiter 30
% maxitercg 10
% cg_eps 0.001000
% ABSTOL 0.000000
% RELTOL 0.000000
% rho 0.500000
% alpha 1.600000
% tau 2.000000
% tau_max 20.000000
% mu 3.000000
Rec_CompgB0_C_RSS=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);
%
% Rec_CompgB0_C_RSS_FatP=Rec_CompgB0_C_RSS;
Rec_CompgB0_RSS_M=cat(4,Rec_CompgB0_C_RSS);
% Rec_CompgB0_RSS_MX=perm43(squeeze(sum(Rec_CompgB0_RSS_M.*gpermute(CompsP,[4 1]),6)));
Rec_CompgB0_RSS_MX=perm43(squeeze(sum(Rec_CompgB0_RSS_M.*CompsP,6)));
Rec_CompgB0_RSS_MX=Rec_CompgB0_RSS_MX.*perm43(squeeze(TSCxPMedOnlyB0S));
Rec_CompgB0_RSS_MX=single(Rec_CompgB0_RSS_MX);
%
WhichTSToUseTV=round(nTSMed.*0.15):round(nTSMed.*0.85);
% WhichTSToUseTV=1:(nTSMed-2);
% WhichTSToUseTV=15:60;
[PDBase_TV7, UpdatedB0Map_TV7, UpdatedT2SMap_ms_TV7, s_vals_TV7, Fitted0_TV7, PDBase0_TV7]=...
    FitToModel_MPBD1CSf(perm43(squeeze(Rec_CompgB0_RSS_MX)),WhichTSToUseTV,dT_Med_ms,TE0_ms+dT_Med);
%         FitToModel_MPBD1CSf(perm43(squeeze(Rec_CompgB0_RSS_MX)),WhichTSToUseTV,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
Fitted0_TV7=perm43(single(Fitted0_TV7));
% fgmontagex(UpdatedB0Map_TV7,[-20 20]);ylabel(RepSetX)
fgmontagex(UpdatedT2SMap_ms_TV7,[0 200]);colormap hot;ylabel(num2str(RepSetX,'%d '));xlabel(mainP,'Interpreter','None');title(RecStr);
title([num2str(size(CurData,3)) 'shots ' num2str(size(Rec_CompgB0_RSS_MX,3)) 'slices ' num2str(size(Rec_CompgB0_RSS_MX,4)) 'TS ' num2str(size(CompsP,6)) 'comps. mat ' num2str(gsize(Rec_CompgB0_RSS_MX,1:2),'%d ') ' ' RecStr]);

% fgmontagex(Rec_CompgB0_RSS_MX(:,:,1,11))
%%
% fgmontagex(Rec_CompgB0_RSS_MX(:,:,:,15));
% fgmontagex(Rec_CompgB0_RSS_MX(:,:,5:4:end,15));caxis(caxis*1.2);
% fgmontagex(Fitted0_TV7(:,:,2:2:end,15));
% fgmontagex(Rec_CompgB0_RSS_MX(:,:,12:7:end-10,11));
% fgmontagex(Fitted0_TV7(:,:,11:5:end-10,15));caxis(caxis/1.5);
% fgmontagex(UpdatedT2SMap_ms_TV7(:,:,3:2:end-2),[0 200]);colormap hot
Idxs=4:9:nSliIs;
% Idxs=7:3:nSliIs-9;
CurMsks=imresize(RefMsk(:,:,SliIs(Idxs)),RSz);
fgmontagex(CurMsks.*UpdatedT2SMap_ms_TV7(:,:,Idxs),[0 200]);colormap hot
xlabel(mainP,'Interpreter','None');
title([num2str(size(CurData,3)) 'shots ' num2str(size(Rec_CompgB0_RSS_MX,3)) 'slices ' num2str(size(Rec_CompgB0_RSS_MX,4)) 'TS ' num2str(size(CompsP,6)) 'comps. mat ' num2str(gsize(Rec_CompgB0_RSS_MX,1:2),'%d ') ' ' RecStr]);
%%
SmallVarT=50e6;
getSmallVars;
save([mainP filesep 'SmallVars2.mat'],SmallVars{:});
dir([mainP filesep 'SmallVars2.mat'])
%%
Ord=[2:2:nSlices 1:2:nSlices];

ShotsStr=['sh' GroupToStr( gmat2cell(num2str(BaseRepSetX.'),1))];
MatStr=['m' GroupToStr( gmat2cell(num2str(RSz.'),1))];
SliStr=['sl' num2str(Ord(SliIs(1))) '-' num2str(Ord(SliIs(end)))];
TSStr=['TS' num2str(size(Rec_CompgB0_RSS_MX,4))];
CompStr=['comp' num2str(size(CompsP,6))];

RecFN=[mainP filesep 'SubspaceRec_' SliStr '_' ShotsStr '_' MatStr '_' TSStr '_' CompStr '.mat'];
%%
save(RecFN,'RecStr','Rec_CompgB0_RSS_MX');
%%



%% Motion test
% RecStr=[' -m -S -d 2 -R W:3:0:0.01 '];
RecStrMotion=[' -m -S -d 2 -R W:3:0:0.1 '];
if(HighRes)
    RecStrMotion=[' -m -S -d 2 -R W:3:0:0.6 '];
end
nccy=15;
nSliIs=1;

RepToCalc=7:40;
Rec_CompgB0_RSS_MXCMS=cell(max(RepToCalc),nSlices);
for SliI=1:nSlices
    SliIs=SliI;
    
    CurSens=perm93(RefSensCCS(:,:,SliIs,1:nccToUse));
    CurSens=imresize(CurSens,RSz);

    TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(RefB0xSr(:,:,SliIs)).*TimePoints_Med_ms3/1e3));
for RepSetX=RepToCalc
    disp(['SliI ' num2str(SliI) ' RepSetX ' num2str(RepSetX)]);

clear CurKerns CurSTraj
CurKerns=single(zeros([RSz*2 1 1 1 1 nTSMed 1 nSliIs]));
CurSTraj=single(zeros([3 size(STraj3MMed,2) size(RepSetX,2) 1 1 1 1 1 nSliIs]));
for i=1:nSliIs
    CurSTraj(:,:,:,:,:,:,:,:,i)=STraj3MMed(:,:,RepsmIxS(SliIs(i),RepSetX(i,:)));
    CurKerns(:,:,:,:,:,:,:,1,i)=sum(KernsPMMedy(:,:, RepsmIxS(SliIs(i),RepSetX(i,:)) ,:,:,:,:),3);
end
CurData=single(zeros([1,nTrajToUse,size(RepSetX,2),nccToUse,1,1,1,1,nSliIs]));
for i=1:nSliIs
    SlicesToRead=SliIs(i);
    SliI=SlicesToRead;
    
    sccmtx=sccmtxS(:,:,SlicesToRead);

    ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepSetX(i,:),:,:,:,:,:,:,:,:);
    ADataIsL=permute(ADataIsL,[1 2 9 11 3:8 10 12:20]);
    ADataIsL=CombineDims(ADataIsL,[4 1]);
    
    ADataIsLCC=MultTensorMat1(ADataIsL,sccmtx(:,1:nccToUse));
%     ADataIsLCC=single(zeros([size(ADataIsL,1) nccToUse size(ADataIsL,3)]));
%     for c=1:nccToUse
%         ADataIsLCC(:,c,:)=sum(ADataIsL.*permute(sccmtx(:,c),[3 1 4 5 6 7 8 9 2]),2);
%     end
    DataC=permute(ADataIsLCC,[1 3 2]);
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(1,:,RepsmIxS(SliIs(i),RepSetX(i,:)));
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    CurData(1,:,:,:,1,1,1,1,i)=DataCCP;
end
CurData=single(CurData);

Rec_CompgB0_C_RSS=bart(['picsS ' RecStrMotion ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

Rec_CompgB0_RSS_M=cat(4,Rec_CompgB0_C_RSS);
Rec_CompgB0_RSS_MX=perm43(squeeze(sum(Rec_CompgB0_RSS_M.*gpermute(CompsP,[4 1]),6)));
Rec_CompgB0_RSS_MX=Rec_CompgB0_RSS_MX.*perm43(squeeze(TSCxPMedOnlyB0S));
Rec_CompgB0_RSS_MX=single(Rec_CompgB0_RSS_MX);
% fgmontagex(grmss(Rec_CompgB0_RSS_MX,3:4));title(RecStrMotion);xlabel(mainP,'Interpreter','None');

% Rec_CompgB0_RSS_MXCM{RepSetX}=Rec_CompgB0_RSS_MX;
Rec_CompgB0_RSS_MXCMS{RepSetX,SliI}=Rec_CompgB0_RSS_MX;
end
end % Slice loop
% Rec_CompgB0_RSS_MXCMX=cat(5,Rec_CompgB0_RSS_MXCM{:});
%%
clear Rec_CompgB0_RSS_MXCMX
dispstat('','init');
for SliI=1:nSlices
    for j=1:50
        dispstat(num2str([SliI j],'%d '),'timestamp');
        if(~isempty(Rec_CompgB0_RSS_MXCMS{j,SliI}))
            PerShotCalculated(SliI,j)=1;
            Rec_CompgB0_RSS_MXCMX(:,:,SliI,:,j)=squeeze(Rec_CompgB0_RSS_MXCMS{j,SliI});
        end
    end
%     Rec_CompgB0_RSS_MXCMX1{SliI}=cat(5,Rec_CompgB0_RSS_MXCMS{:,SliI});
end
% Rec_CompgB0_RSS_MXCMX(:,:,:,:,7:50)=cat(3,Rec_CompgB0_RSS_MXCMX1{:});
Rec_CompgB0_RSS_MXCMY=grmss(Rec_CompgB0_RSS_MXCMX,4); % X Y Slcs Shots
%%
Rec_CompgB0_RSS_MXCMX=cat(5,Rec_CompgB0_RSS_MXCMS{:,SliI});
Rec_CompgB0_RSS_MXCMY=grmss(Rec_CompgB0_RSS_MXCMX,4); % X Y Slcs Shots
% Rec_CompgB0_RSS_MXCMX1=Rec_CompgB0_RSS_MXCMX;
%%
save([mainP filesep 'Rec_CompgB0_RSS_MXCMY.mat'],'Rec_CompgB0_RSS_MXCMY');
%%
TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(RefB0xSr(:,:,SliI)).*TimePoints_Med_ms3/1e3));
% tmp=squeeze(Rec_CompgB0_RSS_MXCMX(:,:,SliI,:,5:8)).*squeeze(conj(TSCxPMedOnlyB0S));
% tmp=squeeze(Rec_CompgB0_RSS_MXCMX(:,:,1,:,5:8)).*squeeze(conj(TSCxPMedOnlyB0S));
tmp=squeeze(Rec_CompgB0_RSS_MXCMX(:,:,1,:,:)).*squeeze(conj(TSCxPMedOnlyB0S));
WhichTSToUset=2:16;
[PDBase_TVt, UpdatedB0Map_TVt, UpdatedT2SMap_ms_TVt, s_vals_TVt, Fitted0_TVt, PDBase0_TVt]=...
    FitToModel_MPBD1CSf(tmp,WhichTSToUset,dT_Med_ms,TE0_ms+dT_Med);
% UpdatedB0Map_TVtx=getEMEB0(tmp(:,:,WhichTSToUset,:),[6 6 16 16],dT_Med_ms);
UpdatedB0Map_TVtx2=getEMEB0(tmp(:,:,WhichTSToUset,:),[5 5 12 12],dT_Med_ms);
% UpdatedB0Map_TVtx3=getEMEB0(tmp(:,:,WhichTSToUset,:),[5 5 30 30],dT_Med_ms);
UpdatedB0Map_TVty=getB0_ExpPolyfit(tmp(:,:,WhichTSToUset,:),dT_Med_ms,3);
%% Now recon with added B0 variation estimation
tmpB0=SmoothByW(RefB0x(:,:,Ord(SliIs)),RefSV1x(:,:,Ord(SliIs)));
tmpB0=RefB0x(:,:,Ord(SliIs));
tmpB0=RefB0EMEx(:,:,Ord(SliIs),4);
tmpTSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(tmpB0).*TimePoints_Med_ms3/1e3));

Rec_CompgB0_C_RSS_NoB0V=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,tmpTSCxPMedOnlyB0S,...
    CurKerns,CompsP);
Rec_CompgB0_C_RSS_NoB0VX=squeeze(sum(Rec_CompgB0_C_RSS_NoB0V.*CompsP,6));

ksp_adj=bart(['linopScript -A ' ScriptFN_CompgBo],Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

ScriptFN_CompgB0_OnlyN=[BaseSP 'nuftCompgB0_OnlyN.txt'];

Rec_CompgB0_C_RSS_NoB0V_OnlyN=bart(['picsS -w 0.001063 ' RecStr ScriptFN_CompgB0_OnlyN],...
    Sz16CompgB0S,ksp_adj,...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

grmss(Rec_CompgB0_C_RSS_NoB0V_OnlyN-Rec_CompgB0_C_RSS_NoB0V)/grmss(Rec_CompgB0_C_RSS_NoB0V) % OK
%% get adj
CurKernsA=KernsPMMedy(:,:, RepsmIxS(SliIs(i),RepSetX(i,:)) ,:,:,:,:);
nShots=size(RepSetX,2);

clear ksp_adjA
for s=1:nShots
    tmp=bart(['linopScript -A ' ScriptFN_CompgBo],Sz16CompgB0S,CurData(1,:,s,1:nccy,1,1,1,1,:),...
        CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj(:,:,s),TSBPMed,TSCxPMedOnlyB0S,...
        CurKernsA(:,:,s,1,1,1,:),CompsP);
    ksp_adjA(:,:,s,1,1,:)=tmp;
end
ksp_adjAs=sum(ksp_adjA,3);

grmss(ksp_adj-ksp_adjAs)/grmss(ksp_adj) % OK
%%
% RefSensCC=permute(sum(RefSensx(:,:,:,Ord(SliI)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
% RefSensCC=permute43(RefSensCC);
% RefSensCCS(:,:,SlicesToRead,:)=RefSensCC;
%%
% TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(CurB0).*TimePoints_Med_ms3/1e3));
CurB0WVar=perm93(CurB0)+UpdatedB0Map_TVty(:,:,1:3);
% CurB0WVar=perm93(CurB0)+UpdatedB0Map_TVtx2(:,:,1:3);
TSCxPMedOnlyB0SWB0Var=perm43(perm73(exp(-1i*2*pi*perm43(CurB0WVar).*TimePoints_Med_ms3/1e3)));
clear ksp_adjB
for s=1:nShots
    tmp=bart(['linopScript -A ' ScriptFN_CompgBo],Sz16CompgB0S,CurData(1,:,s,1:nccy,1,1,1,1,:),...
        CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj(:,:,s),TSBPMed,TSCxPMedOnlyB0SWB0Var(:,:,s,1,1,1,:),...
        CurKernsA(:,:,s,1,1,1,:),CompsP);
    ksp_adjB(:,:,s,1,1,:)=tmp;
end
ksp_adjBs=sum(ksp_adjB,3);

Rec_CompgB0_C_RSS_WB0V_OnlyN=bart(['picsS -w 0.001063 ' RecStr ScriptFN_CompgB0_OnlyN],...
    Sz16CompgB0S,ksp_adjBs,...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0SWB0Var,...
    CurKernsA,CompsP);

Rec_CompgB0_C_RSS_WB0V_OnlyNX=squeeze(sum(Rec_CompgB0_C_RSS_WB0V_OnlyN.*CompsP,6));
[PDBase_TVp, UpdatedB0Map_TVp, UpdatedT2SMap_ms_TVp, s_vals_TVp, Fitted0_TVp, PDBase0_TVp]=...
    FitToModel_MPBD1CSf(Rec_CompgB0_C_RSS_WB0V_OnlyNX,WhichTSToUset,dT_Med_ms,TE0_ms+dT_Med);


Rec_CompgB0_C_RSS_NoB0V_OnlyNX=squeeze(sum(Rec_CompgB0_C_RSS_NoB0V_OnlyN.*CompsP,6));
[PDBase_TVa, UpdatedB0Map_TVa, UpdatedT2SMap_ms_TVa, s_vals_TVa, Fitted0_TVa, PDBase0_TVa]=...
    FitToModel_MPBD1CSf(Rec_CompgB0_C_RSS_NoB0V_OnlyNX,WhichTSToUset,dT_Med_ms,TE0_ms+dT_Med);




%%
QQ=uint8(grmss(Rec_CompgB0_RSS_MXCMX,3:4)/(16e-4)*256);
filename = 'Motion.gif'; % Specify the output file name
imwrite(perm43(QQ),gray(256),filename,'gif','LoopCount',Inf,'DelayTime',.25);
%%
figure;
gmontage(double(QQ(:,:,6)),[0 256]);daspect([1 1 1]);removeTicks;
text(10,10,'1','FontSize',12,'Color',[1 1 1]);

for i=1:size(QQ,3)
    gmontage(double(QQ(:,:,i)),[0 256]);daspect([1 1 1]);removeTicks;
    text(10,10,num2str(i),'FontSize',12,'Color',[1 1 1]);
    AA(i)=getframe;
    QQQ(:,:,:,i)=AA(i).cdata;
end
% end
filename = 'Motion.gif'; % Specify the output file name
imwrite(QQQ(:,:,2,:),gray(256),filename,'gif','LoopCount',Inf,'DelayTime',.25);
%% End Motion test

%% More motion stuff: Focus on one and rotate
SliIs=ROrd(30);
% RepSetX=10:13;
RepSetX=14:20;

% ToMove=Rec_CompgB0_RSS_MXCMY(:,:,SliIs,5:8);
[~,Ord]=sort(ROrd);
ToMove=Rec_CompgB0_RSS_MXCMY(:,:,SliIs,RepSetX);
ForRef=RefSV1x(:,:,Ord(SliIs));

RunMcflirtAndFlirt; % gets Moved

fgmontagex(repmat(ForRef,[1 1 numel(RepSetX)]));title('Ref');
fgmontagex(ToMove);title('ToMove');
fgmontagex(Moved);title('Moved');

WPos=get(gcf,'Position');
%% Prepare spaces for translating into comfortable coordinates
RefSpace=imref2d(Sz);

tTranslationToCenterAtOriginb=gTransMat2D(-mean(RefSpace.XWorldLimits),-mean(RefSpace.YWorldLimits));
tTranslationBackToOriginalCenterb=gTransMat2D(mean(RefSpace.XWorldLimits),mean(RefSpace.YWorldLimits));

tTranslationToCenterAtOriginc = gTransMat2D(1,Sz(1));
tTranslationBackToOriginalCenterc = gTransMat2D(-1,-Sz(1));

tTranslationToCenterAtOriginbc=tTranslationToCenterAtOriginb*tTranslationToCenterAtOriginc;
tTranslationBackToOriginalCenterbc=tTranslationBackToOriginalCenterb*tTranslationBackToOriginalCenterc;
%%
TMat=readFSLMatFile([MovFN '_mcf_mean_to_' RefFN '.mat']);
TMats=readFSLMatFolder([MovFN '_mcf.mat']);

TMatsx=MultMatTensor(TMat,TMats);

% Tm1=ApplyTransformAsFLIRT(squeeze(ToMove),TMatsx);

% tmpMatx=tTranslationBackToOriginalCenterbc*tmpMatM*tTranslationToCenterAtOriginbc;
% Tmc=ApplyTransformFromCenter(ForRef,tmpMatx);

TMatsx3=gTMatFLIRTtoM(TMatsx);

TMatsx3c=MultMatTensor(tTranslationBackToOriginalCenterbc,MultTensorMat1(TMatsx3,tTranslationToCenterAtOriginbc));

MovedMC=ApplyTransformFromCenter(squeeze(ToMove),TMatsx3c);

fgmontagex(MovedMC);title('MovedMC');

Rot_deg=-asind(squeeze(TMatsx3(1,2,:)));
Trans=res_mm*squeeze(TMatsx3(3,1:2,:));
figure;plot(Rot_deg,'o-');hold on;plot(Trans.','o-');xlabel(mainP,'Interpreter','None');legend({'Rotation','x','y'});
%%
dx_vox=dx*Sz(1);
dy_vox=dy*Sz(2);

% TBase=gTransMat2D(-dy_vox,-dx_vox);
TBase=gTransMat2D(dy_vox,dx_vox);

TMatsx3cc=MultMatTensor(TBase,TMatsx3c);

RotsD=-asind(squeeze(TMatsx3cc(1,2,:)));
DxTtl_vox=squeeze(TMatsx3cc(3,2,:));
DyTtl_vox=squeeze(TMatsx3cc(3,1,:));
% DxTtl_vox=squeeze(TMatsx3c(3,2,:));
% DyTtl_vox=squeeze(TMatsx3c(3,1,:));
%%
% RepSetX=11;

CurSens=perm93(RefSensCCS(:,:,SliIs,1:nccToUse));
CurSens=imresize(CurSens,RSz);

TSCxPMedOnlyB0S=perm73(exp(-1i*2*pi*perm93(RefB0xSr(:,:,SliIs)).*TimePoints_Med_ms3/1e3));

i=1;
SlicesToRead=SliIs(i);
SliI=SlicesToRead;

sccmtx=sccmtxS(:,:,SlicesToRead);

ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,:,:,:,RepSetX(i,:),:,:,:,:,:,:,:,:);
ADataIsL=permute(ADataIsL,[1 2 9 11 3:8 10 12:20]);
ADataIsL=CombineDims(ADataIsL,[4 1]);

ADataIsLCC=single(zeros([size(ADataIsL,1) nccToUse size(ADataIsL,3)]));
for c=1:nccToUse
    ADataIsLCC(:,c,:)=sum(ADataIsL.*permute(sccmtx(:,c),[3 1 4 5 6 7 8 9 2]),2);
end
DataC=permute(ADataIsLCC,[1 3 2]);
disp(['Prepared for SliI ' num2str(SliI) ' RepSetX ' num2str(RepSetX)]);
%%
disp(['SliI ' num2str(SliI) ' RepSetX ' num2str(RepSetX)]);

nShots=size(RepSetX,2);
i=1;
% RotsD=[0 7.5 5.5 2];
% RotsD=[0 0 0 0];
% RotsD=7.5;

% dx,dy are in ration to FOV
% dx=RotatedLocs(2,1)/FOVx;
% AddDxs=0/Sz(1); % positive is downward
% AddDys=7/Sz(1); % positive is rightward
% DxTtl_vox=dx_vox;
% DyTtl_vox=dy_vox;

% DxTtl_voxm=[dx_vox  0   dx_vox+7    0           dx_vox  0       dx_vox+7    0           dx_vox+7    dx_vox+7    ];
% DyTtl_voxm=[dy_vox  0   0           dy_vox+5    dy_vox  0       0           dy_vox+5    dy_vox+5    dy_vox+5    ];
% RotsDm    =[0       0   0           0           7.5     7.5     7.5         7.5         7.5         -15         ];
% for m=1:numel(DxTtl_voxm)
%     DxTtl_vox=DxTtl_voxm(m);
%     DyTtl_vox=DyTtl_voxm(m);
%     RotsD=RotsDm(m);
%     disp(['---------- ' num2str(m) ' -------------']);
    
clear CurKernsA Curmodx CurSTraj
for k=1:nShots
    disp(k);
    CurIdx=RepsmIxS(SliIs(i),RepSetX(i,k));
    CurTraj=TrajM(CurIdx,:);
    RotD=RotsD(k);

% CurTraj=CurTraj.*exp(1i*2*pi*-10/360);
% RotD=2; % positive is anti-clockwise
% 11: 7.5 ok
% 12: 5.5 ok
% 13: 2 ok

    CurTraj=CurTraj.*exp(1i*2*pi*RotD/360);
% CurTraj=CurTraj.*exp(1i*2*pi*3/360);
% BARTTraj=cat(3,real(TrajM),imag(TrajM),imag(TrajM)*0);
% BARTTrajP=permute(BARTTraj,[3 2 1]);
% kx=BARTTrajP(1,:,:)*2*pi;
% ky=BARTTrajP(2,:,:)*2*pi;
    
%     Curmodx=exp(-1i*2*pi.*(dx*real(CurTraj)+dy*imag(CurTraj)));
%     CurDx=dx+AddDxs(k);
%     CurDy=dy+AddDys(k);
    CurDx=DxTtl_vox(k)/Sz(1);
    CurDy=DyTtl_vox(k)/Sz(1);
%     Curmodx(1,:,k)=exp(-1i*2*pi.*(dx*real(CurTraj)+dy*imag(CurTraj)));
    Curmodx(1,:,k)=exp(-1i*2*pi.*(CurDx*real(CurTraj)+CurDy*imag(CurTraj)));
    
% Curmodx=Curmodx(:,:,CurIdx);

% clear CurKerns CurSTraj
% CurKerns=single(zeros([RSz*2 1 1 1 1 nTSMed 1 nSliIs]));
% CurSTraj=single(zeros([3 size(STraj3MMed,2) size(RepSetX,2) 1 1 1 1 1 nSliIs]));

% for i=1:nSliIs
%     CurSTraj(:,:,:,:,:,:,:,:,i)=STraj3MMed(:,:,RepsmIxS(SliIs(i),RepSetX(i,:)));
%     CurKerns(:,:,:,:,:,:,:,1,i)=sum(KernsPMMedy(:,:, RepsmIxS(SliIs(i),RepSetX(i,:)) ,:,:,:,:),3);
% end


%%%%%%%%%
% CurSTraj=STraj3MMed(:,:,CurIdx);

%     CurSTraj=CTo3Rows(CurTraj(1,TrajPartMed));
    CurSTraj(:,:,k)=CTo3Rows(CurTraj(1,TrajPartMed));
    
% CurKerns=KernsPMMedy(:,:,CurIdx ,:,:,:,:);

%     CurKerns=getKernsFromTrajM(CurTraj(1,TrajPartMed),RSz,TSBMed);
    CurKernsA(:,:,k,1,1,1,:)=getKernsFromTrajM(CurTraj(1,TrajPartMed),RSz,TSBMed);
end
%
%     clear TrajM
%     for i=1:nRepsToUse
%         CurJitterIdx=mod(i-1,5)+1;
%         Traj=interp1((0:(size(kKA,1)-1))*GradDwellTime_us,kKA(:,CurJitterIdx),AcqTimePoints_us);
%         TrajM(i,:)=Traj.*exp(-1i*2*pi*PhiRotPerRep/360*(i-1));
%     end
% 
%     clear STraj3MMed
% for CurRep=1:(nRepsHdr-1)
%     STraj=TrajM(CurRep,TrajPartMed);
%     STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
% end


% KernsPMMedBase=getKernsFromTrajM(TrajM(1:36,TrajPartMed),RSz,TSBMed);
% KernsPMMedBase=single(KernsPMMedBase);
% KernsPMMedy=KernsPMMedBase(:,:,mod((1:(nRepsHdr-1))-1,36)+1,:,:,:,:,:,:,:,:);
%%%%%%%%%%%%%%%%


CurData=single(zeros([1,nTrajToUse,size(RepSetX,2),nccToUse,1,1,1,1,nSliIs]));

% for i=1:nSliIs
%     DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(1,:,CurIdx);
    DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*Curmodx;
    DataCCP=permute(DataPC,[1:3 8 4:7]);
    CurData(1,:,:,:,1,1,1,1,i)=DataCCP;
% end
CurData=single(CurData);
disp('Calculated CurKernsA Curmodx CurSTraj');
%% First per-shot
for k=1:nShots
    Rec_CompgB0_PerShotC{k}=bart(['picsS ' RecStrMotion ScriptFN_CompgBo],...
        Sz16CompgB0S,CurData(1,:,k,1:nccy,1,1,1,1,:),...
        CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj(:,:,k),TSBPMed,TSCxPMedOnlyB0S,...
        CurKernsA(:,:,k,1,1,1,:),CompsP);
end
Rec_CompgB0_PerShotCM=cat(4,Rec_CompgB0_PerShotC{:});
Rec_CompgB0_PerShotCMX=squeeze(sum(Rec_CompgB0_PerShotCM.*CompsP,6));
Rec_CompgB0_PerShotCMX=Rec_CompgB0_PerShotCMX.*perm43(squeeze(TSCxPMedOnlyB0S));
Rec_CompgB0_PerShotCMY=grmss(Rec_CompgB0_PerShotCMX,4);
%% Together
RecStr=[' -m -S -d 4 -R W:3:0:0.01 '];

CurKerns=sum(CurKernsA,3);

Rec_CompgB0_C_RSS=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

Rec_CompgB0_RSS_M=cat(4,Rec_CompgB0_C_RSS);
Rec_CompgB0_RSS_MX=perm43(squeeze(sum(Rec_CompgB0_RSS_M.*gpermute(CompsP,[4 1]),6)));
Rec_CompgB0_RSS_MX=Rec_CompgB0_RSS_MX.*perm43(squeeze(TSCxPMedOnlyB0S));
Rec_CompgB0_RSS_MX=single(Rec_CompgB0_RSS_MX);

Rec_CompgB0_RSS_MXC{m}=Rec_CompgB0_RSS_MX;

% end

fgmontagex(grmss(Rec_CompgB0_RSS_MX,4));
% title(num2str([RepSetX RotD],'%.1f '));
% fgmontagex(Rec_CompgB0_RSS_MXCMX(:,:,SliIs,:,4))
%%
WhichTSToUsea=round(nTSMed.*0.15):round(nTSMed.*0.85);
[PDBase_TV, UpdatedB0Map_TV, UpdatedT2SMap_ms_TV, s_vals_TV, Fitted0_TV, PDBase0_TV]=...
        FitToModel_MPBD1CSf(squeeze(Rec_CompgB0_RSS_MX),WhichTSToUsea,dT_Med_ms,TE0_ms+dT_Med);
Fitted0_TV=single(Fitted0_TV);

%%
fgmontagex(grmss(Rec_CompgB0_RSS_MX0,4));title('Rec_CompgB0_RSS_MX0');
fgmontagex(grmss(Rec_CompgB0_RSS_MXr,4));title('Rec_CompgB0_RSS_MXr');
fgmontagex(grmss(Rec_CompgB0_RSS_MXrx,4));title('Rec_CompgB0_RSS_MXrx');
fgmontagex(grmss(Rec_CompgB0_RSS_MXry,4));title('Rec_CompgB0_RSS_MXry');
%%
SmallVarT=150e6;
getSmallVars;
save([mainP filesep 'SmallVars3.mat'],SmallVars{:});

%% End motion test specific


%% Multi TV test
Lambdas=10.^(-2:0.5:1);

for j=1:numel(Lambdas)
    disp(j);
    disp(datestr(now));
    RecStr=[' -m -S -d 5 -R T:256:0:' num2str(Lambdas(j)) ' ' ];
    Rec_CompgB0_CTV{j}=bart(['picsS ' RecStr ScriptFN_CompgBo],...
        Sz16CompgB0S,CurData,...
        CurSens,CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
        CurKerns,CompsP);
end
%%
clear Rec_CompgB0_C_TVCMX
Rec_CompgB0_C_TVCM=cat(5,Rec_CompgB0_CTV{:});
Rec_CompgB0_C_TVCMX=sum(Rec_CompgB0_C_TVCM.*CompsP,6);
%%
WhichTSToUseTV=2:(nTS_THLR-1);
[PDBase_TV, UpdatedB0Map_TV, UpdatedT2SMap_ms_TV, s_vals_TV, Fitted0_TV, PDBase0_TV]=...
        FitToModel_MPBD1CSf(squeeze(Rec_CompgB0_C_TVCMX(:,:,1,1,6,1,:,1,:)),WhichTSToUseTV,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
Fitted0_TV=single(Fitted0_TV);
[PDBase_TV7, UpdatedB0Map_TV7, UpdatedT2SMap_ms_TV7, s_vals_TV7, Fitted0_TV7, PDBase0_TV7]=...
        FitToModel_MPBD1CSf(squeeze(Rec_CompgB0_C_TVCMX(:,:,1,1,7,1,:,1,:)),WhichTSToUseTV,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
Fitted0_TV7=single(Fitted0_TV7);
save([mainP filesep 'TVRec6.mat'],'UpdatedT2SMap_ms_TV','Fitted0_TV','s_vals_TV');
save([mainP filesep 'TVRec7.mat'],'UpdatedT2SMap_ms_TV7','Fitted0_TV7','s_vals_TV7');
%% Multi lambda,rho test: doesn't change much. 1e-1 and lower good for both
Rhos=10.^(-3:1);
Lambdas=10.^(-3:1);

for i=1:numel(Rhos)
    for j=1:numel(Lambdas)
        disp([i j]);
        disp(datestr(now));
%         RecStr=[' -m -S -d 5 -R T:256:0:' num2str(Lambdas(j)) ' ' ];
        RecStrS=[' -m -S -d 5 -u ' num2str(Rhos(i)) ' -R W:259:0:' num2str(Lambdas(j)) ' '];
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




































%%


%%
SmallVarT=50e6;
getSmallVars;
save([mainP filesep 'SmallVars.mat'],SmallVars{:});
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
for SliI=1:nSlices
    sccmtx=sccmtxS(:,:,SliI);
    RSensCC=permute(sum(RefSensx(:,:,:,Ord(SliI)).*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    RSensCC=permute43(RSensCC);
    RSensCCS(:,:,:,SliI)=RSensCC;
end
SensCCS=permute(RSensCCS,[1 2 5 3 4]);
%%
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
Rhos=10.^(-3:1);RecStr=[' -m -S -d 4 -R T:256:0:1 -R T:3:0:.001 '];
Rec_CompgB0_C_RSS=bart(['picsS ' RecStr ScriptFN_CompgBo],...
    Sz16CompgB0S,CurData(1,:,:,1:nccy,1,1,1,1,:),...
    CurSens(:,:,1,1:nccy,1,1,1,1,:),CurSTraj,TSBPMed,TSCxPMedOnlyB0S,...
    CurKerns,CompsP);

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