FN='/home/a/SpiralData/27Apr18/meas_MID80_gSpiral_BOLD_1Sli_1Rep_1Shot_VD10_6ADCs_25_125_FID15656.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID78_gSpiral_BOLD_1Sli_1Rep_1Shot_VD13_6ADCs_25_125_FID15654.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID85_gSpiral_BOLD_1Sli_1Rep_2Shot_VD10_4ADCs_25_125_FID15660.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID87_gSpiral_BOLD_1Sli_4Rep_2Shot_VD10_4ADCs_25_125_FID15662.dat';

FN='/media/a/DATA1/13May18/Phantom/meas_MID364_gBP_VD11_U19_FID17753.dat';
FN='/media/a/DATA1/13May18/Phantom/meas_MID366_gBP_VD11_U10_FID17755.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID407_gBP_VD11_U10_FID17796.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID405_gBP_VD11_U19_FID17794.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID417_gBP_VD15_U20_FID17806.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID426_gBP_VD11_U19_A35S155_FID17815.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID409_gBP_VD11_U19_7ADCs_FID17798.dat';

% FN='/media/a/DATA/14May18/Ben/meas_MID107_gBP_VD11_U19_FID17942.dat';

FN='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944.dat';

ScanP='/media/a/DATA/13May18/Me/';
BaseFN='meas_MID409_gBP_VD11_U19_7ADCs_FID17798';
RefFldMapP='/media/a/DATA/13May18/Me/meas_MID399_BP_fieldmap_4echos_FID17788/';

% ScanP='/media/a/DATA/14May18/Ben/';
% BaseFN='meas_MID111_gBP_VD11_U19_G35S155_FID17946';
% BaseFN='meas_MID109_gBP_VD11_U19_4min_FID17944';
% RefFldMapP='/media/a/DATA/14May18/Ben/meas_MID123_BP_fieldmap_5echosX_FID17958/';
% 
% BaseP='/media/a/DATA/180628_AK/';
% ScanP='/media/a/DATA/180628_AK/';
% RefFldMapP='/media/a/DATA/180628_AK/meas_MID265_BP_fieldmap_5echosX_FID22460/';
% BaseFN='meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439';

ScanP='/media/a/DATA/11Jul18/RL/';
BaseFN='meas_MID149_gBP_VD11_U19_G35S155_FID23846';
RefFldMapP='/media/a/DATA/11Jul18/RL/meas_MID141_BP_fieldmap_5echosX_FID23838/';

% ScanP='/media/a/DATA/PhantomCAIPI/';
% RefFldMapP='/media/a/DATA/PhantomCAIPI/meas_MID330_BP_fieldmap_5echosX_FID24027/';
% BaseFN='meas_MID342_gBP_Spi_12Sli_MB2_FID24036';
% BaseFN='meas_MID344_gBP_Spi_12Sli_MB2_Cs36_FID24038';

FN=[ScanP BaseFN '.dat'];
%%
mkdir([ScanP BaseFN]);
%% Read raw
AData = mapVBVD(FN);
% ADataI=AData.image();
ADataIx=AData.image(:,:,:,:,:,3,:,:,:,:,:,:,:);
% ADataIs=squeeze(ADataI);
% ADataIs=ADataIs(ADataIs(:,:,:,3,:)); % remove weird dimension with zeros
% ADataIs=squeeze(ADataIs);
ADataIsL=squeeze(ADataIx);

for i=1:numel(AData.hdr.Phoenix.sSliceArray.asSlice)
    SLoc(i)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dTra;
end

nSlices=numel(AData.hdr.Phoenix.sSliceArray.asSlice);
% Ord=[1:2:nSlices 2:2:nSlices];
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);


FOVx=AData.hdr.Meas.ReadFOV;
dFOV=FOVx/1000;

paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{20};
spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{6};
paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{5};
paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
MB=AData.hdr.MeasYaps.sWiPMemBlock.alFree{9};
% MB
if(MB>1)
    spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{14};
    paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
    paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};
    paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
end
%%
save([ScanP BaseFN filesep 'Data.mat'],'ADataIsL','AData');
disp(['Saved ' [ScanP BaseFN filesep 'Data.mat']]);
%%
% save('StatusForTrajDesign1Band_2.mat');
%% given everything but VD and Acc
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

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
disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(SlewBuf);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'Traj'],[]) 
close(gcf);
save([ScanP BaseFN filesep 'Traj.mat'],'BARTTrajMS');
%%
Trajm2=BARTTrajMS(1:2,1:end-2);
Sz128=[128 128];

[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
% save('ForTFNUFT.mat','SN','Kd','P','A','NUbyFS3');
save([ScanP BaseFN filesep 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');
disp('Saved TrajForNUFT');
%%
% Traj=BARTTrajMS(1:2,1:end-2);
% save([ScanP BaseFN filesep 'Traj.mat'],'Traj');
%%
if(MB>1)
    SliOffset=0;
else
    SliOffset=nSlices/2;
end
SensB=load([RefFldMapP 'Sens.mat']);
SensB=SensB.SensB;
SnsSzB=gsize(SensB,1:2);

% B0_HzU=load([RefFldMapP 'B0_HzU.mat']);
% B0_HzU=B0_HzU.B0_HzU;
% B0Q=imresizeBySlices(B0_HzU,SnsSzB);

FirstEcho=load([RefFldMapP 'FirstEcho.mat']);
FirstEcho=FirstEcho.FirstEcho;

FirstEcho=gflip(FirstEcho(:,:,:,SliOffset+(1:nSlices)),1:2);
Mg=grmss(FirstEcho,3);

B0S=load([RefFldMapP 'B0S.mat']);
B0S=B0S.B0S;
disp('ok');

% load('CurStatus_Ben4MinASL_x.mat','SensB','B0Q','Mg'); ;
%%
% SensX=SensB(:,:,:,6+(1:12),:);
SensX=permute(SensB(:,:,:,SliOffset+(1:nSlices),:),[1 2 3 5 4]);
SensX=gflip(SensX,1:2);

% B0Q2=B0Q(:,:,6+(1:12));
B0Q2=B0S(:,:,SliOffset+(1:nSlices));
B0Q2=gflip(B0Q2,1:2);
disp('ok');
%%
ADataIsP=double(permute(ADataIsL,[1 5 2 3 4]));
nReps=size( ADataIsP,5);
nChannels=size( ADataIsP,3);
%
try
    dx=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
catch
    disp('No x shift!');
    dx=0;
end
% dy=-15;
dy=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
% dy=-45;
read_raw=1;

kx=BARTTrajMS(1,:)*2*pi;
ky=BARTTrajMS(2,:)*2*pi;
modx=double(exp(1i*(dx*kx+dy*ky))');


ADataIsPx=reshape(ADataIsP,[paramLongROSamples,nChannels nSlices nReps]);
ADataIsPy=ADataIsPx(1:size(BARTTrajMS,2),:,:,:,:);

ADataIsPy=RepDotMult(ADataIsPy,modx);

clear ADataIsPx ADataIsP ADataIx ADataIsL
disp('ok a')
%%
SliI=3;
%%
nukData=ADataIsPy(:,:,SliI,1).';

nukData=nukData(:,3:end);
nTraj=size(nukData,2);
TimeInMs2=(0:nTraj-1)*2.5/1e3;
T2SCompStr='';
T2SEstMs=20;
T2SEstDecay=exp(-TimeInMs2/T2SEstMs);
% T2SEstComp=exp(TimeInMs2/T2SEstMs);

Sz2=gsize(SensX,1:2);

BARTTrajAct=BARTTrajMS(:,1:size(nukData,2));
%%
nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
% nukData=nukData.*T2SEstComp;
nukDataP=permute(nukData,[3 2 4 1]);

SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

% %% MB
% nukDataP=nukDataP(:,:,1);
% SensP=ones(Sz2);
%%

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));

RecIfTVsMM=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP);
RecIfWsMM=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP);
% ScrfTV=@(x) gScoreImageSimilarity(RecIfTV(x),ARefImg,RefDx,RefDy);

Lambda=1e-6;
% Rec=RecIfTVs(Lambda);
Rec=RecIfWs(Lambda);

RecMM=RecIfWsMM(Lambda);

fgmontage(cat(3,Rec,RecMM(:,:,1)));
xlabel(['BART No n0: Left - Single map, ADMM. Right - with 2 maps, Wavelet Lambda= ' num2str(Lambda)]);

% fgmontage(RecMM)
%
title(['BART No B0 W=' num2str(Lambda)]);
YLbl=['Sli' num2str(SliI,'%02d')];
ylabel(YLbl);
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_BARTRecon_NoB0_W' num2str(Lambda) T2SCompStr],[]) 
close(gcf);
save([ScanP BaseFN filesep YLbl '_BARTRecon_NoB0_W' num2str(Lambda) T2SCompStr '.mat'],'Rec');
disp(['Saved ' [ScanP BaseFN filesep YLbl '_BARTRecon_NoB0_W' num2str(Lambda) T2SCompStr '.mat']]);
%% All timepoint
if(nReps>20)
    RecT=zeros([size(Rec) nSlices nReps]);
    RecIfWsD=@(x,D,SensP) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, D, SensP(:,:,:,:,1));
    for SliI=1:nSlices
        for r=1:nReps
            nukData=ADataIsPy(:,:,SliI,r).';
            nukData=nukData(:,3:end);
            % nukData=nukData.*T2SEstComp;
            nukDataP=permute(nukData,[3 2 4 1]);
            
            SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);
            
            RecBART(:,:,SliI,r)=RecIfWsD(Lambda,nukDataP,SensP);
            disp(['ok ' num2str(SliI) ' r' num2str(r)]);
        end
    end
    save([ScanP BaseFN filesep 'BARTRecon_NoB0_W' num2str(Lambda) T2SCompStr '.mat'],'RecBART');
end
% MCE=abs(mean(RecT(:,:,SliI,2:2:end),4));
% MCO=abs(mean(RecT(:,:,SliI,1:2:end),4));
% D=MCE-MCO;
% fgmontage(rot90(D),[0 7e-3])
% 
% fgmontage(abs(std(RecT(:,:,SliI,1:2:end),[],4)),[0 5e-5])
%%
figure;
subplot(2,2,1);
gmontage(Rec);title('BART recon');
ylabel(YLbl);

subplot(2,2,2);
gmontage(Mg(:,:,SliI));title('FieldMap 1st echo mag');
subplot(2,2,3);
gmontage(SensX(:,:,:,SliI));title('Fieldmap sens');
subplot(2,2,4);
gmontage(B0Q2(:,:,SliI),[-300 300]);title('Fieldmap B0 Hz');colorbar
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_Params'],[]) 
close(gcf);
disp(['printed fig to file' [ScanP BaseFN filesep YLbl '_Params']]);

%%
Lambda=1e-5;
for SliI=1:nSlices
    disp(['--- Sli #' num2str(SliI) ' --- ' datestr(now)]);
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);
%     nukData=nukData.*T2SEstComp;
    nukDataP=permute(nukData,[3 2 4 1]);
    SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

    RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));
    RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));
    RecIfWsMM=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP);

    RecS(:,:,SliI)=RecIfWs(Lambda);
    tmp=RecIfWsMM(Lambda);
    RecSMM(:,:,SliI)=tmp(:,:,1);
end
disp('Finished BART No b0 all slices');
%%
fgmontage(RecS)
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'BARTRecon_AllS_NoB0_W' num2str(Lambda) T2SCompStr],[]) 
close(gcf);
save([ScanP BaseFN filesep 'BARTRecon_AllS_NoB0_W' num2str(Lambda) T2SCompStr '.mat'],'RecS');

fgmontage(RecSMM)
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'BARTRecon2Maps_AllS_NoB0_W' num2str(Lambda) T2SCompStr],[]) 
close(gcf);
save([ScanP BaseFN filesep 'BARTRecon2Maps_AllS_NoB0_W' num2str(Lambda) T2SCompStr '.mat'],'RecSMM');
%%
SliI=3;
%% All B0 effects across time
% Mgc=imresizeBySlices(gflip(Mg(:,:,SliI+6),1:2),Sz2);
% Mgc=imresizeBySlices(Mg(:,:,SliI+6),Sz2);
Mgc=imresizeBySlices(Mg(:,:,SliI),Sz2);
Mskc=Mgc>7e-5;

% B0M2=-B0Q2(:,:,SliI);
% B0M2(~Mskc)=0;
% 
% SymMskC=abs(B0M2-gflip(B0M2,1))>230;
% SymMskC(1:30,:,:)=true;
% 
% B0M2(SymMskC & B0M2>150)=-20;
% fgmontage(B0M2)
% 
% cx=caxis
% fgmontage(X2,cx)
% 
% B0M2=X2;
B0M2=B0Q2(:,:,SliI);

% B0M2=-B0M2;
% B0M2(B0M2>150)=-20;
% fgmontage(B0M2)
%%
% B0M2=B0RealEx;

[U_TimeInMs2, IA_TimeInMs2, IB_TimeInMs2]=unique(TimeInMs2);
nU_TimeInMs2=numel(U_TimeInMs2);

AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
MgcN=Mgc./grmss(Mgc);
WE=Col(MgcN);
% WE=Col(MgcN)*0+1;

WeightedE=WE.*E;
%% Fessler time segmentation
% nTS=7;
clear ErrTS
TS_Thresh=1e-5;
for nTS=8:15
    disp(nTS)
    FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    TSC2=reshape(TSC,prod(Sz2),nTS);
    WTSC2=WE.*TSC2;
%     tic
%     TSB=(E.')/(TSC2.');% W in both sides
    TSB=(WeightedE.')/(WTSC2.');% W in both sides
    
%     ErrTS(nTS)=grmss(E-TSC2*(TSB.')); %include W
    ErrTS(nTS)=grmss(WeightedE-WTSC2*(TSB.')); %include W
    
    disp([datestr(now) 'nTS ' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
    if(ErrTS(nTS)<TS_Thresh)
        disp(['Stopped at #TS=' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
        break;
    end
end
figure(87234);clf;plot(log10(ErrTS),'-*')
%% GPU TS
% Sens=imresizeBySlices( squeeze(SensP2),Sz2);
osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

% TSC_MB=ones([NTrg1 1 nBands]);

% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),NTrg1,osf,wg,sw,ones(1,size(BARTTrajMS,2),nBands),TSC_MB,SensX(:,:,:,SliI));

NTrg1=Sz2;
nBands=1;
nShots=paramLongInterleaves;

TSBF=TSB(IB_TimeInMs2,:).';
Sens=squeeze(SensX(:,:,:,1,SliI));

GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajAct,NTrg1,osf,wg,sw,TSBF,TSC,Sens);
GOP_MC=GOP_MCSMBMS;
nTrajP2=nU_TimeInMs2;
disp('ok');
%%
SliIP=[ScanP BaseFN filesep YLbl filesep];
mkdir(SliIP);
save([SliIP 'TSBF_TSC_Sens.mat'],'TSBF','TSC','Sens');
disp(['Saved ' SliIP 'TSBF_TSC_Sens.mat']);
%%
% GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);
x = randn(Sz2) + 1j*randn(Sz2);
y = randn([size(Sens,3) nTrajP2]) + 1j*randn([size(Sens,3) nTrajP2]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
%%
% T2SCompStr='_T2S20Comp';
%%
% if(nShots==1)
%     TVOP=TVOP_MSlice;
% else
%     TVOP=TVOP_MTC_W([1 1 0 1e1]);
% end
% TVOP=Wavelet;
% XFM=1;
% XFM = Wavelet('Daubechies',4,4);	% Wavelet
% XFM = Wavelet('Daubechies_TI',4,4);	% Wavelet

nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
% nukData=nukData.*T2SEstComp;
    
DataP=nukData;

% QQ=load('/media/a/H1/TFTSNUFTOut.mat');
% WithTSBR=double(QQ.Out)/100;
% DataP=WithTSBR.';

AOdd = GOP_MCSMBMS;
% AOdd = GOP_MC;

% TVW=0.1;
TVW=1e-7;


% filterType:   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies','Symmlet', 'Vaidyanathan','Battle'
% Use suffix _TI for translation invariance, for example, 'Daubechies_TI'
% filterSize: related to the support and vanishing moments of the particular wavelet (See MakeONFilter in wavelab)
% wavScale: 	scallest scale of wavelet decomposition

XFMStr='Daubechies';

filterSize=4;
wavScale=4;

if(strcmp(XFMStr,'None'))
    XFM=1;
    xfmWeight=0;
else
    XFM = Wavelet(XFMStr,filterSize,wavScale);
    xfmWeight = 1e-7;	% Weight for Transform L1 penalty
end


param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);

% XFM=SubSpaceOp(V,SdT,1);
% xfmWeight=1e-1;
% % param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);
% % param=ExtendStruct(struct('pNorm',2,'TVWeight',0,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',1,'xfmWeight',xfmWeight),init);
% param=ExtendStruct(struct('pNorm',2,'TVWeight',0.01,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);


param.data =     DataP;

nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
if(nShots>1)
    res=repmat(resA,[1 1 1 nShots]);
    res=res+randn(size(res))*max(abs(resA(:)))/20;
end

FigH=4000;
figure(FigH);close(FigH);

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
param.Verbose=false;
clear ObjConv Score
for n=1:nfnlCgIters
% for n=1:2
%     disp([Slis WhichRep n]);
    [res, CurObj] = fnlCg(res,param);
%     res=repmat(mean(res,4),[1 1 1 2]);
    ObjConv(n)=CurObj;
%     (ObjConv(end)-ObjConv(end-1))*2/(ObjConv(end)+ObjConv(end-1))
    im_res = param.XFM'*res;
%     figure(FigH), gmontage(abs(gflip(im_res,1))), drawnow;% title(qq)\
%     Score(n)=grmss(CurMBMSIs-im_res);
    if(param.ShowFig)
        figure(FigH); subplot(1,3,1);
        gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,3,2);
        gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
        subplot(1,3,3);
        plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
        
%         subplot(2,3,4);
%         gmontage(abs(CurMBMSIs-im_res),[-100 100]);
        
%         if(nShots>1)
%             subplot(2,3,5);
%             gmontage(RepDotMult(Msks(:,:,1:nBands),abs(cat(4, diff(CurMBMSIs,[],4),diff(im_res,[],4)))),[-100 100]);
%         end
        
%         subplot(2,3,6);
%         plot(Score);setYaxis([0 Score(n)*3]);if(n>1), setXaxis([1 n]);end
    end
%     t=toc;
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

fgmontage(im_res,[0 7e-3]);
%%
XFMStrFull=['[' XFMStr ',' num2str(filterSize) ',' num2str(wavScale) ',' num2str(xfmWeight) ']'];
%     XFMStr ' Size=' num2str(filterSize) ' Scale=' num2str(wavScale) ' W=' num2str(xfmWeight)
XLbl=['L' num2str(param.pNorm) ',TVW=' num2str(param.TVWeight) ',' XFMStrFull ];
xlabel(XLbl)
YLbl=['Sli' num2str(SliI,'%02d')];
ylabel(YLbl);

gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_' XLbl T2SCompStr],[]) 
close(gcf);
save([ScanP BaseFN filesep YLbl '_' XLbl T2SCompStr '.mat'],'im_res');
disp(['Saved ' ScanP BaseFN filesep YLbl '_' XLbl T2SCompStr '.mat']);
%% ESPIRIT
% GOP_TS=ggpuNUFT_TS(GOP2,TSBF);
% 
% 
% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajAct,NTrg1,osf,wg,sw,TSBF,TSC,Sens);
% 
% GOP=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTraj(1:2,:,s),Sz)'/(2*pi),ones(nTrajP,1),osf,wg,sw,Sz,Sens_TS(:,:,:,b));
% GOP2=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,TSC);
% 
% QQ=GOP_MC'*DataP;
% 
% QQ=GOP2'*DataP(1,:).';
% 
% QA=GOP_TS*rand(128);
% QB=GOP_TS'*QA;
% 
% QB2=GOP_TS'*DataP;
%%
GOP2=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,double(TSC));
GOP_TS=ggpuNUFT_TS(GOP2,double(TSBF));
GOP_TS_ESP=ggpuNUFT_TS_MCnosum(GOP_TS,[1 1 32]);

CurSens=SensX(:,:,:,:,SliI);
weights=ones(gsize(CurSens,[1 2 4]));
ESP = ESPIRiT(CurSens,weights);

DataP2=double(permute(DataP,[3 2 1]));
%%
XOP = Wavelet('Daubechies_TI',4,6);

splitWeight=0.5;
lam=1e-6;
% nIterSplit=15;
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, im_res*0, GOP_TS_ESP, ESP, 5,XFM,lam,splitWeight,nIterSplit);
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, weights*0, GOP_TS_ESP, ESP, 5,XOP,lam,splitWeight,nIterSplit);
% This works but we skiup and go directly to CC
[resL1ESPIRiT] = cgL1ESPIRiT(DataP2, weights*0, GOP_TS_ESP, ESP, 5,XOP,lam,splitWeight,2e-3);
disp('ok cgL1ESPIRiT');
%%
% save('CurStatusForL1ESPIRIT_MM_Ben4Min.mat')
%%
resL1ESPIRiT1=resL1ESPIRiT(:,:,1);
fgmontage(resL1ESPIRiT1,[0 7e-3])
title('resL1ESPIRiT 2 maps, B0');
ylabel(YLbl);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);

gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr ],[]) 
close(gcf);
save([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '.mat'],'resL1ESPIRiT1');
disp(['Saved ' ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '.mat']);
%%
nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
%%
[U,S,sccmtx] = svd(nukData.','econ');
ncc=13;

CurSens=SensX(:,:,:,:,SliI);
weights=ones(gsize(CurSens,[1 2 4]));
ESP = ESPIRiT(CurSens,weights);

SensCC=permute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2 4])),[2 3 1 4]);
weightsCC=ones(gsize(CurSens,[1 2 4]));
ESPCC = ESPIRiT(SensCC,weightsCC);

GOP_TS_ESPCC=ggpuNUFT_TS_MCnosum(GOP_TS,[1 1 ncc]);

% GOP_MC_CC = ggpuNUFT_TS_MCx(BARTTrajAct,Sz2,osf,wg,sw,TSB.',TSC,SensCC(:,:,:,1));


    
DataCC=(sccmtx(:,1:ncc).'*nukData);

DataPCC=double(permute(DataCC,[3 2 1]));

XOP = Wavelet('Daubechies_TI',4,6);

splitWeight=0.5;
lam=1e-8;
nIterSplit=2e-3;
% nIterSplit=15;
% [resL1ESPIRiT] = cgL1ESPIRiT(DataP2, im_res*0, GOP_TS_ESP, ESP, 5,XFM,lam,splitWeight,nIterSplit);
[resL1ESPIRiTCC] = cgL1ESPIRiT(DataPCC, weightsCC*0, GOP_TS_ESPCC, ESPCC, 5,XOP,lam,splitWeight,nIterSplit);

disp('ok cgL1ESPIRiTCC');
%%
resL1ESPIRiTCC1=resL1ESPIRiTCC(:,:,1);
fgmontage(cat(3,resL1ESPIRiT1, resL1ESPIRiTCC1),[0 7e-3])
title(['resL1ESPIRiT 2 maps, B0. Right - with CC -> ' num2str(ncc)]);
ylabel(YLbl);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC' num2str(ncc)],[]) 
close(gcf);
save([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat'],'resL1ESPIRiTCC1');
disp(['Saved ' ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat']);
%% Clean save
% BaseOutLoc='/media/a/DATA/';
SensMsk=single(grmss(SensCC(:,:,:),3)>0.01);
BaseOutLoc=ScanP;
CurOurP=[BaseOutLoc BaseFN filesep];
mkdir(CurOurP);
CurOurPSli=[BaseOutLoc BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
mkdir(CurOurPSli);
save([CurOurPSli 'Sens.mat'],'CurSens');
save([CurOurPSli 'SensCC.mat'],'SensCC','sccmtx','SensMsk');
save([CurOurPSli 'B0TS.mat'],'TSBF','TSB','TSC','osf','wg','sw','Mgc','TimeInMs2');
save([CurOurPSli 'TrajAndRealData.mat'],'BARTTrajAct','DataP2','DataPCC');
%% Now calc on all for NN
AllImWithPhaseComplexSingle=load('/media/a/H1/All32kImWithPhaseComplexSingleX128x128.mat');
AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle.AllImWithPhaseComplexSingle;

nI=size(AllImWithPhaseComplexSingle,1);
nTrain=nI;
AllData=single(zeros(nTrain,nTrajP2,ncc));
AllData(1)=0+1i*eps;
%%
% poolobj = gcp;
% % myFiles = {'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/gUM/@ggpuNUFT_TS_MCnosum/mtimes.m'};
% % addAttachedFiles(poolobj, myFiles);
% poolobj.addAttachedFiles(fullfile(pwd, '@ggpuNUFT_TS_MCnosum'));
% poolobj.addAttachedFiles(fullfile(pwd, '@ggpuNUFT_TS'));
% poolobj.addAttachedFiles('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/gpuNUFFT-master/gpuNUFFT/@gpuNUFFT');
    
StartTime=now;
for i=1:nTrain
    if(mod(i,100)==1)
        DeltaTime=now-StartTime;
        FullTime=DeltaTime*nTrain/i;
        EndEstTime=StartTime+FullTime;
        disp([num2str(i,'%05d') ' Cur ' datestr(now) ' est finish ' datestr(EndEstTime)]);
    end
    CurI=double(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
    CurIWChCC=CurI.*SensCC(:,:,:,1);
    CurIData=GOP_TS_ESPCC*CurIWChCC;
%     CurIData=(GOP_MC_CC*CurI).';
    AllData(i,:,:)=single(CurIData);
end
disp('Finished');
%%
i=4455;
fgmontage(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
[resL1ESPIRiTCC] = cgL1ESPIRiT(double(AllData(i,:,:)), weightsCC*0, GOP_TS_ESPCC, ESPCC, 5,XOP,lam,splitWeight,nIterSplit);
%%
weightsCC1=weightsCC(:,:,1);
ESPCC1 = ESPIRiT(SensCC(:,:,:,1),weightsCC1);

ESPCCx = ESPIRiT(SensCC(:,:,1:nChToUseInNN,1),weightsCC1);
GOP_TS_ESPCCx=ggpuNUFT_TS_MCnosum(GOP_TS,[1 1 nChToUseInNN]);

i=4455;
% fgmontage(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
% [resL1ESPIRiTCC1] = cgL1ESPIRiT(double(AllData(i,:,:)), weightsCC1*0, GOP_TS_ESPCC, ESPCC1, 5,XOP,lam,splitWeight,nIterSplit);
[resL1ESPIRiTCCx] = cgL1ESPIRiT(double(AllData(i,:,1:nChToUseInNN)), weightsCC1*0, GOP_TS_ESPCCx, ESPCCx, 5,XOP,lam,splitWeight,nIterSplit);

%%
% save([CurOurPSli 'AllDataCCPhantom' num2str(ncc) '.mat'],'-v7.3','AllData');
save([CurOurPSli 'AllDataCC' num2str(ncc) '.mat'],'-v7.3','AllData');
disp(['Saved to ' CurOurPSli 'AllDataCC' num2str(ncc) '.mat']);
%%
nTrajAct=nTrajP2;
%%
Msk=grmss(SensCC(:,:,:,1),3)>0.01;
MskP=permute(Msk,[3 1 2]);
%
TFDataP=[CurOurPSli 'TF/'];
mkdir(TFDataP);
% nData=5000;
nData=nI;
ChunkSize=100;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);
%%
clear DataPR LabelsP
StartTime=now;
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:).*MskP;
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    
    DataPR=NaN(CurChunkSize,nTrajAct*ncc*2);

    for i=1:CurChunkSize
        CurIData=squeeze(AllData(CurIs(i),:,:));
        CurIDataV=Row(CurIData); % Changed order in AllData!!!
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        DataPR(i,:)=single(CurIDataVR);
    end
    
    % Data to TFRecords
    Data=single(DataPR);
    Labels=single(LabelsP);
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('~/HomeA/TF/CurChunk.mat','Data','Labels','FNs','-v6')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRec.py ~/HomeA/TF/CurChunk.mat ' TFDataP]);
%     if(mod(i,100)==1)
        DeltaTime=now-StartTime;
        FullTime=DeltaTime*numel(ChunkStartI)/k;
        EndEstTime=StartTime+FullTime;
        disp([num2str(k,'%03d') ' Cur ' datestr(now) ' est finish ' datestr(EndEstTime)]);
%     end
end
%%
CurBartTraj=BARTTrajAct;
% CurBartTraj=BARTTrajMS;
nTraj=size(CurBartTraj,2);
figure;plot(CurBartTraj(1,:),CurBartTraj(2,:),'.');
kMax=ceil(max(max(abs(CurBartTraj),[],2)));
xlabel(nTraj);
Acc=(kMax*2)^2/nTraj;
title(['kMax: ' num2str(kMax) ' Acc (for kMax): ' num2str(Acc)]);
%%

osfForNbrhd=1.3;
osN=ceil(kMax*osfForNbrhd)*2+1;

C=linspaceWithHalfStep(-kMax,kMax,osN);
%% Generate Idx mat of neighbors
nNeighbors=12;
NMap=NaN(osN,osN,nNeighbors);
for i=1:osN
    for j=1:osN
        CurLoc=C([i j]).';
        D=CurBartTraj(1:2,:)-CurLoc;
        R=grmss(D,1);
        [~, Idx]=sort(R);
        NMap(i,j,:)=Idx(1:nNeighbors);
    end
end
%%
nChToUseInNN=8;
NMapC=RepDotMult(ones(size(NMap))*nTrajAct,permute(0:nChToUseInNN-1,[1 4 3 2]));
NMapC=NMapC+repmat(NMap,[1 1 1 nChToUseInNN]);
NMapCX=CombineDims(NMapC,[3 4]);
NMapCR=int32(cat(3,NMapCX,NMapCX+nTrajAct*ncc)-1);

% NMapCRC=cat(4,NMapCX,NMapCX+nTrajAct*nChToUseInNN)-1;

NMapFN=[ScanP BaseFN filesep 'NMapIndCB0X8.mat'];
save(NMapFN,'NMapCR');

% NMapCRB=repmat(permute(NMapCR,[4 1 2 3]),[16 1 1 1]);
% NMapCRB=NMapCRB+ int32(permute(0:15,[2 1 3 4])*nTrajAct*ncc*2);
% NMapFNB='/media/a/H1/NMapIndCB0XB.mat';
% save(NMapFNB,'NMapCRB');

%%
tmp=squeeze(DataPCC);
% RealDataFac=grmss(AllData(5656:34:6767,:,:))/grmss(tmp)
RealDataFac=100;
CurIDataV=Row(tmp)*RealDataFac;
CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
Data=repmat(single(CurIDataVR),[16 1]);
RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
save(RealDataFN,'Data');
disp('Saved real data');
%%
% TBaseP='~/HomeA/TF/srez/';
TBaseP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';

DataH=size(NMap,1);
DataW=size(NMap,2);
DataCh=size(NMap,3)*nChToUseInNN*2;
LabelsH=Sz2(1);
LabelsW=Sz2(2);
LabelsCh=2;

ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
  'dataset',TFDataP,'learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0,'WL2_Lambda',0,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'DiscStartMinute',500,...
  'ShowRealData',1,'CmplxBias',0,...
  'InputMode','RegridTry1',...
  'NetMode','RegridTry1C2_TS',...
  'SessionNameBase','RegridTry1C2_TS',...
  'ImgMode','Cmplx',...
  'nTimeSegments',7,...
  'UseSharedWightesInRelaxedFT',1,...
  'WPhaseOnly',0.001,...
  'NMAP_FN',NMapFN,...
  'RealDataFN',RealDataFN);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
%%
ParamsS=ParamsSDefaults;
ParamsS.WPhaseOnly=0.01;
ParamsS.train_time=120;
Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
system('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');

system('sudo bash -c ''echo "I am $USER, with uid $UID"''')
system('sudo -H -u a bash -c ''echo "I am $USER, with uid $UID"''')

sendTFMail(TBaseP,ParamsS,Txt);
%%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0;
ParamsS.WL2_Lambda=0;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);







%%
V=load('/media/a/H1/All32kImWithPhaseComplexSingleX128x128_V.mat');V=V.V;
Sd=load('/media/a/H1/All32kImWithPhaseComplexSingleX128x128_Sd.mat');Sd=Sd.Sd;
%%
SdT=Sd.';
Scores=reshape(im_res,[1 128*128])*V;
NScores=Scores./SdT;
Scores(1,13000:end)=0;
SoftTScores=max(0,(abs(Scores)-0.0005)).*exp(1i*angle(Scores));
Re=SoftTScores*(V');
% Re=Scores*(V');
im_resX=reshape(Re ,[128 128]);
fgmontage(cat(3,im_res,im_resX))
%%
WhichSlices=nSlices:-1:1
%% Multi slice
for SliI=WhichSlices
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);

    [U,S,sccmtx] = svd(nukData.','econ');
    
    sccmtxS(:,:,SliI)=sccmtx;
    
    ncc=13;
    
    CurSens=SensX(:,:,:,:,SliI);
%     weights=ones(gsize(CurSens,[1 2 4]));
%     ESP = ESPIRiT(CurSens,weights);

    SensCC=permute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2 4])),[2 3 1 4]);
    SensCCS(:,:,:,:,SliI)=SensCC;
% weightsCC=ones(gsize(CurSens,[1 2 4]));
% ESPCC = ESPIRiT(SensCC,weightsCC);
end
disp('SensCCS ok');
%%
for SliI=WhichSlices
    CurOurPSli=[BaseOutLoc BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
    mkdir(CurOurPSli);
    
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);
    
    DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
    
    DataPCC=double(permute(DataCC,[3 2 1]));
    
    tmp=squeeze(DataPCC);
%     RealDataFacS(SliI)=grmss(AllData(5656:34:6767,:,:))/grmss(tmp);
    CurIDataV=Row(tmp)*100; %*RealDataFacS(SliI);
    CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
    
    Data=repmat(single(CurIDataVR),[16 1]);
    RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    save(RealDataFN,'Data');
end
disp('Save slices real data for NN');
% save([BaseOutLoc BaseFN filesep 'RealDataFacS.mat'],'RealDataFacS');
%% All reps
CurRealDataP=[ScanP BaseFN filesep 'RealData' filesep];
mkdir(CurRealDataP);
for SliI=WhichSlices
    for r=1:nReps
        disp([SliI r]);
        nukData=ADataIsPy(:,:,SliI,r).';
        nukData=nukData(:,3:end);
        nukDataCC=MultMatTensor(sccmtxS(:,1:ncc,SliI).',nukData);
        
        CurIDataV=Row(nukDataCC.')*100; % *RealDataFac;
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
        Data=repmat(single(CurIDataVR),[16 1]);
        RealDataFN=[CurRealDataP 'Sli' num2str(SliI) '_r' num2str(r,'%02d') '.mat'];
        %     RealDataFN=['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RealData/b_Ben14May_Sli5_r' num2str(r,'%02d') '.mat'];
        save(RealDataFN,'Data');
    end
end
disp('Saved slices real data for NN all reps');

CurRealDataOutP=[ScanP BaseFN filesep 'RealDataOut' filesep];
mkdir(CurRealDataOutP);
for SliI=WhichSlices
    mkdir([CurRealDataOutP 'Sli' num2str(SliI,'%02d')]);
end
%%
% MgcS=imresizeBySlices(Mg(:,:,(1:nSlices)+6),Sz2);
MgcS=imresizeBySlices(Mg,Sz2);
MskcS=MgcS>7e-5;

% B0M2S=-B0Q2;
% B0M2S(~MskcS)=0;
% 
% SymMsk=abs(B0M2S-gflip(B0M2S,1))>230;
% SymMsk(1:30,:,:)=true;
% B0M2S(SymMsk & B0M2S>150)=-20;

B0M2S=B0Q2;

% fgmontage(B0M2S,[-500 500])

% WhichSlices=1:nSlices;
WhichSlices=nSlices:-1:1;
%%
clear TSBS TSCS
%%
nTS=12;
for SliI=WhichSlices
    B0M2=B0M2S(:,:,SliI);
    Mgc=MgcS(:,:,SliI);
    
    AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
    MgcN=Mgc./grmss(Mgc);
    WE=Col(MgcN);
    % WE=Col(MgcN)*0+1;
    WeightedE=WE.*E;
    % Fessler time segmentation
    % nTS=7;
    clear ErrTS
    TS_Thresh=1e-5;
%     for nTS=13:15
%     nTS=12;
%         disp(nTS)
        FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
        TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
        TSC2=reshape(TSC,prod(Sz2),nTS);
        WTSC2=WE.*TSC2;
    %     tic
    %     TSB=(E.')/(TSC2.');% W in both sides
        TSB=(WeightedE.')/(WTSC2.');% W in both sides

    %     ErrTS(nTS)=grmss(E-TSC2*(TSB.')); %include W
        CurErrTS=grmss(WeightedE-WTSC2*(TSB.')); %include W
        ErrTSS(nTS,SliI)=CurErrTS;

        disp([datestr(now) ' Sli #' num2str(SliI) ' nTS ' num2str(nTS) ' err=' num2str(CurErrTS)]);
%         if(ErrTSS(nTS)<TS_Thresh)
%             disp(['Stopped at #TS=' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
%             break;
%         end
%     end
%     figure(87234);clf;plot(log10(ErrTS),'-*')
    TSBS(:,:,SliI)=TSB;
    TSCS(:,:,:,SliI)=TSC;
end
disp('ok TSBS TSCS');
%%
save([ScanP BaseFN filesep 'TSBS_TSCS.mat'],'TSBS','TSCS','ErrTSS');
disp('Saved TSBS TSCS');
%% GPU TS
for SliI=WhichSlices
    TSBF=squeeze(TSBS(IB_TimeInMs2,:,SliI)).';
    Sens=squeeze(SensX(:,:,:,1,SliI));
    
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

for SliI=WhichSlices
    disp(SliI);
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);
    % nukData=nukData.*T2SEstComp;
    
    DataP=nukData;

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

    FigH=4000;
    figure(FigH);close(FigH);

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
BaseOutLoc='/media/a/DATA/';
CurOurP=[BaseOutLoc BaseFN filesep];
mkdir(CurOurP);
save([CurOurP 'im_resS.mat'],'im_resS');
%% L1cgESPIRIT_CC on all:
for SliI=WhichSlices
    SensCC=squeeze(SensCCS(:,:,:,:,SliI));    
    TSBF=squeeze(TSBS(IB_TimeInMs2,:,SliI)).';
    GOP2S{SliI}=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,double(TSCS(:,:,:,SliI)));
    GOP_TSS{SliI}=ggpuNUFT_TS(GOP2S{SliI},double(TSBF));
    GOP_TS_ESPCCS{SliI}=ggpuNUFT_TS_MCnosum(GOP_TSS{SliI},[1 1 ncc]);
    
    ESPCCS{SliI} = ESPIRiT(SensCC,weightsCC);
end
disp('ok');
%%
for SliI=WhichSlices
    disp(['Sli #' num2str(SliI) ' ' datestr(now)]);
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);
    DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
    DataPCC=double(permute(DataCC,[3 2 1]));

    XOP = Wavelet('Daubechies_TI',4,6);
    splitWeight=0.5;
    lam=1e-5;
    nIterSplit=5e-4;
    % nIterSplit=15;
    [resL1ESPIRiTCCS(:,:,:,SliI)] = cgL1ESPIRiT(DataPCC, weightsCC*0, GOP_TS_ESPCCS{SliI}, ESPCCS{SliI}, 5,XOP,lam,splitWeight,nIterSplit);
end
disp('ok resL1ESPIRiTCCS');
%%
resL1ESPIRiTCCS1=squeeze(resL1ESPIRiTCCS(:,:,1,:));
fgmontage(resL1ESPIRiTCCS1,[0 7e-3])
title(['resL1ESPIRiT 2 maps, B0, with CC -> ' num2str(ncc)]);
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
        nukData=nukData(:,3:end);
        DataCC=(sccmtxS(:,1:ncc,SliI).'*nukData);
        DataPCC=double(permute(DataCC,[3 2 1]));
        
        XOP = Wavelet('Daubechies_TI',4,6);
        splitWeight=0.5;
        lam=1e-5;
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
%% Now calc on all for NN
nI=size(AllImWithPhaseComplexSingle,1);
nTrain=nI;
AllData=single(zeros(nTrain,nTrajP2,ncc));
AllData(1)=0+1i*eps;
%%
for SliI=WhichSlices
    CurOurPSli=[BaseOutLoc BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
    mkdir(CurOurPSli);

    StartTime=now;
    
    SensCC=squeeze(SensCCS(:,:,:,:,SliI));
    
    TSBF=squeeze(TSBS(IB_TimeInMs2,:,SliI)).';
    CurGOP2=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTrajAct,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,double(TSCS(:,:,:,SliI)));
    CurGOP_TS=ggpuNUFT_TS(CurGOP2,double(TSBF));
    GOP_TS_ESPCC=ggpuNUFT_TS_MCnosum(CurGOP_TS,[1 1 ncc]);
    
    for i=1:nTrain
        if(mod(i,100)==1)
            DeltaTime=now-StartTime;
            FullTime=DeltaTime*nTrain/i;
            EndEstTime=StartTime+FullTime;
            disp(['SliI:' num2str(SliI,'%02d') ', ' num2str(i,'%05d') ' Cur ' datestr(now) ' est finish ' datestr(EndEstTime)]);
        end
        CurI=double(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
        CurIWChCC=CurI.*SensCC(:,:,:,1);
        CurIData=GOP_TS_ESPCC*CurIWChCC;
        %     CurIData=(GOP_MC_CC*CurI).';
        AllData(i,:,:)=single(CurIData);
    end
    
    save([CurOurPSli 'AllDataCC' num2str(ncc) '.mat'],'-v7.3','AllData');
    disp(['Saved to ' CurOurPSli 'AllDataCC' num2str(ncc) '.mat']);
end
disp('Finished');
%%
nData=nI;
ChunkSize=100;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);
%%
for SliI=WhichSlices
    CurOurPSli=[BaseOutLoc BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
    clear AllData
    load([CurOurPSli 'AllDataCC' num2str(ncc) '.mat']);

    Msk=grmss(SensCCS(:,:,:,1,SliI),3)>0.01;
    MskP=permute(Msk,[3 1 2]);
    
%     TFDataP=[CurOurPSli 'TF/'];
    TFDataP=[CurOurPSli 'TFPhantom/'];
    mkdir(TFDataP);
    clear DataPR LabelsP
    StartTime=now;
    for k=1:numel(ChunkStartI)
        CurIs=ChunkStartI(k):ChunkEndI(k);
        CurChunkSize=numel(CurIs);
        
        AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:).*MskP;
        LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
        
        DataPR=NaN(CurChunkSize,nTrajAct*ncc*2);
        
        for i=1:CurChunkSize
            CurIData=squeeze(AllData(CurIs(i),:,:));
            CurIDataV=Row(CurIData); % Changed order in AllData!!!
            CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
            DataPR(i,:)=single(CurIDataVR);
        end
        
        % Data to TFRecords
        Data=single(DataPR);
        Labels=single(LabelsP);
        FNs=strcat(num2str(CurIs.','%05d'),'asd');
        save('~/HomeA/TF/CurChunk.mat','Data','Labels','FNs')
        
        system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRec.py ~/HomeA/TF/CurChunk.mat ' TFDataP]);
        %     if(mod(i,100)==1)
        DeltaTime=now-StartTime;
        FullTime=DeltaTime*numel(ChunkStartI)/k;
        EndEstTime=StartTime+FullTime;
        disp(['SliI: ' num2str(SliI) ' ' num2str(k,'%03d') ' Cur ' datestr(now) ' est finish ' datestr(EndEstTime)]);
        %     end
    end
end
%% Make folder structure in scratch disk
BaseScratchP='/media/a/H1/';
ScratchP=[BaseScratchP BaseFN filesep];
mkdir(ScatchP);
for SliI=WhichSlices
    mkdir([ScratchP 'Sli' num2str(SliI)]);
    fileattrib([ScratchP 'Sli' num2str(SliI)],'+w','a');
end
%% End Multislice
%%
AllData=load('/media/a/DATA/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/Sli08/AllDataCC13.mat');
AllDatax=AllData.AllData(1:10000,:,1:8);
save([ScratchP 'AllData_Sli' num2str(SliI) '.mat'],'-v7.3','AllDatax');
save([ScratchP 'AllData_Sli' num2str(SliI) '_1k.mat'],'-v7.3','AllDatax');
%%
NNN=100;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');

NNN=200;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');

NNN=500;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');

NNN=1000;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');

NNN=2000;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');

NNN=4000;
AllDatax=AllData.AllData(1:NNN,:,1:8);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
save([ScratchP 'AllData_Sli' num2str(SliI) '_' num2str(NNN,'%5d') '.mat'],'-v7.3','AllDatax');
%%
NNN=6000;
AllDatax=AllData.AllData(1:NNN,:,1:8);
% AllDatax=permute(AllDatax,[1 3 2]);
AllDatax=reshape(AllDatax,NNN,[]);
AllDatax=[real(AllDatax) imag(AllDatax)];
AllDatax=permute(AllDatax,[2 1]);
% save([ScratchP 'AllData_Sli' num2str(SliI) '_1k.mat'],'-v7.3','AllDatax');
save([ScratchP 'AllData_Sli' num2str(SliI) '_6k.mat'],'-v7.3','AllDatax');
save([ScratchP 'AllData_Sli' num2str(SliI) '_6k_nc.mat'],'-v7.3','-nocompression','AllDatax');
%

AllLh5=cat(4,real(AllImWithPhaseComplexSingle(1:NNN,:,:)),imag(AllImWithPhaseComplexSingle(1:NNN,:,:)));
AllLh5=permute(AllLh5,[4 3 2 1]);
save('/media/a/H1/First1kImWithPhaseComplexSingle_h5.mat','-v7.3','AllLh5');
%%
AllLh5=cat(4,real(AllImWithPhaseComplexSingle),imag(AllImWithPhaseComplexSingle));
AllLh5=permute(AllLh5,[4 3 2 1]);
save('/media/a/H1/AllImWithPhaseComplexSingle_h5.mat','-v7.3','AllLh5');
%%
S=load('/media/a/DATA/180628_AK/meas_MID265_BP_fieldmap_5echosX_FID22460/Sens.mat');
S=S.SensB;
Msk=grmss(S(:,:,:,14,1),3)>0.01;
Msk=single(gflip(Msk,1:2));
Msk=repmat(Msk,[1 1 2]);
save('/media/a/DATA/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/Sli08/Msk.mat','Msk');

%%
A=rand(2,3,4)
%%
B=permute(A,[3 2 1]);
save('/media/a/H1/A.mat','-v7.3','B');

%%



























%%
resL1ESPIRiT1e5=resL1ESPIRiT;
%%
































MC=ESP*im_res;
OneCh=ESP'*MC;


ImgMC=GOP_TS_ESP'*DataP2;

size(ImgMC)

size(Datax)

DMC=GOP_TS_ESP*ImgMC;

%%
[resL1ESPIRiT] = cgL1ESPIRiT(DATAc, resESPIRiT*0, FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);

%%
size(DataP)


II=GOP_TS_ESP'*DataP2;

DD=GOP_TS_ESP*II;

AOdd
ESP = ESPIRiT(maps,weights);


[resL1ESPIRiT] = cgL1ESPIRiT(DATAc, resESPIRiT*0, FT, ESP, 5,XFM,lambda,splitWeight,nIterSplit);


% ESPIRiT CG reconstruction with 1 map
disp('Performing SENSE reconstruction from 1 set of maps')
SNS = ESPIRiT(maps(:,:,:,end));

%% Reconsturctions
% ESPIRiT CG reconstruction with soft-sense and 1 sets of maps
disp('Performing L1-ESPIRiT reconstruction from 2 maps')
tic
[resL1ESPIRiT] = cgL1ESPIRiT(DATAc, resESPIRiT*0, FT, ESP, nIterCG,XOP,lambda,splitWeight,nIterSplit);
toc

XOP = Wavelet('Daubechies_TI',4,6);
FT = p2DFT(mask,[sx,sy,Nc]);

disp('Performing ESPIRiT reconstruction from 2 maps')
tic; [reskESPIRiT, resESPIRiT] = cgESPIRiT(DATAc,ESP, nIterCG*3, 0.01,DATAc*0); toc

disp('Performing SENSE reconstruction from 1 set of maps')
tic;[reskSENSE, resSENSE] = cgESPIRiT(DATAc, SNS, nIterCG*3,0.01 ,DATAc*0);toc

%%







%% Fit B0:
Sz1=Sz2;
Area=[5 5];
Area=[7 7];

I=Mgc;
B0Real=B0M2;

TmpMaskA=imfillholesBySlices(I>5e-5);
TmpMaskE=TmpMaskA;
% se = strel('disk',1);
% TmpMaskE=imerode(TmpMaskA,se);

Z=zeros(Sz1);
MidSz=floor(Sz1/2);
Reg=false(Sz1);
Reg((MidSz(1)-Area(1)):(MidSz(1)+Area(1)),(MidSz(2)-Area(2)):(MidSz(2)+Area(2)))=true;
F=find(Reg);
n=numel(F);
ParamToField=@(x) abs(fft2c(PutValsInIdxs(Z,x(1:n)+1i*x(n+(1:n)),F)));
GetValsInLoc=@(A,B) A(B);
FMsk=find(TmpMaskE);


ToAdd=-min(B0Real(:))+20;


% CostFunc=@(x) grmss(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+200,FMsk));

W=GetValsInLoc(I,FMsk).';
CostFunc=@(x) W*abs(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+ToAdd,FMsk));

options = optimset('Display','iter');
B0M=(B0Real+ToAdd).*TmpMaskE;
B0M(~isfinite(B0M))=0;
C=ifft2c(B0M);
X0=[real(C(Reg)).' imag(C(Reg)).'];
[BestX, BestCost]=fminsearch(CostFunc,X0,options);

% fgmontage(ParamToField(BestX)-ToAdd,[-200 300]);colorbar
% fgmontage((B0Real).*TmpMaskE,[-200 300]);colorbar
% 
% fgmontage(ParamToField(BestX).*TmpMaskA,[0 500]);colorbar
% 
% fgmontage(ParamToField(BestX),[0 500]);colorbar

B0RealEx=ParamToField(BestX)-ToAdd;
%%
fgmontage(B0Real,[-400 300])
fgmontage(B0RealEx,[-400 300])

fgmontage(B0Real-B0RealEx,[-100 100])
%%
fo = fitoptions('Weights',z*0+1);
    
load franke
sf = fit([x, y],z,'poly23',fo)
figure;plot(sf,[x,y],z)
%%
% ToFit=B0Q2(:,:,SliI);
ToFit=B0M2;
W=abs(Mgc);
[x,y]=meshgrid(1:128,1:128);
xF=x(:);
yF=y(:);
wF=W(:);
zF=ToFit(:);

DForFit=5;
% 5: 1.14 + 4.19 sec
x=xF(1:DForFit:end);
y=yF(1:DForFit:end);
z=zF(1:DForFit:end);
w=wF(1:DForFit:end);
fo = fitoptions('Weights',w,'Method','LowessFit','Span',0.01);
tic
sf = fit([x, y],z,'lowess',fo)
% sf = fit([x, y],z,fo)
toc
% figure;plot(sf,[x,y],z)
tic
X=sf([xF(1:1:end),yF(1:1:end)]);
toc
X2=reshape(X,Sz2);

fgmontage(ToFit,[-800 400])
fgmontage(X2,[-800 400])