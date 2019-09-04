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

ScanP='/media/a/DATA/ASLSubjData/S04/';
BaseFN='meas_MID149_gBP_VD11_U19_G35S155_FID23846';
RefFldMapP='/media/a/DATA/ASLSubjData/S04/meas_MID141_BP_fieldmap_5echosX_FID23838/';
% BaseFN='meas_MID164_gBP_VD11_U27_G35S155_2ADCs_FID23861';

% ScanP='/media/a/DATA/ASLSubjData/S03/';
% BaseFN='meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439';
% RefFldMapP='/media/a/DATA/ASLSubjData/S03/meas_MID265_BP_fieldmap_5echosX_FID22460/';

% ScanP='/media/a/DATA/ASLSubjData/S02/';
% BaseFN='meas_MID109_gBP_VD11_U19_4min_FID17944';
% RefFldMapP='/media/a/DATA/ASLSubjData/S02/meas_MID123_BP_fieldmap_5echosX_FID17958/';

% ScanP='/media/a/DATA/ASLSubjData/S01/';
% BaseFN='meas_MID409_gBP_VD11_U19_7ADCs_FID17798';
% RefFldMapP='/media/a/DATA/ASLSubjData/S01/meas_MID399_BP_fieldmap_4echos_FID17788/';

ScanP='/media/a/DATA/FC/';
% BaseFN='meas_MID203_gBP_ep2d_bold_multiecho_ASL_SMS_Spic_OutInOut_FID24589';
% RefFldMapP='/media/a/DATA/FC/meas_MID197_BP_fieldmap_5echosX_FID24583/';
BaseFN='meas_MID158_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_oblique_FID24551';
RefFldMapP='/media/a/DATA/FC/meas_MID156_BP_fieldmap_5echosX_oblique_FID24549/';

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
    SLoc(i,1)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dSag;
    SLoc(i,2)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dCor;
    SLoc(i,3)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dTra;
end

asSlice=AData.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[AData.hdr.Phoenix.sSliceArray.asSlice{:}];
end

nSlices=numel(AData.hdr.Phoenix.sSliceArray.asSlice);
% nSlices=numel(asSlice);
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

RotMat = transpose(Quat2RotMat(AData.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;

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
% if(MB>1)
if(isempty(spBW))
    spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{14};
    paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
    paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};
    paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
end
disp('Read data');
%%
save([ScanP BaseFN filesep 'Data.mat'],'ADataIsL','AData');
disp(['Saved ' [ScanP BaseFN filesep 'Data.mat']]);
%%
% save('StatusForTrajDesign1Band_2.mat');
%% given everything but VD and Acc
% OutInOut=false;
OutInOut=paramLongROSamples>10000;
if(OutInOut)
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/3,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);
else
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
        paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
end
WhichE=1;
if(OutInOut)
    BaseLen=size(kTraj,1)/3;
    WhichE=2;
    TrajIdxs=(1:BaseLen)+BaseLen*(WhichE-1);
    kTraj=kTraj(TrajIdxs,:);
    GradBuf=GradBuf(TrajIdxs,:);
end


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
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
%%
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'Traj'],[]) 
close(gcf);
save([ScanP BaseFN filesep 'Traj.mat'],'BARTTrajMS');
disp('Saved traj fig');
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
% dy=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;

kx=BARTTrajMS(1,:)*2*pi;
ky=BARTTrajMS(2,:)*2*pi;

dy=RotatedLocs(1,1)/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
modx=double(exp(1i*(dx*kx+dy*ky))');

ADataIsPx=reshape(ADataIsP,[paramLongROSamples,nChannels nSlices nReps]);
if(OutInOut)
    ADataIsPx=PartitionDim(ADataIsPx,1,3);
    ADataIsPx=ADataIsPx(:,:,:,:,WhichE);
end
ADataIsPy=ADataIsPx(1:size(BARTTrajMS,2),:,:,:,:);

ADataIsPy=RepDotMult(ADataIsPy,modx);

% IsOblique=isfield(AData.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal,'dCor');
% if(~IsOblique)
%     modx=double(exp(1i*(dx*kx+dy*ky))');
%     
%     ADataIsPx=reshape(ADataIsP,[paramLongROSamples,nChannels nSlices nReps]);
%     if(OutInOut)
%         nReps=3;
%         ADataIsPx=PartitionDim(ADataIsPx(:,:,:,1),1,3);
%     end
%     ADataIsPy=ADataIsPx(1:size(BARTTrajMS,2),:,:,:,:);
%     
%     ADataIsPy=RepDotMult(ADataIsPy,modx);
% else
%     Vnor=[AData.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dCor AData.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal.dTra];
%     Vpar=[Vnor(2) -Vnor(1)].';
%     for s=1:nSlices
%         Vpos=[AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor AData.hdr.Phoenix.sSliceArray.asSlice{SliI}.sPosition.dTra];
%         dyS(s)=Vpos*Vpar/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
%         modxS(s,:)=double(exp(1i*(dx*kx+dyS(s)*ky))');
%     end
%     ADataIsPx=reshape(ADataIsP,[paramLongROSamples,nChannels nSlices nReps]);
%     if(OutInOut)
%         nReps=3;
%         ADataIsPx=PartitionDim(ADataIsPx(:,:,:,1),1,3);
%     end
%     ADataIsPy=ADataIsPx(1:size(BARTTrajMS,2),:,:,:,:);
%     
%     ADataIsPy=ADataIsPy.*permute(modxS,[2 3 1 4]);
% end

% clear ADataIsPx ADataIsP ADataIx ADataIsL
disp('ok a')
%%
SliI=6;
%%
nukData=ADataIsPy(:,:,SliI,2).';

nukData=nukData(:,3:end);
nTraj=size(nukData,2);
TimeInMs2=(0:nTraj-1)*2.5/1e3;
T2SCompStr='';
T2SEstMs=20;
T2SEstDecay=exp(-TimeInMs2/T2SEstMs);
% T2SEstComp=exp(TimeInMs2/T2SEstMs);

Sz2=gsize(SensX,1:2);

BARTTrajAct=BARTTrajMS(:,1:size(nukData,2));
disp('ok');
%%
nukData=ADataIsPy(:,:,SliI,2).';
nukData=nukData(:,3:end);
% nukData=nukData.*T2SEstComp;
nukDataP=permute(nukData,[3 2 4 1]);

SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

% %% MB
% nukDataP=nukDataP(:,:,1);
% SensP=ones(Sz2);
disp('ok');
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
% close all;fgmontage(Mg(:,:,SliI));fgmontage(Rec);
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
gmontage(SensX(:,:,:,1,SliI));title('Fieldmap sens');
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
SliI=6;
%% All B0 effects across time
% Mgc=imresizeBySlices(gflip(Mg(:,:,SliI+6),1:2),Sz2);
% Mgc=imresizeBySlices(Mg(:,:,SliI+6),Sz2);
Mgc=imresizeBySlices(Mg(:,:,SliI),Sz2);
Mskc=Mgc>7e-5;
MskcE=imdilate(imfill(Mskc,'holes'),strel('disk',5,8));

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

B0M2(~MskcE)=0;

% B0M2=-B0M2;
% B0M2(B0M2>150)=-20;
% fgmontage(B0M2)
disp('ok');
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

disp('ok a');
%% Fessler time segmentation
% nTS=7;
clear ErrTS
TS_Thresh=1e-5;
for nTS=8:15
    disp(nTS)
    FesTimePoints=linspace(TimeInMs2(1),TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    TSC2=reshape(TSC,prod(Sz2),nTS);
    WTSC2=WE.*TSC2;
%     tic
%     TSB=(E.')/(TSC2.');% W in both sides
    TSB=(WeightedE.')/(WTSC2.');% W in both sides
    
%     ErrTS(nTS)=grmss(E-TSC2*(TSB.')); %include W
    ErrTS(nTS)=grmss(WeightedE-WTSC2*(TSB.')); %include W
    
    disp([datestr(now) ' nTS ' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
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

nukData=ADataIsPy(:,:,SliI,2).';
nukData=nukData(:,3:end);
% nukData=nukData.*T2SEstComp;
    
DataP=nukData;

% QQ=load('/media/a/H1/TFTSNUFTOut.mat');
% WithTSBR=double(QQ.Out)/100;
% DataP=WithTSBR.';

AOdd = GOP_MCSMBMS;
% AOdd = GOP_MC;

% TVW=0.1;
TVW=1e-5;


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
    xfmWeight = 1e-5;	% Weight for Transform L1 penalty
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
% res=double(MLN*grmss(im_resA)/grmss(MLN));
% res=im_resA;
% res=XFM*double(MLN*grmss(im_resA)/grmss(MLN));
% res=XFM*im_resA;
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

disp('ok im_res');

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
nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
%%
[U,S,sccmtx] = svd(nukData.','econ');
ncc=13;

nMapsToUse=size(SensX,4);
% nMapsToUse=1;

CurSens=SensX(:,:,:,1:nMapsToUse,SliI);
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
gprint(get(gcf,'Number'),[ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC' num2str(ncc)],[]) 
close(gcf);
save([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat'],'resL1ESPIRiTCC1');
disp(['Saved ' ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0_lam' num2str(lam) T2SCompStr '_CC.mat']);
%% Clean save
% BaseOutLoc='/media/a/DATA/';
SensMsk=single(grmss(SensCC(:,:,:),3)>0.01);
CurOurP=[ScanP BaseFN filesep];
mkdir(CurOurP);
CurOurPSli=[ScanP BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
mkdir(CurOurPSli);
save([CurOurPSli 'Sens.mat'],'CurSens');
save([CurOurPSli 'SensCC.mat'],'SensCC','sccmtx','SensMsk');
save([CurOurPSli 'B0TS.mat'],'TSBF','TSB','TSC','osf','wg','sw','Mgc','TimeInMs2');
save([CurOurPSli 'TrajAndRealData.mat'],'BARTTrajAct','DataP2','DataPCC');
disp('Saved');
%%
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
    CurOurPSli=[ScanP BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
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
nTS=11;
for SliI=WhichSlices
    B0M2=B0M2S(:,:,SliI);
    Mgc=MgcS(:,:,SliI);
    
    Mskc=Mgc>7e-5;
    MskcE=imdilate(imfill(Mskc,'holes'),strel('disk',5,8));
    B0M2(~MskcE)=0;

    
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
% load([ScanP BaseFN filesep 'TSBS_TSCS.mat'])
save([ScanP BaseFN filesep 'TSBS_TSCS.mat'],'TSBS','TSCS','ErrTSS');
disp('Saved TSBS TSCS');
%% save TSBF TSC Per slice
for SliI=WhichSlices
    YLbl=['Sli' num2str(SliI,'%02d')];
    
    SliIP=[ScanP BaseFN filesep YLbl filesep];
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
    CurOurP=[ScanP BaseFN filesep];
    CurOurPSli=[ScanP BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];
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

    ParamsSDefaults=getParamsStructFromFN('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/');
    ParamsS=ParamsSDefaults;
    ParamsS.SessionNameBase=['RegridTry3C2_7TS_' ScanP(end-3:end-1) '_Sli' num2str(SliI,'%02d')];
    ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    ParamsS.BaseTSDataP=CurOurPSli;
    ParamsS.BaseNUFTDataP=[ScanP BaseFN filesep];
    Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    
    system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
end
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

FigH=4000;
    figure(FigH);clf;
    
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
    nukData=nukData(:,3:end);
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
        nukData=nukData(:,3:end);
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
    system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
end