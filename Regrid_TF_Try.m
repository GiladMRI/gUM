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

% figure;plot(linspace(-kMax,kMax,osN+1),1,'.k');hold on;plot(C,1,'.r')
%% Load HCP dataset
% DatasetP='/media/a/fd6c4d95-0129-4d66-b27b-4267b50519fd/HCPDataset/';
DatasetP='/media/a/H1/HCPDataset/';
% OutP='/home/a/TF/srezx/dataset';
D=dir([DatasetP '*.mat']);
DFNs={D.name};
Ni=numel(DFNs);
% Ni=100;
HCPData=int16(zeros(256,256,Ni));
for i=1:Ni
    if(mod(i,1000)==1), disp(i), end
    tmp=load([DatasetP DFNs{i}]);
    tmp=int16(tmp.CurIc);
    HCPData(:,:,i)=tmp;
end
%%
save('/media/a/H1/HCPData_256x256_int16.mat','-v7.3','HCPData')
%%
% RandFx=rand(Ni,1);
% RandFy=rand(Ni,1);
% RandFp=rand(Ni,1);
% %%
% ToFlip=RandFx>0.5;
% HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),1);
% ToFlip=RandFy>0.5;
% HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),2);
% ToFlip=RandFp>0.5;
% HCPData(:,:,ToFlip)=permute(HCPData(:,:,ToFlip),[2 1 3]);
%%
save('/media/a/H1/HCPDataR_256x256_int16.mat','-v7.3','HCPData')
%%
CurSens=SensX(:,:,:,SliI);
% CurSens=squeeze(SensB(:,:,:,12));
ncc=size(CurSens,3);
SensCC=CurSens;
sccmtx=eye(ncc);
%%
DataSmallGoodChannels=CombineDims(CurSens,[1 2]);

[U,S,sccmtx] = svd(DataSmallGoodChannels,'econ');
ncc=8;

SensCC=permute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2])),[2 3 1]);
%% Try BART
nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],CurBartTraj(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(CurSens,[1 2 4 3]));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],CurBartTraj(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(CurSens,[1 2 4 3]));

Lambda=1e-9;
Rec=RecIfTVs(Lambda);
% Rec=RecIfWs(Lambda);

ShowAbsAngle(Rec)
%% Try BART on CC
nukData=ADataIsPy(:,:,SliI,3).';
nukData=nukData(:,3:end);
nukDataCC=MultMatTensor(sccmtx(:,1:ncc).',nukData);
nTrajAct=size(nukDataCC,2);
CurBartTrajAct=CurBartTraj(:,1:nTrajAct);
%%
RecIfTVsCC=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],CurBartTrajAct, permute(nukDataCC,[3 2 4 1]), permute(SensCC,[1 2 4 3]));
RecIfWsCC=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],CurBartTrajAct, permute(nukDataCC,[3 2 4 1]), permute(SensCC,[1 2 4 3]));

Lambda=1e-9;
RecCC=RecIfTVsCC(Lambda);
% RecCC=RecIfWsCC(Lambda);

ShowAbsAngle(RecCC)
%% Generate Idx mat of neighbors
nNeighbors=12;
NMap=NaN(osN,osN,nNeighbors);
for i=1:osN
    for j=1:osN
        CurLoc=C([i j]).';
        D=CurBartTrajAct(1:2,:)-CurLoc;
        R=grmss(D,1);
        [~, Idx]=sort(R);
        NMap(i,j,:)=Idx(1:nNeighbors);
    end
end
NMapC=RepDotMult(ones(size(NMap))*nTrajAct,permute(0:ncc-1,[1 4 3 2]));
NMapC=NMapC+repmat(NMap,[1 1 1 ncc]);
NMapCX=CombineDims(NMapC,[3 4]);
NMapCR=cat(3,NMapCX,NMapCX+nTrajAct*ncc)-1;

NMapCRC=cat(4,NMapCX,NMapCX+nTrajAct*ncc)-1;
% small test
CurIData=GOP_MC*CurI;
CX=gmat2cell(CurIData,1);
clear CurIDataC
for c=1:ncc
    CurIDataC(:,:,:,c)=CX{c}(NMap);
end
CurIDataCX=CombineDims(CurIDataC,[3 4]);
CurIDataCXR=single(cat(3,real(CurIDataCX),imag(CurIDataCX)));

CurIDataV=Row(CurIData.');
CurIDataVR=[real(CurIDataV) imag(CurIDataV)];

CurIDataCXR2=CurIDataVR(NMapCR);

grmss(CurIDataCXR2-CurIDataCXR)

% save('~/HomeA/TF/NMapIndTesta.mat','NMapCR','CurIDataVR','CurIDataCXR');

NMapFN='/media/a/H1/NMapIndCB0.mat';
save(NMapFN,'NMapCR');
% save('~/HomeA/TF/NMapIndCB0.mat','NMapCRC');
%% GPU TS
Sz2=gsize(SensCC,1:2);

TSB=ones(1,nTrajAct);
TSC=ones(Sz2);
% w=ones(nTrajP2,1);
osf = 1.5; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

GOP_MC = ggpuNUFT_TS_MCx(CurBartTrajAct,Sz2,osf,wg,sw,TSB,TSC,SensCC);

x = randn(Sz2) + 1j*randn(Sz2);
y = randn([ncc nTrajAct]) + 1j*randn([ncc nTrajAct]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))

%%
nI=size(HCPData,3);
% P=randperm(nI);
%%
nTrain=nI;
nTrain=20000;
% nTrain=10000;
ChunkSize=100;
%%
CurIDataC=NaN([gsize(NMap,1:3) ncc]);
%%
% SzTrg=[kMax kMax]*2;
SzTrg=[128 128];
%%
BaseDatasetsP='/media/a/H1/TFDatasets/';
dataFoldName='HCP128x128ImagesWithPhase';
mkdir([BaseDatasetsP dataFoldName '/']);
system(['chmod -R 777 ' BaseDatasetsP dataFoldName '/']);
[status,msg,msgID] = fileattrib([BaseDatasetsP dataFoldName filesep],'+w','a');

nData=nTrain;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);
%% Just the forward NUFT part
AllImWithPhaseComplexSingle=single(zeros([nTrain SzTrg]));
AllData=single(zeros([nTrain-20000 ncc nTrajAct]));
RandVecsForPhase=rand(nTrain,11);
% RandVecsForPhase(20001:nTrain,:)=rand(nTrain-20000,11);
%%
for i=1:nTrain
    if(mod(i,100)==1), disp(i),end
    CurI=double(HCPData(:,:,i));
    CurI=imresize(CurI,Sz2);
    Mx=max(1,max(CurI(:)));
    CurI=CurI/Mx;
    GPhi=GenerateRandomSinPhase(Sz2,5,0.1,RandVecsForPhase(i,:));
    CurI=CurI.*GPhi;
%     AllImWithPhaseFullSz(i,:,:)=
    
    if(all(SzTrg==Sz2))
        CurITrg=CurI;
    else
%     F=fft2cg(CurI);
%     FC=crop(F,SzTrg);
%     CurITrg=ifft2c(FC);
    end
    
    AllImWithPhaseComplexSingle(i,:,:)=single(CurITrg);
end
    %%
% end
% %%
for i=1:nTrain
    if(mod(i,100)==1), disp(i),end
    CurI=double(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
    CurIData=GOP_MC*CurI;
    AllData(i-20000,:,:)=CurIData;
end
%%
save('/media/a/H1/All32kImWithPhaseComplexSingleX128x128.mat','-v7.3','AllImWithPhaseComplexSingle');
save('/media/a/H1/AllDataHCPNoB08ChB0.mat','-v7.3','AllData');
save('/media/a/H1/Rest20kPlusDataHCPNoB08ChB0.mat','-v7.3','AllData');
%%
AllIm128x128MagSingle=abs(AllImWithPhaseComplexSingle);
save('/media/a/H1/AllIm128x128MagSingle.mat','-v7.3','AllIm128x128MagSingle');

First10kIm128x128MagSingle=AllIm128x128MagSingle(1:10000,:,:);
save('/media/a/H1/First10kIm128x128MagSingle.mat','-v7.3','First10kIm128x128MagSingle');

First3kIm128x128MagSingle=AllIm128x128MagSingle(1:3000,:,:);
save('/media/a/H1/First3kIm128x128MagSingle.mat','-v7.3','First3kIm128x128MagSingle');

save('/media/a/H1/First3kIm128x128MagSinglex.mat','First3kIm128x128MagSingle');
%%
AllData=load('AllDataHCPNoB08Ch.mat');AllData=AllData.AllData;
% AllImWithPhaseComplexSingle=load('AllImWithPhaseComplexSingle.mat');AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle.AllImWithPhaseComplexSingle;
%% save Images
clear LabelsP
for k=201:205 %numel(ChunkStartI)
    CurIs=(ChunkStartI(k):ChunkEndI(k)).';
    CurChunkSize=numel(CurIs);
    
    AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:);
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    Labels=single(LabelsP);
    
    FNs=strcat(num2str(CurIs,'%05d'),'asd');
    
    CurIs=single(CurIs);
    save('/media/a/H1/CurChunk.mat','Labels','FNs','CurIs')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRecx.py /media/a/H1/CurChunk.mat ' BaseDatasetsP dataFoldName '/ '])
end
%%
% k=1
% 
% CurIs=(1:3).';
% CurChunkSize=numel(CurIs);
% 
% AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:);
% LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
% Labels=single(LabelsP);
% 
% FNs=strcat(num2str(CurIs,'%05d'),'asd');
% 
% CurIs=single(randi(7,3,1));
% save('/media/a/H1/CurChunk.mat','Labels','FNs','CurIs')
% 
% system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRecx.py /media/a/H1/CurChunk.mat ' BaseDatasetsP dataFoldName '/ '])

%%
BaseDatasetsP='/media/a/H1/TFDatasets/';
dataFoldName='Traj5118Ben4minASL32ch';
mkdir([BaseDatasetsP dataFoldName '/']);
system(['chmod -R 777 ' BaseDatasetsP dataFoldName '/']);
[status,msg,msgID] = fileattrib([BaseDatasetsP dataFoldName filesep],'+w','a');
%% save Data
clear DataPR
for k=201:205 %numel(ChunkStartI)
    CurIs=(ChunkStartI(k):ChunkEndI(k)).';
    CurChunkSize=numel(CurIs);
    
    CurData=AllData(CurIs-20000,:,:);
    DataPR=cat(4,single(real(CurData)),single(imag(CurData)));
    
    FNs=strcat(num2str(CurIs,'%05d'),'asd');
    CurIs=single(CurIs);
    save('/media/a/H1/CurChunk.mat','DataPR','FNs','CurIs')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRecx.py /media/a/H1/CurChunk.mat ' BaseDatasetsP dataFoldName '/ '])
end
%%
%% save Data
clear DataPR
for k=20001:20500 %numel(ChunkStartI)    
    CurData=squeeze(AllData(k-20000,:,:));
    DataPR=cat(3,single(real(CurData)),single(imag(CurData)));
    
    CurFN=strcat(num2str(k,'%05d'),'asd');
    save([BaseDatasetsP dataFoldName filesep CurFN '.mat'],'DataPR','k')
    disp(k)
end
%%












%%
for i=1:nTrain
    if(mod(i,100)==1), disp(i),end
    CurI=double(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
    CurI=imresize(CurI,Sz2);
    Mx=max(1,max(CurI(:)));
    CurI=CurI/Mx;
    GPhi=GenerateRandomSinPhase(Sz2,5,0.1,RandVecsForPhase(i,:));
    CurI=CurI.*GPhi;
%     AllImWithPhaseFullSz(i,:,:)=
    
    F=fft2cg(CurI);
    FC=crop(F,SzTrg);
    CurITrg=ifft2c(FC);
    AllImWithPhaseComplexSingle(i,:,:)=CurITrg;
% end
% %%
% for i=1:nTrain
%     CurI=double(squeeze(AllImWithPhaseComplexSingle(i,:,:)));
    CurIData=GOP_MC*CurI;
    AllData(i,:,:)=CurIData;
end

%% save as complex
clear DataPR LabelsP
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:);
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    
    DataPR=NaN(CurChunkSize,nTrajAct*ncc*2);

    for i=1:CurChunkSize
        CurIData=squeeze(AllData(CurIs(i),:,:));
        CurIDataV=Row(CurIData.');
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        DataPR(i,:)=single(CurIDataVR);
    end
    
    % Data to TFRecords
    Data=single(DataPR);
    Labels=single(LabelsP);
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('~/HomeA/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRec.py ~/HomeA/TF/CurChunk.mat ~/HomeA/TFDatasets/' dataFoldName '/ '])
%     disp(ChunkEndI(k));
end
%%
clear DataPR LabelsP
for k=numel(ChunkStartI):numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:);
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    
    DataPR=NaN(CurChunkSize,nTrajAct*ncc*2);

    for i=1:CurChunkSize
        CurIData=squeeze(AllData(CurIs(i),:,:));
        CurIDataV=Row(CurIData.');
        CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        DataPR(i,:)=single(CurIDataVR);
    end
    
    % Data to TFRecords
    Data=single(DataPR);
    Labels=single(LabelsP);
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('~/HomeA/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRec.py ~/HomeA/TF/CurChunk.mat ~/HomeA/TFDatasets/' dataFoldName '/ '])
%     disp(ChunkEndI(k));
end
%%
Msk=single(grmss(SensCC,3)*sqrt(ncc)>0.5);

F=fft2cg(Msk);
FC=crop(F,SzTrg);
Msk=ifft2c(FC);
Msk=single(repmat(abs(Msk)>0.5,[1 1 2]));
        
save('~/HomeA/TF/Msk.mat','Msk');
%%
CurIDataV=Row(nukDataCC.')*60;
CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
Data=repmat(single(CurIDataVR),[16 1]);
RealDataFN='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RealData/b_Ben14May_Sli5.mat';
save(RealDataFN,'Data');
%% 
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    AllData=NaN([CurChunkSize,gsize(NMap,1:2),nNeighbors*ncc]);
    AllLabel=NaN([CurChunkSize SzTrg]);
    CurIDataC=NaN([gsize(NMap,1:3) ncc]);

    for i=1:CurChunkSize
        CurI=HCPData(:,:,P(CurIs(i)));
        CurI=imresize(CurI,Sz2);
        Mx=max(1,max(CurI(:)));
        CurI=CurI/Mx;
        GPhi=GenerateRandomSinPhase(Sz2,5,0.1);
        CurI=CurI.*GPhi;

        F=fft2cg(CurI);
        FC=crop(F,SzTrg);
        CurITrg=ifft2c(FC);

        CurIData=GOP_MC*CurI;
        C=gmat2cell(CurIData,1);
        for c=1:ncc
            CurIDataC(:,:,:,c)=C{c}(NMap);
        end
        CurIDataCX=CombineDims(CurIDataC,[3 4]);
        AllData(i,:,:,:)=CurIDataCX;
        AllLabel(i,:,:)=CurITrg;
    end
    % collect couples
    DataPR=cat(4,single(real(AllData)),single(imag(AllData)));
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    
    % Data to TFRecords
    Data=single(DataPR);
    Labels=single(LabelsP);
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('~/HomeA/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/HomeA/b.sh ~/HomeA/TF/Mat2TFRec.py ~/HomeA/TF/CurChunk.mat ~/HomeA/TFDatasets/' dataFoldName '/ '])
%     disp(ChunkEndI(k));
end
%%
TBaseP='~/HomeA/TF/srez/';

DataH=size(NMap,1);
DataW=size(NMap,2);
DataCh=size(NMap,3)*ncc*2;
LabelsH=SzTrg(1);
LabelsW=SzTrg(2);
LabelsCh=2;

ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
  'dataset',dataFoldName,'learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0,'WL2_Lambda',0,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'DiscStartMinute',500,...
  'ShowRealData',1,'CmplxBias',0,...
  'Mode','RegridTry1C2',...
  'SessionNameBase','RegridTry1C2',...
  'NMAP_FN',NMapFN,...
  'RealDataFN',RealDataFN);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');

%% from HandleNewData
%%
B0Q2=B0Q(:,:,6+(1:12));
B0Q2=gflip(B0Q2,1:2);
%% All B0 effects across time
Mgc=imresizeBySlices(gflip(Mg(:,:,SliI+6),1:2),Sz2);
Mskc=Mgc>7e-5;

B0M2=-B0Q2(:,:,SliI);
B0M2(~Mskc)=0;

TimeInMs2=(0:nTrajAct-1)*2.5/1e3;
%%
% B0M2=B0RealEx;

[U_TimeInMs2, IA_TimeInMs2, IB_TimeInMs2]=unique(TimeInMs2);
nU_TimeInMs2=numel(U_TimeInMs2);

AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);

MgcN=Mgc./grmss(Mgc);
WE=Col(MgcN)*0+1;

WeightedE=WE.*E;
%% Fessler time segmentation
% nTS=7;
clear ErrTS
TS_Thresh=1e-5;
for nTS=15:15
    disp(nTS)
    FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    
    TSC2=reshape(TSC,prod(Sz2),nTS);
    WTSC2=WE.*TSC2;
    tic
%     TSB=(E.')/(TSC2.');% W in both sides
    TSB=(WeightedE.')/(WTSC2.');% W in both sides
    toc
%     ErrTS(nTS)=grmss(E-TSC2*(TSB.')); %include W
    ErrTS(nTS)=grmss(WeightedE-WTSC2*(TSB.')); %include W
    if(ErrTS(nTS)<TS_Thresh)
        disp(['Stopped at #TS=' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
        break;
    end
end
figure;plot(log10(ErrTS),'-*');title('Time Segmentation Error');xlabel('number of time segments');ylabel('Error');
%% GPU TS
% Sens=imresizeBySlices( squeeze(SensP2),Sz2);
osf = 1.5; % oversampling: 2 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

% TSC_MB=ones([NTrg1 1 nBands]);

% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),NTrg1,osf,wg,sw,ones(1,size(BARTTrajMS,2),nBands),TSC_MB,SensX(:,:,:,SliI));

NTrg1=Sz2;
nBands=1;
nShots=paramLongInterleaves;

TSBF=TSB(IB_TimeInMs2,:).';
% Sens=squeeze(SensX(:,:,:,SliI));

GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(CurBartTrajAct,NTrg1,osf,wg,sw,TSBF,TSC,SensCC);
GOP_MC=GOP_MCSMBMS;
nTrajP2=nU_TimeInMs2;
%%
% GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);

x = randn(Sz2) + 1j*randn(Sz2);
y = randn([ncc nTrajP2]) + 1j*randn([ncc nTrajP2]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))

%%
if(nShots==1)
    TVOP=TVOP_MSlice;
else
    TVOP=TVOP_MTC_W([1 1 0 1e1]);
end

nukData=ADataIsPy(:,:,SliI,13).';
nukData=nukData(:,3:end);
nukDataCC=MultMatTensor(sccmtx(:,1:ncc).',nukData);

DataP=nukDataCC;

DataP=squeeze(AllData(3843,:,:));

% AOdd = GOP_MCSMBMS;
AOdd = GOP_MC;

% TVW=0.1;
TVW=1e-5;

param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',1,'TV',TVOP,'xfmWeight',0),init);

param.data =     DataP;

param.pNorm=1;
% param.TVWeight=1e-7;
% param.TVWeight=0.00001;
% param.TVWeight=0.0001;
% param.TVWeight=0.01*2;

nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
% res=resC;
% res=resB;
if(nShots>1)
    res=repmat(resA,[1 1 1 nShots]);
    res=res+randn(size(res))*max(abs(resA(:)))/20;
end
% max(abs(resA(:)))


FigH=4000;
figure(FigH);close(FigH);

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
param.Verbose=false;
RunFnlViewAmp=2;
clear ObjConv Score
%
for n=1:nfnlCgIters
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