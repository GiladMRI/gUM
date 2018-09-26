% SMASH
AllImWithPhaseComplexSingle=load('AllImWithPhaseComplexSingle.mat');AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle.AllImWithPhaseComplexSingle;
AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle(:,7:end-6,:);
Szx=gsize(AllImWithPhaseComplexSingle,2:3);
%%
Brain8Ch=load('brain_8ch.mat');
%%
Sens8=RunESPIRiTForSensMaps(fft2cg(Brain8Ch.DATA),12,Szx);
ShowAbsAngle(Sens8)
%%
i=345;
CurI=squeeze(AllImWithPhaseComplexSingle(i,:,:));
CurIS=CurI.*Sens8;
F=fft2cg(CurIS);
F2=F(:,1:2:end,:);
FN=cat(4,F2(1:end-2,1:end-1,:),F2(2:end-1,1:end-1,:),F2(3:end,1:end-1,:),F2(1:end-2,2:end,:),F2(2:end-1,2:end,:),F2(3:end,2:end,:));
FNP=permute(FN,[1 2 4 3]);
FNPx=CombineDims(FNP,[3 4]);
FNPxC=cat(4,real(FNPx),imag(FNPx));
%%
nData=3000;
ChunkSize=100;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);

dataFoldName='dataNeighborhoodSMASH6F';
mkdir(['~/HomeA/TFDatasets/' dataFoldName '/']);
system(['chmod -R 777 ' '~/HomeA/TFDatasets/' dataFoldName '/']);
[status,msg,msgID] = fileattrib(['~/HomeA/TFDatasets/' dataFoldName filesep],'+w','a');

%% save as complex
clear DataPR LabelsP
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    AllLabel=AllImWithPhaseComplexSingle(CurIs,:,:);
    LabelsP=cat(4,single(real(AllLabel)),single(imag(AllLabel)));
    
    DataPR=NaN(CurChunkSize,98,55,48*2);

    CurIS=AllLabel.*permute(Sens8,[4 1 2 3]);
    F=fft1cg(fft1cg(CurIS,2),3);
    F2=F(:,:,1:2:end,:);
    FN=cat(5,F2(:,1:end-2,1:end-1,:),F2(:,2:end-1,1:end-1,:),F2(:,3:end,1:end-1,:),F2(:,1:end-2,2:end,:),F2(:,2:end-1,2:end,:),F2(:,3:end,2:end,:));
    FNP=permute(FN,[1 2 3 5 4]);
    FNPx=CombineDims(FNP,[4 5]);
    FNPx=ifft1cg(FNPx,2);
    DataPR=cat(4,real(FNPx),imag(FNPx));

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

DataH=size(DataPR,2);
DataW=size(DataPR,3);
DataCh=size(DataPR,4);
LabelsH=Szx(1);
LabelsW=Szx(2);
LabelsCh=2;

ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
  'dataset',dataFoldName,'learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0,'WL2_Lambda',0,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'DiscStartMinute',500,...
  'ShowRealData',0,'CmplxBias',0,...
  'Mode','SMASHTry1',...
  'SessionNameBase','SMASHTry1');

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');

%%

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCC_dataNeighborhoodSMASH6__2018-06-11_13-34-16_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_14-37-10_train/';
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_14-39-47_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_14-59-31_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_15-03-38_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_15-29-25_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_15-38-15_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_15-42-43_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_16-03-18_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_16-10-42_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_16-19-24_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_16-55-48_train/';
D=dir([VP 'Tra*.mat']);
Q=load([VP D(end).name]);
%%
Flds=fieldnames(Q);
SFlds=sort(Flds);
for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},45) num2str(size(Q.(SFlds{i})),'% 9d         ')]);
end
%%
FTH=double(Q.gene_GEN_L003_M2D_MC_weightR_0 + 1i*Q.gene_GEN_L003_M2D_MC_weightI_0);

GCC=double(Q.gene_GEN_L005_einsum_weightR_0 + 1i*Q.gene_GEN_L005_einsum_weightI_0);
%%
FSens8=fft1cg(Sens8,1);
FSens8=FSens8(2:end-1,:,:);

XSens8=MultMatTensor(FTH',FSens8);

GCCXSens8=squeeze(sum(XSens8.*permute(GCC,[1 4 2 3]),3));

GCCFSens8=MultMatTensor(inv(FTH'),GCCXSens8);

GCCISens8=ifft1cg(GCCFSens8,1);

ShowAbsAngle(GCCISens8)

%%
GCC=double(Q.gene_GEN_L004_einsum_weightR_0 + 1i*Q.gene_GEN_L004_einsum_weightI_0);

% FTH=double(Q.gene_GEN_L006_M2D_MC_weightR_0 + 1i*Q.gene_GEN_L006_M2D_MC_weightI_0);

%%
XSens8=ifft1cg(FSens8,1);

GCCXSens8=squeeze(sum(XSens8.*permute(GCC,[1 4 2 3]),3));
ShowAbsAngle(GCCXSens8)