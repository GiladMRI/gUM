DatasetP='~/HomeA/HCPDataset/';
% OutP='/home/a/TF/srezx/dataset';
D=dir([DatasetP '*.mat']);
DFNs={D.name};
Ni=numel(DFNs);
% Ni=100;
HCPData=NaN(256,256,Ni);
for i=1:Ni
    if(mod(i,1000)==1), disp(i), end
    tmp=load([DatasetP DFNs{i}]);
    HCPData(:,:,i)=tmp.CurIc;
end

ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),1);
ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),2);
ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=permute(HCPData(:,:,ToFlip),[2 1 3]);
%% crop to 128*128
Szk=[128 128];
x=(256-Szk(1))/2;
y=(256-Szk(2))/2;
%%
Idxs=getKrandomSamples(Ni,20000);
%% collect couples
DataPR=permute(HCPData(x:(x+Szk(1)-1),y:(y+Szk(2)-1),Idxs),[3,1,2]);
DataPR=RepDotMult(DataPR,1./max(max(abs(DataPR),[],2),[],1));
LabelsP=DataPR;
%%
dataFoldName='dataSAE';
mkdir(['/home/a/TF/srez/' dataFoldName '/']);
system(['chmod -R 777 ' '/home/a/TF/srez/' dataFoldName '/']);
[status,msg,msgID] = fileattrib(['/home/a/TF/srez/' dataFoldName filesep],'+w','a');

nData=size(DataPR,1);
ChunkSize=100;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
        
    Data=single(DataPR(CurIs,:,:,:));
    Labels=single(LabelsP(CurIs,:,:,:));
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('/home/a/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/b.sh /home/a/TF/Mat2TFRec.py /home/a/TF/CurChunk.mat /home/a/TF/srez/' dataFoldName '/ '])
end
%%
TBaseP='/home/a/TF/srez/';

DataH=size(Data,2);
DataW=size(Data,3);
DataCh=size(Data,4);
LabelsH=size(Labels,2);
LabelsW=size(Labels,3);
LabelsCh=size(Labels,4);

ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
  'dataset',dataFoldName,'SessionNameBase','kKick','learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0.0001,'WL2_Lambda',0.0001,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'Mode','SAE','DiscStartMinute',5);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
%% Run several
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=5;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);
%% Visualization
Q=load('/home/a/TF/srez/kKick_dataSAE__2018-05-03_11-57-40_train/TrainSummary_005340.mat');


Q=load('/home/a/TF/srez/kKick_dataSAE__2018-05-03_15-01-09_train/TrainSummary_004310.mat');
WL=double(Q.gene_GEN_L007_C2D_weight_0);
WLN=RepDotMult(WL,1./grms(WL,1:2));