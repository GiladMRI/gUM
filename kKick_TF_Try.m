DatasetP='/home/a/HCPDataset/';
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
%% Simply save all these images
dataFoldName='dataHCPImages';
mkdir(['/home/a/TF/srez/' dataFoldName '/']);
system(['chmod -R 777 ' '/home/a/TF/srez/' dataFoldName '/']);
[status,msg,msgID] = fileattrib(['/home/a/TF/srez/' dataFoldName filesep],'+w','a');

nData=size(HCPData,3);
ChunkSize=100;
ChunkStartI=1:ChunkSize:nData;
ChunkEndI=min(ChunkStartI+ChunkSize-1,nData);
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
        
    Data=single(permute(HCPData(:,:,CurIs),[3 1 2]));
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('/home/a/TF/CurChunk.mat','Data','FNs')
    
    system(['~/b.sh /home/a/TF/Mat2TFRec.py /home/a/TF/CurChunk.mat /home/a/TF/srez/' dataFoldName '/ '])
end
%% Run TF
TBaseP='/home/a/TF/srez/';

DataH=size(HCPData,1);
DataW=size(HCPData,2);
DataCh=1;
LabelsH=256;
LabelsW=256;
LabelsCh=4;

% DataH=size(Data,2);
% DataW=size(Data,3);
% DataCh=size(Data,4);
% LabelsH=size(Labels,2);
% LabelsW=size(Labels,3);
% LabelsCh=size(Labels,4);

ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
  'dataset',dataFoldName,'SessionNameBase','kKick','learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0.0001,'WL2_Lambda',0.0001,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'Mode','kKick','DiscStartMinute',5);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
%% Run several
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);




%%
%%
ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),1);
ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),2);
ToFlip=rand(Ni,1)>0.5;
HCPData(:,:,ToFlip)=permute(HCPData(:,:,ToFlip),[2 1 3]);
%%
HCPDatax=HCPData(:,:,getKrandomSamples(Ni,1000));

Nis=5000;
% HCPDatax=RepDotMult(HCPDatax,1./max(max(abs(HCPDatax),[],2),[],1));

%% crop to 120*180
Szk=[120 180];
x=(256-Szk(1))/2;
y=(256-Szk(2))/2;
%% collect couples
GPhasex=NaN([Szk 2]);
LabelsA=NaN([Szk Nis 2]);
for i=1:Nis
    if(mod(i,100)==1), disp(i), end
    Cur2Idxs=getKrandomSamples(Ni,2);
    CurIs=HCPData(x:(x+Szk(1)-1),y:(y+Szk(2)-1),Cur2Idxs);
    CurIs=RepDotMult(CurIs,1./max(max(abs(CurIs),[],2),[],1));
    for k=1:2
        GPhasex(:,:,k)=GenerateRandomBrokenPhase(Szk,5,1,2);
    end
    CurIs=CurIs.*GPhasex;
    
    LabelsA(:,:,i,:)=CurIs;
end
%% Data
DataA=NaN([Szk Nis]);
for i=1:Nis
    if(mod(i,100)==1), disp(i), end
    CurIs=LabelsA(:,:,i,:);
    FCurIs=fft1cg(CurIs,2);
    FCurIsc=FCurIs(:,1:120,1)+FCurIs(:,end-119:end,2);
%     Data(:,:,i)=ifft1cg(FCurIsc,2);
    DataA(:,:,i)=ifft1cg(padBoth(FCurIsc,30,2),2);
end
%%
GoodLabels=squeeze(gsum(~isfinite(LabelsA),[1 2 4])==0);
LabelsA=LabelsA(:,:,GoodLabels,:);
DataA=DataA(:,:,GoodLabels,:);
%%
DataP=permute(DataA,[3,1,2]);
LabelsP=permute(abs(LabelsA),[3,1,2,4]);
DataPR=cat(4,real(DataP),imag(DataP));
%%
dataFoldName='dataHCPkKick';
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
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'Mode','kKick','DiscStartMinute',5);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
%% Run several
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);
%%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=20;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-7; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=1e-6;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=0.00001;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=1e-6;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=1e-6;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=1e-6;
ParamsS.DiscStartMinute=90;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset=dataFoldName;
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=1e-8; % '0.0'; % 0.000001;
ParamsS.WL2_Lambda=1e-6;
ParamsS.DiscStartMinute=3;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);