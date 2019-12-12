%%
DatasetP='/home/a/TF/srez/part of dataset/';
OutP='/home/a/TF/srezx/dataset';
D=dir([DatasetP '*.jpg']);
DFNs={D.name};
Ni=18000;
% Ni=100;
for i=1:Ni
    if(mod(i,1000)==1), disp(i), end
    TFData(:,:,:,i)=imread([DatasetP DFNs{i}]);
end
%% crop to 128*128, downscale to 64x64
x=(218-128)/2;
y=(178-128)/2;

TFDatax=TFData(x:(x+127),y:(y+127),:,:);

TFDatax1=mean(TFDatax,3);

StartPs=1:100:Ni;
TFDatac=NaN([64 64 1 Ni]); 
for i=1:numel(StartPs)
    TFDatac(:,:,:,(StartPs(i):(StartPs(i)+99)))=imresizeBySlices(TFDatax1(:,:,:,(StartPs(i):(StartPs(i)+99))),[64 64]);
    disp(i)
end
TFDataP=permute(TFDatac,[4 1 2 3]);

TFDatac=min(TFDatac,255);
TFDatacN=double(TFDatac)/255;
%
TFDatacM=grmss(TFDatacN,3);
%%
save('/home/a/TF/Imgs.mat','TFDatacM')
%%

Msk=imfillholesBySlices(grmss(SensCC,3)*sqrt(ncc)>0.3);

SWarning=warning('off','gpuNUFFT:adj:sens');
%%
GPhase=NaN(N,N,Ni);
for i=1:Ni
%     GPhase(:,:,i)=GenerateRandomBrokenPhase(N);
    GPhase(:,:,i)=GenerateRandomBrokenPhase(N,5,1,4);
    if(mod(i,1000)==0), disp(i),end
end
%%
ImgMskedPhased=NaN(N,N,Ni);
DataMskedPhased=NaN(ncc,nTrajP2,Ni);
for i=1:Ni
    ImgMskedPhased(:,:,i)=TFDatacM(:,:,i).*Msk.*GPhase(:,:,i);
    DataMskedPhased(:,:,i)=GOP_MC*ImgMskedPhased(:,:,i);
    if(mod(i,100)==0), disp(i),end
end
%%
dataFoldName='dataFaceP4C';
mkdir(['/home/a/TF/srez/' dataFoldName '/']);
system(['chmod -R 777 ' '/home/a/TF/srez/' dataFoldName '/']);

[status,msg,msgID] = fileattrib(['/home/a/TF/srez/' dataFoldName '/'],'+w','a');

ChunkSize=100;
ChunkStartI=1:ChunkSize:size(DataMskedPhased,3);
ChunkEndI=min(ChunkStartI+ChunkSize-1,size(DataMskedPhased,3));
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    CurDataMP=permute(DataMskedPhased(:,:,CurIs),[3 1 2]);
    CurDataMPFPR=cat(3,real(CurDataMP), imag(CurDataMP));
    
%     CurDataMP=permute(DataMskedPhased(:,:,CurIs),[2 1 3]);
%     CurDataMPF=reshape(CurDataMP,prod(gsize(CurDataMP,1:2)),size(CurDataMP,3));
%     CurDataMPFP=permute(CurDataMPF,[2 1]);
%     CurDataMPFPR=[real(CurDataMPFP) imag(CurDataMPFP)];

    Data=single(CurDataMPFPR);

    
    LabelMskedPhasedP=permute(ImgMskedPhased(:,:,CurIs),[3 1 2]);
    LabelsD=cat(4,real(LabelMskedPhasedP),imag(LabelMskedPhasedP));
    Labels=single(LabelsD(:,:,:,:));
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('/home/a/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/b.sh /home/a/TF/Mat2TFRec.py /home/a/TF/CurChunk.mat /home/a/TF/srez/' dataFoldName '/ '])
end
%%
% RealDataMP=repmat(DataP2.',[1 1 16]);
% RealDataMPF=reshape(RealDataMP,prod(gsize(RealDataMP,1:2)),size(RealDataMP,3));
% RealDataMPFP=permute(RealDataMPF,[2 1]);
% RealDataMPFPR=[real(RealDataMPFP) imag(RealDataMPFP)];

RealDataMPFP=repmat(permute(DataP2,[3 1 2]),[16 1 1]);
RealDataMPFPR=cat(3,real(RealDataMPFP),imag(RealDataMPFP));
Data=single(RealDataMPFPR*30);
save('/home/a/TF/srez/RealData/a.mat','Data')
%%
TBaseP='/home/a/TF/srez/';

DataH=size(Data,2);
DataW=size(Data,3);
ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',1,'LabelsH',64,'LabelsW',64,'channelsOut',2,...
  'dataset','dataKnee','SessionNameBase','Acc3_2','learning_rate_start',0.002,...
  'learning_rate_half_life',30,... % in minutes if <1000
  'summary_period',0.5,'checkpoint_period',20,...
  'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
  'WL1_Lambda',0.0001,'WL2_Lambda',0.0001,'noise_level',0.03,...
  'QuickFailureTimeM',3,'QuickFailureThresh',0.3,);

ParamsS=ParamsSDefaults;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
%% Run several
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda="0.0";
ParamsS.WL1_Lambda=0.000001;
ParamsS.WL2_Lambda=0.0001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');

sendTFMail(TBaseP,ParamsS,Txt);
%%
CurTrainDirbase=[ ParamsS.SessionNameBase '_' ParamsS.dataset '__'];
LastDir=getMostRecentDirWithBase(TBaseP,CurTrainDirbase);
RD=load([TBaseP LastDir filesep 'OnRealData.mat']);
%%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0.000001;
ParamsS.WL2_Lambda=0.0001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0.000001;
ParamsS.WL2_Lambda=0.0001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0.000001;
ParamsS.WL2_Lambda=0.0001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0.0000001;
ParamsS.WL2_Lambda=0.0001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda="0.0";
ParamsS.WL1_Lambda=0.0000001;
ParamsS.WL2_Lambda=0.00001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4C';
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=60;
ParamsS.WL1_Lambda=0.0000001;
ParamsS.WL2_Lambda=0.000001;
Txt=gStruct2txt(ParamsS,'/home/a/TF/Params.txt');
system('/home/a/RunTFForMatlab.sh');
sendTFMail(TBaseP,ParamsS,Txt);
%%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=120;
ParamsS.dataset='dataFaceP4';
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=10000;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');
%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4';
ParamsS.learning_rate_start=0.0005;
ParamsS.learning_rate_half_life=10000;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4';
ParamsS.learning_rate_start=0.001;
ParamsS.learning_rate_half_life=10000;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

%
%%
ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.learning_rate_start=0.0005;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.learning_rate_start=0.005;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=7;
ParamsS.train_time=60;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4';
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=5;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4';
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');

ParamsS=ParamsSDefaults;
ParamsS.MapSize=7;
ParamsS.train_time=60;
ParamsS.dataset='dataFaceP4';
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');
%%
ps aux | grep RunTFForMatlab
%%
GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);
GOP1 = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,ones(Sz2));

Q1=GOP1'*(A(1,:));
Q2=GOP1'*(A(1,:).');

ShowAbsAngle(ImgMskedPhased(:,:,1))

Partition by TSB
%%
PMain='/home/a/TF/srez/';
DMain=dir(PMain);
DMain=DMain([DMain.isdir]);
DMain=DMain(3:end);
DMainN={DMain.name};
for i=1:numel(DMainN)
    system(['chmod -R 777 ' PMain DMainN{i}]);
end