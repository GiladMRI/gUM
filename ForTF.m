Data=rand(5000,28,28,2);
Labels=rand(5000,28,28,2);
save('/home/a/tensorflowtest/MyData_Img.mat','Data','Labels')

%%
DatasetP='/home/a/TF/srez/part of dataset/';
OutP='/home/a/TF/srezx/dataset';
D=dir([DatasetP '*.jpg']);
DFNs={D.name};
Ni=18000;
for i=1:Ni
    if(mod(i,1000)==1), disp(i), end
    TFData(:,:,:,i)=imread([DatasetP DFNs{i}]);
end
% 218   *178     *3
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

TFDatacN=double(TFDatac)/255;
%%
TFDatacM=squeeze(TFDatacN);
%%
save('/home/a/TF/Imgs.mat','TFDataP')
%%
Data=permute(TFDatacN,[4 1 2 3]);
Labels=repmat(mean(Data,4),[1 1 1 1]);
save('/home/a/TF/Imgs.mat','Data','Labels')
%%
M=ones(1,64,64,3);
M(:,1:32,1:32,:)=0;
save('/home/a/TF/M.mat','M')
%%
N=size(TFDataP,2);
NN=N*N;
Idxs1=((-N/2):((N/2)-1))+0.5;

Trajx=repmat(Idxs1,[1 N]);
Trajy=kron(Idxs1,ones([1 N]));
Traj=[Trajx.' Trajy.'];
Coeffs=gFTCoeffs2D([N N],Traj);
CoeffsM=reshape(Coeffs,[N*N size(Traj,1)]);
a=squeeze(mean(double(TFDataP(3,:,:,:)),4));

Data=gsum(RepDotMult(a,Coeffs),1:2);

R=squeeze(Data).'*(CoeffsM');
%%
DataP=permute(Data,[1 3 2]);
DataR=cat(2,real(DataP),imag(DataP));

AdjM=CoeffsM.';
% AdjMr=real(AdjM);
% AdjMi=imag(AdjM);
% 
% AdjMR=[    AdjMr -AdjMi ;  AdjMi AdjMr];

AdjMR=MatComplexToRealAndImag(AdjM);

RR=DataR*AdjMR;

RC=RR(1:NN) + 1i*RR(NN+1:end);

RCI=reshape(RC,[N N]);
ShowAbsAngle(RCI)
%%
AdjMRx=repmat(permute(AdjMR,[3 4 1 2]),[16 1 1 1]);
save('/home/a/TF/AdjMRx.mat','AdjMRx')
%%
save('/home/a/TF/AdjMR.mat','AdjMR')
%%
B=rand(1,NN)>0.5;
BB=[B B];
nChosen=sum(B);
AdjMRB=MatComplexToRealAndImag(AdjM(B,:));
save('/home/a/TF/AdjMRB.mat','AdjMRB','nChosen');
%%
for i=1:Ni
    disp(i)
    BWAll=permute(mean(double(TFDataP(i,:,:,:)),4),[2 3 4 1]);
    Data=permute(gsum(RepDotMult(BWAll,Coeffs),1:2),[4 1 3 2]);
    DataR(i,:,:,:)=cat(3,real(Data),imag(Data));
end
%%
Data=permute(double(TFDatac),[4 1 2 3])/255;
Labels=repmat(mean(Data,4),[1 1 1 1]);

Data=permute(DataR(:,:,:,BB),[1 2 4 3])/255;
save('/home/a/TF/ImgsX.mat','Data','Labels')
%%
DataP=permute(Data,[1 3 2]);
%%
TFDataPF=double(TFDataP)/255;
Labels=grmss(TFDataPF,4);
save('/home/a/TF/ImgsMC.mat','Data','Labels')
%% pseudo EPI prob
Labels=grmss(TFDataPF,4);
PLabels=permute(Labels,[2 3 1]);
F=fft2cg(PLabels);
F2=repmat(F,[1 1 1 2]);
F2(2:2:end,:,:,1)=0;
F2(1:2:end,:,:,2)=0;

L=linspace(-pi,pi,size(F2,2));
for i=1:size(F2,3)
    Bias=(rand*2*pi-pi)/3;
    Scale=(rand*10-5)/3;
    F2(:,:,i,2)=RepDotMult(F2(:,:,i,2),exp(1i*(Bias+L*Scale)));
end
Data=permute(ifft2c(F2),[3 1 2 4]);
Data=cat(4,real(Data),imag(Data));
save('/home/a/TF/ImgsMC.mat','Data','Labels')
%%
F3=ifft2c(sum(F2,4));
%%
TFDatacM=grmss(TFDatacN,3);

Msk=imfillholesBySlices(grmss(Sens,3)>0.3);

GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);
for j=1:nTS
    GOPC{j}=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTraj2,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,TSC(:,:,j));
end
GOP=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTraj2,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,ones(Sz2));

SWarning=warning('off','gpuNUFFT:adj:sens');

TSBF=TSB(IB_TimeInMs2,:);
%%
% 
%CurIR=NaN(Sz2(1),Sz2(2),nTS*ncc,Ni);
CurDataM=NaN(ncc,nTrajP2,Ni);
for i=1:Ni
    disp(i);
    CurI=TFDatacM(:,:,i).*Msk;
    CurData=GOP_MC*CurI;
    CurDataM(:,:,i)=CurData;
end
%%
CurDataMR=repmat(CurDataM,[2 1 1 1]);
CurDataMR(1:2:end,:,:,:)=real(CurDataM);
CurDataMR(2:2:end,:,:,:)=imag(CurDataM);

Data=permute(CurDataMR,[3 4 2 1]);

Labels=permute(RepDotMult(TFDatacM,Msk),[3 1 2]);
save('/home/a/TF/ImgsMCD.mat','Data','Labels')
%%
CurDataMP=permute(CurDataM,[2 1 3]);
CurDataMPF=reshape(CurDataMP,prod(gsize(CurDataMP,1:2)),size(CurDataMP,3));

CurDataMPFP=permute(CurDataMPF,[2 1]);

CurDataMPFPR=[real(CurDataMPFP) imag(CurDataMPFP)];

Data=CurDataMPFPR;
Labels=permute(RepDotMult(TFDatacM,Msk),[3 1 2]);
save('/home/a/TF/ImgsMCD.mat','Data','Labels')
%%
CurDataMR=repmat(CurDataMP,[2 1 1 1]);
CurDataMR(1:2:end,:,:,:)=real(CurDataMP);
CurDataMR(2:2:end,:,:,:)=imag(CurDataMP);

%%
i=325;
CurI=TFDatacM(:,:,i).*Msk;
CurData=GOP_MC*CurI;

Coeffs=gFTCoeffs2D(Sz2,BARTTraj2.');
CoeffsWithB0=Coeffs.*AllB0C(:,:,IB_TimeInMs2);
CoeffsWithB0Sens=RepDotMult(CoeffsWithB0,permute(Sens,[1 2 4 3]));

CurData2=squeeze(gsum(RepDotMult(CurI,CoeffsWithB0Sens),1:2)).';

CoeffsWithB0SensM=reshape(CoeffsWithB0Sens,prod(Sz2),nTrajP2*ncc);

% piCoeffsWithB0SensM=pinv(CoeffsWithB0SensM);
%%
N=64;
Dx4D=zeros([N N N-1 N]);
Dy4D=zeros([N N N N-1]);
for i=1:N-1
    for j=1:N
        Dx4D(i,j,i,j)=-1;
        Dx4D(i+1,j,i,j)=1;
        Dy4D(j,i,j,i)=-1;
        Dy4D(j,i+1,j,i)=1;
    end
end
DMx=reshape(Dx4D,N*N,(N-1)*N);
DMy=reshape(Dy4D,N*N,(N-1)*N);
DM=[DMx DMy];
%%
CoeffsWithB0SensMWithD=[CoeffsWithB0SensM DM*0.03];

tic
piCoeffsWithB0SensMWithD=pinv(CoeffsWithB0SensMWithD);
toc

piMDx=piCoeffsWithB0SensMWithD(1:(ncc*nTrajP2),:);
%%
figure;
subplot(2,2,1);
gmontage(CurI,[0 1]);
% subplot(2,2,3);
% X=Row(CurData.')*piCoeffsWithB0SensM;
% XX=reshape(X,Sz2);
% gmontage(XX,[0 1])
subplot(2,2,4);
YD=Row(CurData.')*piMDx;
YYD=reshape(YD,Sz2);
gmontage(YYD,[0 1])
%%
GPhase=NaN(N,N,Ni);
for i=1:Ni
%     GPhase(:,:,i)=GenerateRandomBrokenPhase(N);
    GPhase(:,:,i)=GenerateRandomBrokenPhase(N,5,1,4);
    if(mod(i,1000)==0), disp(i),end
end
%%
ImgMskedPhased=NaN(N,N,Ni);
DataMskedPhased=NaN(ncc,3024,Ni);
for i=1:Ni
    ImgMskedPhased(:,:,i)=TFDatacM(:,:,i).*Msk.*GPhase(:,:,i);
    DataMskedPhased(:,:,i)=GOP_MC*ImgMskedPhased(:,:,i);
    if(mod(i,100)==0), disp(i),end
end
%%
dataFoldName='dataFaceP4';
mkdir(['/home/a/TF/srez/' dataFoldName '/']);
system(['chmod -R 777 ' '/home/a/TF/srez/' dataFoldName '/']);

[status,msg,msgID] = fileattrib('/home/a/TF/srez/dataFaceP4/','+w','a');

ChunkSize=100;
ChunkStartI=1:ChunkSize:size(DataMskedPhased,3);
ChunkEndI=min(ChunkStartI+ChunkSize-1,size(DataMskedPhased,3));
for k=1:numel(ChunkStartI)
    CurIs=ChunkStartI(k):ChunkEndI(k);
    CurChunkSize=numel(CurIs);
    
    CurDataMP=permute(DataMskedPhased(:,:,CurIs),[2 1 3]);
    CurDataMPF=reshape(CurDataMP,prod(gsize(CurDataMP,1:2)),size(CurDataMP,3));
    
    CurDataMPFP=permute(CurDataMPF,[2 1]);
    
    CurDataMPFPR=[real(CurDataMPFP) imag(CurDataMPFP)];
    
    LabelMskedPhasedP=permute(ImgMskedPhased(:,:,CurIs),[3 1 2]);
    LabelsD=cat(4,real(LabelMskedPhasedP),imag(LabelMskedPhasedP));
    
    Data=single(CurDataMPFPR(:,:));
    Labels=single(LabelsD(:,:,:,:));
    FNs=strcat(num2str(CurIs.','%05d'),'asd');
    save('/home/a/TF/CurChunk.mat','Data','Labels','FNs')
    
    system(['~/b.sh /home/a/TF/Mat2TFRec.py /home/a/TF/CurChunk.mat /home/a/TF/srez/' dataFoldName '/ '])
end
%%
CurDataMP=permute(DataMskedPhased,[2 1 3]);
CurDataMPF=reshape(CurDataMP,prod(gsize(CurDataMP,1:2)),size(CurDataMP,3));

CurDataMPFP=permute(CurDataMPF,[2 1]);

CurDataMPFPR=[real(CurDataMPFP) imag(CurDataMPFP)];


ImgMskedPhasedP=permute(ImgMskedPhased,[3 1 2]);
LabelsD=cat(4,real(ImgMskedPhasedP),imag(ImgMskedPhasedP));

Data=single(CurDataMPFPR(1:10000,:));
Labels=single(LabelsD(1:10000,:,:,:));
save('/home/a/TF/ImgsMCD3x4.mat','Data','Labels')

QQ=dir('/home/a/TF/ImgsMCD3x.mat')
QQ.bytes/1e6
%%

%%
% RealDataMP=repmat(permute(DataP2,[2 4 1 3]),[1 1 16]);
RealDataMP=repmat(DataP2,[1 1 16]);
RealDataMPF=reshape(RealDataMP,prod(gsize(RealDataMP,1:2)),size(RealDataMP,3));

RealDataMPFP=permute(RealDataMPF,[2 1]);

RealDataMPFPR=[real(RealDataMPFP) imag(RealDataMPFP)];
Data=single(RealDataMPFPR*30);
save('/home/a/TF/srez/RealData/a.mat','Data')
%%
i=13;
CurI=TFDatacM(:,:,i).*Msk.*GPhase(:,:,i);
CurData=GOP_MC*CurI;

YD=Row(CurData.')*piMDx;
YYD=reshape(YD,Sz2);

figure;
subplot(2,2,1);
gmontage(abs(CurI),[0 1]);
subplot(2,2,2);
gmontage(angle(CurI),[-pi pi]);
subplot(2,2,3);
gmontage(abs(YYD),[0 1]);
subplot(2,2,4);
gmontage(angle(YYD),[-pi pi]);

%%
for i=1:Ni
    disp(i);
    CurI=TFDatacM(:,:,i).*Msk;
    CurData=GOP_MC*CurI;
    for c=1:ncc
        CurDataB=RepDotMult(CurData(c,:),TSBF');
        CurIR(:,:,(c-1)*nTS+(1:nTS),i)=GOP'*(CurDataB(:,:).');
%         for j=1:nTS
%             dd(:,:,j)=GOP'*(CurDataB(j,:)');
% %             CurIR(:,:,(j-1)*ncc+c,i)=GOPC{j}'*(.');
%         end
    end
end

CurIRR=repmat(CurIR,[1 1 2 1]);
CurIRR(:,:,1:2:end,:)=real(CurIR);
CurIRR(:,:,2:2:end,:)=imag(CurIR);
Data=permute(CurIRR,[4 1 2 3]);
Labels=permute(RepDotMult(TFDatacM,Msk),[3 1 2]);
save('/home/a/TF/ImgsMC.mat','Data','Labels')
%%
% x = randn(Sz2) + 1j*randn(Sz2);
% y = randn([size(Sens,3) nTrajP2]) + 1j*randn([size(Sens,3) nTrajP2]);
% Ax = GOP_MC*x;
% Aty = GOP_MC'*y;
% Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
% %%
% Sens_TS=reshape(RepDotMult(TSC,permute(Sens,[1 2 4 3])),Sz2(1),Sz2(2),[]);
% 
% GOP=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTraj2,Sz2)'/(2*pi),ones(nTrajP2,1),osf,wg,sw,Sz2,Sens_TS);
% GOP2=gpuNUFFT(BART2Fes_NUFT_Idxs(BARTTraj,Sz)'/(2*pi),ones(nTrajP,1),osf,wg,sw,Sz,TSC);
% 
% Ax = GOP*x;
% 
% I=

Len=5;
SessionName='WithpMatAcc3x';
save('/home/a/TF/MatlabParams.mat','Len','SessionName');

SessionName='OnlyManifoldAfterpMat';
save('/home/a/TF/MatlabParams.mat','Len','SessionName');

piMDy=reshape(piMDx,size(piMDx,1),N,N);
piMDz=reshape(permute(piMDy,[1 3 2]),size(piMDx,1),N*N);
piMDR=MatComplexToRealAndImag(piMDz);
% piMDR=piMDR(:,1:4096);
% piMDR(:,MIdxs)=0;
piMDR(:,[MIdxs MIdxs+N*N])=0;
piMDRA=piMDR;
piMDR=piMDRA;
piMDR(:,1:2:end)=piMDRA(:,1:4096);
piMDR(:,2:2:end)=piMDRA(:,4096+(1:4096));
save('/home/a/TF/Consts.mat','piMDR');
%%
TFRes=load(['/home/a/TF/srez/' SessionName '_train/TrainSummary.mat']);
F=find(TFRes.G_LossV(5:end)==0,1);
figure;plot(TFRes.G_LossV(1:F))
%%
D=dir(['/home/a/TF/srez/' SessionName '_train/TrainSummary_*.mat']);
TFRes=load(['/home/a/TF/srez/' SessionName '_train/' D(end).name]);
%%
TFRes=load(['/home/a/TF/srez/' SessionName '_train/TrainSummary_020000.mat']);

i=34;
A=Data(i,:)*TFRes.gene_GEN_L003_weight_0+TFRes.gene_GEN_L003_bias_0;
%%
for i=1:Ni
    A=Data(i,:)*TFRes.gene_GEN_L003_weight_0+TFRes.gene_GEN_L003_bias_0;
    Recd1(:,:,i)=reshape(A,[64 64]).';
end

Data=permute(abs(Recd1),[3 1 2]);
Labels=permute(RepDotMult(TFDatacM,Msk),[3 1 2]);
save('/home/a/TF/ImgsMC1.mat','Data','Labels')

for i=1:Ni
    CurData=squeeze(Recd1(:,:,i));
    save(['/home/a/TF/srez/DataAfterpMat/' num2str(i) '.mat'],'CurData')
end

%%
for i=1:Ni
    A=load(['/home/a/TF/srez/OnlyManifoldAfterpMat_OutMat/' num2str(i) '_out.mat']);
    Recd2(:,:,i)=A.X;
end
%%
i=142;
AA=TFDatacM(:,:,i).*Msk;
fgmontage(cat(3,Recd1(:,:,i),Recd2(:,:,i),AA),'Size',[1 3]);

fgmontage(cat(3,Recd1(:,:,i)-AA,Recd2(:,:,i)-AA),'Size',[1 2]);
%%
MIdxs=find(~Row(Msk.'));
piMDRt=piMDR;
piMDRt(:,MIdxs)=0;
aa=CurDataMPFPR(200,:);
bb=aa*piMDRt;
fgmontage(cat(3,reshape(bb,[64 64]).',squeeze(Labels(200,:,:))),[0 0.01])

%%
ParamsSDefaults=struct('DataH',36288,'DataW',1,'channelsIn',1,'LabelsH',64,'LabelsW',64,'channelsOut',2,...
  'dataset','dataKnee','SessionNameBase','Acc3','learning_rate_start',0.002,'learning_rate_half_life',5000,...
  'MapSize',3,'train_time',120);

ParamsS=ParamsSDefaults;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
%% Run several
ParamsS=ParamsSDefaults;
ParamsS.MapSize=3;
ParamsS.train_time=60;
gStruct2txt(ParamsS,'/home/a/TF/Params.txt')
system('/home/a/RunTFForMatlab.sh');
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