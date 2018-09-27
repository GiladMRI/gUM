addpath(genpath('DLMRI_v6'));
%%
I=rgb2gray(imread('DLMRI_v6/BMR2.JPG'));
FI=fft2cg(I);
Q=load('DLMRI_v6/SamplingMasksDLMRI_TMI/Figure1/Q1.mat');Q=Q.Q1;Q=gfftshift(Q,1:2);

n=49;
K2=n+1;
num=20;
DLMRIparams=struct('num',num,'n',n,'K2',K2,'N',200*K2,'T0',round((0.2)*n),'Lambda',140,'KSVDopt',2,'thr',(0.023)*[2 2 2 2 1.4*ones(1,num-4)],'numiterateKSVD',20,'r',1);


[Iout1,param1] = DLMRIDemoComplex(I,Q,0,0,DLMRIparams);
%%
I=rgb2gray(imread('DLMRI_v6/BMR2.JPG'));
I=imresize(I,[256 256]);
FI=fft2cg(I);
Q=load('DLMRI_v6/SamplingMasksDLMRI_TMI/Figure1/Q1.mat');
Q=Q.Q1;
Q=gfftshift(Q,1:2);
Q=crop(Q,[256 256]);

n=25;
K2=n+1;
num=20;
DLMRIparams=struct('num',num,'n',n,'K2',K2,'N',200*K2,'T0',round((0.2)*n),'Lambda',140,'KSVDopt',2,'thr',(0.023)*[2 2 2 2 1.4*ones(1,num-4)],'numiterateKSVD',20,'r',1);


[Iout1,param1] = DLMRIDemoComplex(I,Q,0,0,DLMRIparams);
DD=reshape(param1.Dictionary,[5 5 26]);
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03')
% RecG3=bart(['pics -S -m -R T:3:0:' num2str(0.0003) ' -W /tmp/RR -t'],BARTTraj, DataP, SensP);
tic
RecBase=bart(['pics -S -m -R T:3:0:' num2str(0.0003) ' -p '],Q, FI.*Q, FI*0+1);
t1=toc;
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03ConvSparse')
% RecG3=bart(['pics -S -m -R T:3:0:' num2str(0.0003) ' -W /tmp/RR -t'],BARTTraj, DataP, SensP);
nDict=K2;
IndicatorCollapse=diag(ones(1,nDict),1);
IndicatorCollapse=IndicatorCollapse(1:end-1,1:end-0);
% IndicatorCollapse(1)=0;
IndicatorCollapse=permute(IndicatorCollapse,[8 7 6 1 5 4 3 2]);
% IndicatorCollapse=repmat(IndicatorCollapse,[256 256]);
writecfl('/tmp/ICollapse',IndicatorCollapse);

OneS=eye(nDict+1);
OneS=permute(OneS,[9 8 7 6 5 4 3 2 1]);
writecfl('/tmp/CRTensorOneS',OneS);
OneSB=eye(nDict+1);
OneSB(1,1)=0;
OneSB=permute(OneSB,[9 8 7 6 5 4 3 2 1]);
writecfl('/tmp/CRTensorOneSB',OneSB);

OneSC=zeros(nDict+1);
OneSC(1,1)=1;
OneSC=permute(OneSC,[9 8 7 6 5 4 3 2 1]);
writecfl('/tmp/CRTensorOneSC',OneSC);

CRTensor=DD;
CRTensor(256,256,1)=0;
CRTensor=gfftshift(CRTensor,1:2);
CRTensor=ifft2cg(CRTensor);
% CRTensor=gflip(CRTensor,1:2);
CRTensor=cat(3,-ones(256,256),CRTensor);
% CRTensor=permute(CRTensor,[1 2 8 7 6 5 4 3]);
CRTensor=permute(CRTensor,[1 2 9 8 7 6 5 4 3]);
writecfl('/tmp/CRTensor',CRTensor);

SensP=FI*0+1;
SensP(1,1,1,1,1,1,1,nDict+1)=0;
DataP=FI.*Q;
DataP(1,1,1,1,1,1,1,nDict+1)=0;
Pattern=Q;
% Pattern(1,1,1,1,1,1,1,nDict+1)=0;
% SensP=SensP;

Pattern=repmat(Q,[1 1 1 1 1 1 1 nDict+1]);

WarmStart=repmat(Iout1,[1 1 1 1 1 1 1 nDict+1]);
writecfl('/tmp/WarmStart',WarmStart);

tic
CSRatio=1e-3;
CSBase=1e-2;
% RecI=bart(['pics -S -m -R T:3:0:' num2str(0.0003) ' -p '],Pattern, DataP, SensP);
% RecI=bart(['pics -S -m -R I:3:0:' num2str(0.0003) ' -p '],Pattern, DataP, SensP);
% RecC=bart(['pics -S -m -R C:3:0:' num2str(0.0003) ' -p '],Pattern, DataP, SensP);
% RecIC=bart(['pics -S -m -R C:3:0:' num2str(CSBase) ' -R I:3:0:' num2str(CSBase*CSRatio) ' -p '],Pattern, DataP, SensP);
RecIC=bart(['pics -S -m -R C:3:8:' num2str(CSBase) ' -R CI:3:0:' num2str(CSBase) ' -R I:3:0:' num2str(CSBase*CSRatio) ' -p '],Pattern, DataP, SensP);
% RecIC=bart(['pics -S -m -W /tmp/WarmStart -R C:3:0:' num2str(3) ' -R I:3:0:' num2str(0.00003) ' -p '],Pattern, DataP, SensP);
t1=toc;
% RecIS=squeeze(RecI);
% RecCS=squeeze(RecC);
RecICS=squeeze(RecIS);
fgmontage(RecIC(:,:,1));colorbar;
xlabel([CSBase CSRatio]);
fgmontage(RecIC(:,:,2:10));colorbar;
%%
FRecC=fft2cg(RecC);