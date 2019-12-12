%% use Native bart
filename='/media/a/DATA1/2018_01_25/meas_MID147_spiral_vas_T1_OUT_minTE_FatSat_FID291.dat';
% filename='/media/a/DATA/2018_01_25/meas_MID139_spiral_vas_T1_IN_TE15_FatSat_FID283.dat';
dt=2.5E-6;
dg=1E-5; % 10us grad raster
nit=5;
lambda=.0;
dx=-25;
dy=0;
read_raw=1;

%% Read raw
if read_raw
    sTwix = mapVBVD(filename);
   % save sTwix sTwix
else
    load sTwix;
end

% set parameters
NCol=sTwix.image.NCol;
NSli=sTwix.image.NSli;
NPar=sTwix.image.NPar;
NLin=sTwix.image.NLin;
NCha=sTwix.image.NCha;
NEco=sTwix.image.NEco;
NSet=sTwix.image.NSet;
NSeg=sTwix.image.NSeg;
NRep=sTwix.image.NRep;
NAve=sTwix.image.NAve;


Amp=sTwix.image.freeParam(2)/10;
Slew=sTwix.image.freeParam(1)*100;
FOV=sTwix.image.freeParam(3)/10;
NRead=NSeg*NCol;

% make raw double
data=double(sTwix.image());

%perform averaging (DIM=6)
data = mean(data, 6);


% k-space trajectory
[kx,ky,kxr,kyr,N] = SiemensSpiral(FOV,NRead,NLin,Amp,Slew,dt,0,dg);
if(0) % for Spiral In
    kx=kxr;ky=kyr;
end
% [kx,ky,~,~,N] = SiemensSpiral(FOV,NReadx,NLin,Amp,Slew,dt,0,dg);

N=double(floor(N));
kxFull=double(kx*FOV/N)*2*pi;
kyFull=double(ky*FOV/N)*2*pi;

%%
ShotsIdxs=kron(1:NLin,ones(1,NRead));
WhichShots=[1:1:NLin];
nWhichShots=numel(WhichShots);
B=ismember(ShotsIdxs,WhichShots);
kx=kxFull(B);
ky=kyFull(B);
modx=double(exp(1i*(dx*kx+dy*ky))');

A = nuFTOperator([kx ky],[N N]);



NRep = 1;
%%
kx1=kx;
ky1=ky;
N1=N;
%%
N=64;

%%
Sz1=[N1 N1];
BARTTraj=FES2BART_NUFT_Idxs([kx1 ky1],Sz1);

%%
nTrajP=size(BARTTraj,2);
% Sz2=[72 72];
Sz2=[N N];
% Sz2=[128 128];
% Sz2=[256 256];
TrajR=sqrt(BARTTraj(1,:).^2+BARTTraj(2,:).^2);
FirstOutPoint=find(TrajR>(Sz2(1)/2),1);
BTraj=repmat([true(1,FirstOutPoint) false(1,NRead-FirstOutPoint)], [1,NLin]);

BTrajAllShots=BTraj;
BARTTraj2AllShots=BARTTraj(:,BTraj);
BTraj(ismember(mod(1:nTrajP,4),[0 2]))=false;
ShotB=ismember(ShotsIdxs,1:3:NLin);
BTraj=BTraj & ShotB;

BARTTraj2=BARTTraj(:,BTraj);

nTrajP2=size(BARTTraj2,2);
Acc=prod(Sz2)/nTrajP2;

figure;plot(BARTTraj2(1,:),BARTTraj2(2,:),'.');axis([-Sz2(1) Sz2(1) -Sz2(2) Sz2(2)]/2);
title(num2str([Sz2 Acc]))

TimeInMs=repmat(linspace(0,20,NRead),[1 nWhichShots]);
%% Arrange multichannel data
SliI=1;
clear projC
for c=1:NCha
    temp=squeeze(sum(data(:,c,WhichShots,1,SliI,1,1,1,:,1,:,1,1,1,1,1),9));
    projC(:,:,:,c)=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;
end
projC=permute(projC,[2 1 3 4]);
%%
BadChannelsI=[6 11 29];
GoodChannelsI=setdiff(1:NCha,BadChannelsI);
nGoodChannels=numel(GoodChannelsI);
DataSmallGoodChannelsAllShots=squeeze(conj(projC(:,BTrajAllShots,:,GoodChannelsI)));
%% All channels
clear recoC
for c=1:nGoodChannels
    disp(c);
    recoCs(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.001 -t ',BARTTraj2AllShots, DataSmallGoodChannelsAllShots(:,c).', ones(Sz2));
end

SensB=RunESPIRiTForSensMaps(recoCs,42);
ShowAbsAngle(SensB)

recoCx=CalcSENSE1f(recoCs,SensB);
ShowAbsAngle(recoCx)

%%
DataSmallGoodChannels=squeeze(conj(projC(:,BTraj,:,GoodChannelsI)));

[U,S,V] = svd(DataSmallGoodChannels,'econ');
sccmtx=V;
ncc=8;
% sccmtx=eye(ncc);
% SensB=ApplySCCBySlices(SensBA(:,:,GoodChannelsI),sccmtx,ncc);
DataSmallCCA=DataSmallGoodChannelsAllShots*  sccmtx(:,1:ncc);
DataSmallCC=DataSmallGoodChannels*  sccmtx(:,1:ncc);
SensCC=permute(MultMatTensor(sccmtx(:,1:ncc).',permute(SensB,[3 1 2])),[2 3 1]);
%% Small
clear recoCsCCA recoCsCC
for c=1:ncc
    disp(c);
    recoCsCCA(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.001 -t ',BARTTraj2AllShots, DataSmallCCA(:,c).', ones(Sz2));
    recoCsCC(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.001 -t ',BARTTraj2, DataSmallCC(:,c).', ones(Sz2));
end

recoCsCCAx=CalcSENSE1f(recoCsCCA,SensCC);
recoCsCCx=CalcSENSE1f(recoCsCC,SensCC);

ShowAbsAngle(recoCsCCAx)
ShowAbsAngle(recoCsCCx)
%%
SensP=permute(SensCC,[1 2 4 3]);

PICSprm='-m -r:0.0001 -R W:7:0:0.0001 -t';
PICSprm='-m -r:0.0001 -R T:7:0:0.0001 -t';
% recoX = bart('pics -r:0.0001 -R T:7:0:0.001 -t ',TrajectoryP, DataP, SensP);
recoXA = bart(['pics ' PICSprm],BARTTraj2AllShots, permute(DataSmallCCA,[3 1 4 2]), SensP);
recoX = bart(['pics ' PICSprm],BARTTraj2, permute(DataSmallCC,[3 1 4 2]), SensP);
fgmontage(recoX);title(PICSprm);
axis equal
%% Now NUFFT and TS
load('CurStatus11.mat','B0RealExs')
%% B0 map
B0M2=imresizeBySlices( B0RealExs,Sz2);

TimeInMs2=TimeInMs(BTraj)*1;
%% All B0 effects across time
[U_TimeInMs2, IA_TimeInMs2, IB_TimeInMs2]=unique(TimeInMs2);
nU_TimeInMs2=numel(U_TimeInMs2);

AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
%% Fessler time segmentation
% nTS=7;
clear ErrTS
TS_Thresh=1e-6;
for nTS=3:15
    disp(nTS)
    FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    
    TSC2=reshape(TSC,prod(Sz2),nTS);
    tic
    TSB=(E.')/(TSC2.');
    toc
    ErrTS(nTS)=grmss(E-TSC2*(TSB.'));
    if(ErrTS(nTS)<TS_Thresh)
        disp(['Stopped at #TS=' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
        break;
    end
end
figure;plot(log10(ErrTS),'-*')
%% GPU TS
Sens=squeeze(SensP);

% w=ones(nTrajP2,1);
osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);

x = randn(Sz2) + 1j*randn(Sz2);
y = randn([size(Sens,3) nTrajP2]) + 1j*randn([size(Sens,3) nTrajP2]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
%%
DataP=DataSmallCC.';

AOdd = GOP_MC;

TVW=0.1;

param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',1,'TV',TVOP_MSlice,'xfmWeight',0),init);

param.data =     DataP;

param.pNorm=1;
param.TVWeight=0.0001;

nfnlCgIters=15;
RunFnlViewAmp=2;
res=single(zeros(Sz2));
res=zeros(Sz2);

FigH=4000;
figure(FigH);close(FigH);

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
for n=1:nfnlCgIters
%     disp([Slis WhichRep n]);
    res = fnlCg(res,param);
    im_res = param.XFM'*res;
%     figure(FigH), gmontage(abs(gflip(im_res,1))), drawnow;% title(qq)
    if(param.ShowFig)
        figure(FigH); subplot(1,2,1);gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,2,2);gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
%         xlabel([Slis WhichRep n]);
    end
%     t=toc;
    disp(['Iter #' num2str(n) ' ' datestr(now)]);
end
% t=toc;
%%
figure;
subplot(1,3,1);
gmontage(grmss(A.x,3),[0 0.7]);title('TS-NN');axis equal
subplot(1,3,2);
gmontage(recoX,[0 0.5]);title(PICSprm);title('BART. No B0');axis equal
subplot(1,3,3);
gmontage(im_res,[0 0.03]);title('Iterative TS-gpuNUFFT with TV');axis equal
%%
fgmontage(cat(3,grmss(A.x,3)/0.7,abs(recoX/0.5),abs(im_res/0.03)),'Size',[1 3])
title('TS-NN                                                BART. No B0                       Iterative TS-gpuNUFFT with TV');
%% Now TF