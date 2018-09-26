clear

run ( '/home/a/irt/setup.m');


BaseP='/home/a/gUM/4GL/';

addpath([BaseP 'Spiral_recon_T1/'])
addpath(genpath([BaseP 'Spiral_recon_T1/raw_header']))
addpath([BaseP 'Spiral_recon_T1/misc'])
addpath([BaseP 'Spiral_recon_T1/io'])

addpath(genpath('/home/a/SPENOnline'))

addpath(genpath('/home/a/gUM/EPFLSpiral'))

addpath(genpath('/home/a/gpuNUFFT-master'))
%%
% addpath(genpath('/home/a/gpuNUFFT-master/CUDA'))

%% Input parameters
filename='/media/a/DATA1/2018_01_25/meas_MID147_spiral_vas_T1_OUT_minTE_FatSat_FID291.dat';
% filename='/media/a/DATA1/2018_01_25/meas_MID139_spiral_vas_T1_IN_TE15_FatSat_FID283.dat';
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

%% 2D spiral recon
imc=zeros(N,N,NCha,NSli,NPar,NEco,NSet,NRep);
for r=1:NRep
    for n=1:NSet
        for e=1:NEco
            for p=1:NPar
                for s=1 %:NSli
                    for c=1:NCha
                        temp=squeeze(data(:,c,WhichShots,p,s,1,1,e,r,n,:,1,1,1,1,1));
                        proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;
%                         temp=squeeze(data(:,c,WhichShots,p,s,1,1,e,r,n,1,1,1,1,1,1));
%                         proj=reshape(permute(temp,[1 3 2]),[1 NReadx*nWhichShots]).*mod;
                        imc(:,:,c,s,p,e,n,r) = regularizedReconstruction(A,proj',@L2Norm,lambda,'maxit', nit);
                        progresscounter(NEco*NSli*NCha*NRep*NPar*NSet,'NuFFalready: ');
                    end
                end
            end
        end
    end
end

if(NPar>1)
    imc=fff(imc,5);
    im=sqrt(squeeze(sum(abs(imc).*abs(imc),3)));%.*cylind(N,NPar,N/2);
else
    im=sqrt(squeeze(sum(abs(imc).*abs(imc),3))); %.*cylind(N,NSli,N/2);  
end
fgmontage(im(:,:,1));
%%
Z=CalcSENSE1f(imc(:,:,:,1),SensB(1:end-1,1:end-1,:));
% show3d(im);
%save im im;
%% Again, organized
c=3;
s=1;
temp=squeeze(sum(data(:,c,WhichShots,1,s,1,1,1,:,1,:,1,1,1,1,1),9));
proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;

X = regularizedReconstruction(A,proj',@L2Norm,lambda,'maxit', nit);
%% Arrange multichannel data
clear projC
for c=1:NCha
    temp=squeeze(sum(data(:,c,WhichShots,1,s,1,1,1,:,1,:,1,1,1,1,1),9));
    projC(:,:,:,c)=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;
end
projC=permute(projC,[2 1 3 4]);
%%
kx1=kx;
ky1=ky;
N1=N;
%%
AA = nuFTOperator([kx1 ky1],[N1 N1]);

OneIter=AA'*proj';
fgmontage(OneIter)
%%
Sz1=[N1 N1];
BARTTraj=FES2BART_NUFT_Idxs([kx1 ky1],Sz1);

if(mod(Sz1(1),2)==1)
    BARTTraj(1,:)=BARTTraj(1,:)-0.5;
end
if(mod(Sz1(2),2)==1)
    BARTTraj(2,:)=BARTTraj(2,:)-0.5;
end

Out=bart('nufft -a',BARTTraj,conj(proj));
% 
% if(mod(Sz1(1),2)==1)
%     Out=Out.*exp(1i* 2*( (BARTTraj(1,:)+0.5)*pi/Sz1(1)));
% end
% if(mod(Sz1(2),2)==1)
%     Out=Out.*exp(1i* 2*( (BARTTraj(2,:)+0.5)*pi/Sz1(2)));
% end

% NUbyBA3=gBARTnuFT(BARTTraj,proj,[N1 N1]);
%%
Outx=Out;%gflip(Out,1:2);
figure;
subplot(3,2,1);
imagesc(abs(OneIter),[0 4e-3]);
subplot(3,2,2);
imagesc(angle(OneIter),[-pi pi]);
subplot(3,2,3);
imagesc(abs(Outx),[0 4e-3]);
subplot(3,2,4);
imagesc(angle(Outx),[-pi pi]);
D=Outx-OneIter;
subplot(3,2,5);
imagesc(abs(D));
subplot(3,2,6);
imagesc(angle(D),[-pi pi]);
%%

%% use Native bart
setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensorClean')

Sz1=ceil(Sz1/2)*2;
SensP=ones(Sz1);
% SensP=gflip(conj(Sens(:,:,3)),1:2);

nChannels=size(I,3);

DataP=conj(proj);
TrajectoryP=BARTTraj;
reco = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',TrajectoryP, DataP, SensP);
fgmontage(reco)
% reco = bart('pics -S -r:0.01 -RT:1024:0:0.0001 -p Weights  -t ',Trajectory, Data, SensP);
%%
clear recoC
for c=1:NCha
    disp(c);
    DataP=conj(projC(:,:,:,c));
    recoC(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',TrajectoryP, DataP, SensP);
end
%% Small
clear recoC
for c=1:NCha
    disp(c);
    DataPC=conj(projC(:,BTraj,:,c));
    recoCs(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',BARTTraj2, DataPC, SensP);
end

%%
SensB=RunESPIRiTForSensMaps(recoC,30);

recoCx=CalcSENSE1f(recoC,SensB);
%%
projCA=projC;
SensBA=SensB;
%%
BadChannelsI=[6 11 29];
GoodChannelsI=setdiff(1:NCha,BadChannelsI);

calibDATA = reshape(projCA(:,:,:,GoodChannelsI),[size(projCA,2),numel(GoodChannelsI)]);
[U,S,V] = svd(calibDATA,'econ');
sccmtx=V;
ncc=6;
% sccmtx=eye(ncc);
% SensB=ApplySCCBySlices(SensBA(:,:,GoodChannelsI),sccmtx,ncc);
projC=gpermute(ApplySCCBySlices(gpermute(projCA(:,:,:,GoodChannelsI),[4 3]),sccmtx,ncc),[4 3]);

GoodChannelsI=1:ncc;

%%
clear recoC
for c=1:ncc %NCha
    disp(c);
    DataP=conj(projC(:,:,:,c));
    recoC(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',TrajectoryP, DataP, SensP);
end
%%
setenv('TOOLBOX_PATH','/home/a/bart-0.4.03')
%%
SensB=RunESPIRiTForSensMaps(recoC,30);

recoCx=CalcSENSE1f(recoC,SensB);


%% Now Multichannel
WhichShots=[1:1:NLin];
% WhichShots=getKrandomSamples(NLin,20);
nWhichShots=numel(WhichShots);
B=ismember(ShotsIdxs,WhichShots);


setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensorClean')

SensP=permute(SensB,[1 2 4 3]);
% SensP=gflip(SensP,1:2);

DataP=gpermute( conj(projC(:,B,:,GoodChannelsI)),[1 2 3 4]);
TrajectoryP=BARTTraj(:,B);
% PICSprm='-m -r:0.0001 -R T:7:0:0.0001 -t';
PICSprm='-m -r:0.0001 -R W:7:0:0.000001 -t';
% PICSprm='-m -r:0.0001 -R T:7:0:0.1 -t';
% recoX = bart('pics -r:0.0001 -R T:7:0:0.001 -t ',TrajectoryP, DataP, SensP);
recoX = bart(['pics ' PICSprm],TrajectoryP, DataP, SensP(:,:,:,GoodChannelsI));
fgmontage(recoX);title(PICSprm);
axis equal
% reco = bart('pics -S -r:0.01 -RT:1024:0:0.0001 -p Weights  -t ',Trajectory, Data, SensP);
%% Now small
% DataP2=DataP(:,BTraj,:,:);
clear recoC2
for c=1:ncc %NCha
    disp(c);
    DataP=conj(projC(:,BTraj,:,c));
    recoC2(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',BARTTraj2, DataP, ones(Sz2));
end

SensB2=RunESPIRiTForSensMaps(recoC2,30);

recoCx2=CalcSENSE1f(recoC2,SensB2);

fgmontage(recoCx2);title('recoCx2')

SensP2=permute(SensB2,[1 2 4 3]);
% SensP=gflip(SensP,1:2);
%%
DataP2=gpermute( conj(projC(:,BTraj,:,:)),[1 2 3 4]);
ShiftInHz=0;
DataP2=RepDotMult(DataP2,exp(1i*2*pi*ShiftInHz*TimeInMs2/1000));

PICSprm='-m -r:0.0001 -R W:7:0:0.000001 -t';
recoX2 = bart(['pics ' PICSprm],BARTTraj2, DataP2, SensP2(:,:,:,:));
fgmontage(recoX2);title(PICSprm);%axis equal
% fgmontage(B0M2);title(PICSprm);axis equal
%% Multichannel NUFMAC

TimeInMs=repmat((0:(NRead-1))*dt*1000,[1 nWhichShots]);
CSinkHz=0.44*7/3;
FatAddedPhase=2*pi*CSinkHz*TimeInMs;

clc
setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensor_NUFMAC')
CS_Bit=13;

SensP1=repmat(gpermute(SensB,[4 3]),[ones(1,CS_Bit) 2]);
SensP2=repmat(gpermute(cat(3,BW2,BWF),[CS_Bit+1 3]),[1 1 1 NCha]);
SensP=SensP1.*SensP2;

DataP=conj(projC);
TrajectoryP=repmat(BARTTraj,[ones(1,CS_Bit) 2]);
% SensP=repmat(SensP,[ones(1,CS_Bit) 2]);
CSFmap=gpermute(cat(3,ones(size(FatAddedPhase)),exp(-1i*FatAddedPhase)),[CS_Bit+1 3]);
CSFmap=repmat(CSFmap,[1 1 1 NCha]);
recoDN = bart('pics -m -r:0.000001 -R T:7:0:0.0001 -t ',TrajectoryP, DataP, SensP,CSFmap);
% recoDN = bart('pics -m -r:0.000001 -R W:7:0:0.00000001 -t ',TrajectoryP, DataP, SensP,CSFmap);
% fgmontage(squeeze(recoDN))
fgmontage(sum(recoDN,CS_Bit+1))
title('Simulating FAT CS NUFMAC Multichannel');
%%
AA=readcfl('/home/a/bart-gpu_tensorClean/AA');

%%
setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensor/')

Sens=RunESPIRiTForSensMaps(imc(:,:,:,1),15);

Outxx=CalcSENSE1f(imc(:,:,:,1),SensB);
%%
figure;clf;
ha = tight_subplot(2,2,[.03 .01],[.01 .03],[.01 .01]);
axes(ha(1));
imagesc(abs(Outxx));title('Spiral recon T1');
axis equal
removeTicks;
axes(ha(2));
imagesc(angle(Outxx),[-pi pi]);
axis equal;removeTicks;
axes(ha(3));
imagesc(abs(recoX));title('BART TV regularized');
axis equal;removeTicks;
axes(ha(4));
imagesc(angle(recoX),[-pi pi]);
axis equal;removeTicks;

colormap gray
%%
figure;
ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
          for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
          set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
%%
B=im(:,:,1)>1.5e-4;
B= imdilate(B,strel('disk',5));
B=imfill(B,'holes');

fgmontage(B)
Supportx=B;

Sz=size(B);
%%
TrajFOV=[kx*Sz(1)/(2*pi) ky*Sz(2)/(2*pi) ];

Datax=CombineDims(data(:,:,:,1,1,1,1,1,1,1,:),[11 1]);
Datay=CombineDims(Datax,[3 1]);
%%
% [a,A,P] = Prepare4Recon(Datay, TrajFOV, Sens, Supportx);
[a,A,P] = Prepare4Recon(Datay, TrajFOV, Sens(:,:,1), true(Sz));

pp=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]);
[a,A,P] = Prepare4Recon(pp', TrajFOV, Sens(:,:,1)*0+1, true(Sz));
x = a./P./P; % starting point

[a,A,P] = Prepare4Recon(proj', TrajFOV, Sens(:,:,1)*0+1, true(Sz));
x = a./P./P; % starting point

%%
%% Basic reconstruction
disp('Basic reconstruction');
lambda = 0*1e-3*max(abs(a(:)));
[x,t,d] = ReconCG(a, @(x) A(x) + lambda*x, a, 30, P);

% err = x-ref;
fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
figure(1);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (CG)');
% figure(2);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (CG) in inverted gray levels');
figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');

%%
moving=abs(Ix(:,:,:,1));
fixed=abs(recoX2);
%%
figure; 
imshowpair(moving, fixed);
title('Unregistered');
%%
[optimizer,metric] = imregconfig('multimodal');

% movingRegisteredDefault= imregister(moving, fixed, 'affine', optimizer, metric);
tform = imregtform(moving,fixed,'affine',optimizer,metric)
movingRegisteredDefault= imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
%%
figure; 
% imshowpair(movingRegisteredDefault, fixed);
imshowpair(movingRegisteredDefault, fixed,'montage');
title('A: Default registration');
%%
B0RealEx=B0_Hz;
B0_Hzx(~isfinite(B0_Hzx))=0;
B0_Hzx=B0RealEx;
B0_Hzs = imwarp(B0_Hzx,tform,'OutputView',imref2d(size(fixed)));
B0RealExs=B0_Hzs;