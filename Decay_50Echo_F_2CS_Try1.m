setPath
%%
clc
clear
close all
%%
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';

SliIdx=36;
Nres=128;
% Nres=96;
Sz=[Nres Nres];

load([BaseP 'Sli' num2str(SliIdx) '.mat']);
CurSli192=padarray(padarray(CurSliS(:,1:184,:,:),[4,8,0,0],'pre'),[17,0,0,0],'post');

CurSliSz=imresizeBySlices(CurSli192,Sz);
kdata_full=permute(fft2cg(CurSliSz),[3 1 2 4]); % 50   128   128    12
%
clear Sens SensS
load([BaseP 'Sli' num2str(SliIdx) '_Sens.mat']);
Sens=double(SensS);
Sens=permute43(Sens);
Sensr=imresizeBySlices(Sens,Sz);

CurSliCombinedr=sum(CurSliSz.*conj(Sensr),4);

nChannels=size(Sens,4);

im_full_ref=permute(sos(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);
im_full_phase=permute(sop(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);
disp('Read data');
%% Undersampling and Generate Calibration data
% kdata=DataWithDecay;
param.dt=9.3000e-04;
param.t0=0.0091;

dt = param.dt; % echo spacing
t0 = param.t0; % % time for first echo

[nt,nx,npe,nc]=size(kdata_full);

numberOFecho = 6;
calib_PEsize=[48,48];
kdata_calib=kdata_full(1:numberOFecho,:,:,:);
kdata_calib=crop(kdata_calib,[numberOFecho,calib_PEsize,size(kdata_calib,4)]);
t0_calib = t0; % time before first echo

TEs_GRE_calib=(0:numberOFecho-1)*dt+t0_calib;
TEs_GRE_calib = TEs_GRE_calib(:);
TEs_GRE=(0:nt-1)*dt+t0;
TEs_GRE=TEs_GRE(:);
disp('Calib data');
%%
addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/pattern_test')

Block_size_y=8;
Block_size_z=8;
mask_sample = EPTI_sampling_mask_Regular_YZ_V83(kdata_full,Block_size_y,Block_size_z);
mask_sample2 = circshift(mask_sample,[0,-2,-1,0]);
mask_sample2 = mask_sample+mask_sample2;
mask_sample = mask_sample2;

mask_sampleP=permute(mask_sample,[2 3 4 5 1]);
kdata_fullP=permute(kdata_full,[2 3 4 5 1]);
kdata=mask_sample.*kdata_full;

y=mask_sampleP.*kdata_fullP;
%% Pre-Processing: smap, background phase, B0-inhomogeneity phase
[ sens_map,P_dB,Phase0 ] = Preprocess_SubspaceRec_Use( kdata_calib,TEs_GRE_calib,[nx,npe] );

% figure; imshow(abs(permute(P_dB(:,:,1),[2 1])),[-100 100]);
% figure; imshow((permute(Phase0,[2 1])),[-pi pi]);
% sens_map=squeeze(sens_map);
%
sens_map=squeeze(Sensr);

FirstEchoCalibMag=sum(ifft2cg(padarray(permute(kdata_calib(1,:,:,:),[2 3 4 1]),(Sz-calib_PEsize)/2,'both')).*conj(sens_map),3);
LastEchoCalibMag=sum(ifft2cg(padarray(permute(kdata_calib(end,:,:,:),[2 3 4 1]),(Sz-calib_PEsize)/2,'both')).*conj(sens_map),3);

FirstEchoCalibMag=sum(imresizeBySlices(ifft2cg(permute(kdata_calib(1,:,:,:),[2 3 4 1])),Sz).*conj(sens_map),3);
LastEchoCalibMag=sum(imresizeBySlices(ifft2cg(permute(kdata_calib(end,:,:,:),[2 3 4 1])),Sz).*conj(sens_map),3);
DEst=log(abs(LastEchoCalibMag)./abs(FirstEchoCalibMag))/(size(kdata_calib,1)-1);
DEst=max(min(DEst,0),-0.03);

Phase0=exp(1i*Phase0);
% figure; imshow(angle(permute(Phase0.*(-P_dB*dt),[2 1])),[-pi pi]);
% figure; imshow(angle(permute(Phase0.*(-P_dB*dt),[2 1])),[-pi pi]);

Phase_T=zeros(nx,npe,nt);
for t=1:nt
    Phase_T(:,:,t)=exp(1i*2*pi*P_dB*TEs_GRE(t)).*Phase0;  
end
disp('Base phase');
%%
nEchos=nt;
TEs_ms=TEs_GRE;
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);
DataWithEchos=CurSliCombinedr;
PhaseChangePerEcho=2*pi*P_dB*dt;
ShowAbsAngle(DataWithEchos)

nCS=2;
%% Check that the data is low rank after (maximal) temporal Hankelization
HankelTemporalLen=3;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(nEchos,HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[Sz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[Sz, numel(HankelMat)]).*DeHankelizingMatP,3));

[~, s_vals] = llr_thresh_OnPreparedBlocks(H_for(DataWithEchos), 0);

fgmontage(s_vals)

%%
% FatArea=SmoothBySlices(s_vals(:,:,2)>0.5,[5 5],1)>0.1;
FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.01,[5 5],1)>0.1;

% FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.1,[5 5],1)>0.1;
% %%
% FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>10.01,[5 5],1)>0.1;
CSMask=cat(4,ones(Sz),FatArea);

SensCSMap=sens_map.*CSMask;
disp('ok CS');
%% Create linear operators
F_for=@(x) fft2cg(x);
F_adj=@(x) ifft2cg(x);
M_for=@(x) x.*mask_sampleP;
M_adj=@(x) x.*mask_sampleP;
S_for=@(x) sum(x.*SensCSMap,4);
S_adj=@(x) sum(x.*conj(SensCSMap),3);

A_for=@(x) M_for(F_for(S_for(x)));
A_adj=@(x) S_adj(F_adj(M_adj(x)));
Aop=OpFromFunc(A_for,A_adj);

DataSize=size(y);
DataSize(3)=size(SensCSMap,3);
ImSize=DataSize;
ImSize(3)=1;
ImSize(4)=nCS;

OperatorTest2(Aop,ImSize,DataSize);

PFac=1;
PFac=20;
% M = Repmat([1, 1, 1, 1, nEchos, 1]);
T=repmat(permute(eye(2),[3 4 1 2]),[1 1 1 1 nEchos 1]);
TD=permute(eye(2),[3 4 1 2]).*permute(NTEs,[5 4 3 2 1]);
M = fmac(ones(nEchos,1), T,3,[4 5]);
P = fmac(NTEs(1:30), 1i*cat(3,T*PFac,TD),3,4:5);
D = fmac(NTEs, TD,3,4:5);
disp('ok operators');
%% Get initial solution and phase wraps
% m0=abs(permute(y(:,:,:,:,1,:),[1 2 3 4 6 5]));
m0=abs(FirstEchoCalibMag);
m0=m0*0+mean(FirstEchoCalibMag(:));
m0=real(m0);
% m0=abs(CurSliCombinedr(:,:,1));
p0=cat(3,angle(Phase_T(:,:,1))/PFac,PhaseChangePerEcho);

m0=cat(3,m0,m0*0);
p0=cat(3,p0(:,:,1),p0(:,:,1)*0,p0(:,:,2),p0(:,:,2));

d0=m0*0-0.02;
% d0=cat(3,d1,d1*0-0.02);

Mm0=M*m0;
expPp0 = exp(P * p0);
expDd0 = exp(D * d0);

Est=Mm0.*expPp0.*expDd0;

MPD=RefValueMPD;
MPD.m=m0;
MPD.p=p0;
MPD.d=d0;

disp('ok initials');

W=[];
%% Create proximal operators
lambdam = repmat(0.003,[1 nCS])*5;
lambdap = kron([0.05 0.05],ones(1,nCS));
lambdad = repmat(0.003,[1 nCS])/4;

WaveletLevel=3;
Pm = wave_threshMS_real('db4', WaveletLevel, lambdam);
Pp = wave_threshMS_real('db6', WaveletLevel, lambdap);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
disp('ok proximals');

dohogwild = 1;

%%
P = fmac(NTEs(1:15), 1i*cat(3,T*PFac,TD),3,4:5);
%%
MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, W, 50, [1 1 1], dohogwild, 1,CurSliCombinedr);
%%
p=MPD.p;
p(:,:,1:2)=p(:,:,1:2)*PFac;
clear W
for i=1:20
    tmp=angle(exp(1i*(p + (rand*2-1)*pi)))-p;
    tmp(:,:,1:2)=tmp(:,:,1:2)/PFac;
    W{i}=tmp;
end
WW=cat(4,W{:});
disp('ok W');
%%
% MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, 10, 50, [1 1 1]*2, dohogwild, 0,CurSliCombinedr);
% %%
% MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, 10, 50, [1 5 1]*2, dohogwild, 0,CurSliCombinedr);
% %%
% MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, 10, 200, [1 1 1]*2, dohogwild, 0,CurSliCombinedr);

%% Proposed phase regularized reconstruction with phase cycling
niter = 1;
ninneriter = [0 0 10];

MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, W, niter, ninneriter, dohogwild, 1);
%%
niter = 1;
ninneriter = [0 10 3];
doplot = 1;
W=[];
MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, W, niter, ninneriter, dohogwild, 1);

%%
niter = 500;
ninneriter = [10 10 10];
W=[];
MPD = mpdrecon_AopNR(y, Aop, M, P, D, Pm, Pp, Pd, MPD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr);
%%
% MPD = mpdrecon_AopNR_sharedr(y, Aop, M, P, D, Pm, Pp, Pd, MPD, W, 500 , dohogwild, 0,CurSliCombinedr);
%%
m=MPD.m;
% m=RecM;
p=MPD.p;
d=MPD.d;

MmF=M*m;
expPpF = exp(P * p);
expDdF = exp(D * d);

EstF=MmF.*expPpF.*expDdF;
EstF2=sum(EstF.*CSMask,4);
ErrGT=grmss(CurSliCombinedr-squeeze(EstF2));
fgmontage(EstF2(:,:,[1 10 24 45]),[0 3.5]);title(ErrGT)
fgmontage(CurSliCombinedr(:,:,[1 10 24 45]),[0 3.5]);title('GT')
%%
fgmontage(m)
fgmontage(p(:,:,1:2))
fgmontage(p(:,:,3:4))
fgmontage(d,[-0.03 0])
%% Call BART given P,D
% y 128   128    12     1    50
% sens_map      128x128x12
% SensCSMap      128x128x12x2  
% PD 128   128     1     2    50
%  1     1     2     2    50
ImSize=size(m);
ImSize16=FillOnesTo16(ImSize);

expPpF = exp(P * p);
expDdF = exp(D * d);

PD=expPpF.*expDdF;

EstF=sum(m.*T,3).*expPpF.*expDdF;

yEst=Aop*EstF;
S_for=@(x) sum(x.*SensCSMap,4);
F_for=@(x) fft2cg(x);
M_for=@(x) x.*mask_sampleP;

SensCSMapFM=bart('fftmod 3',SensCSMap);
yFM=bart('fftmod 3',y);
BaseSP='/autofs/space/daisy_002/users/Gilad/';
ScriptFN=[BaseSP 'MPBDScript.txt'];
QQ=bart(['linopScript ' ScriptFN],ImSize16,m,T,PD,SensCSMap,mask_sampleP);

RecM=bart(['picsS -R W:3:0:0.05 -m ' ScriptFN],ImSize16,y,T,PD,SensCSMap,mask_sampleP);

%
ScriptFN=[BaseSP 'MPBDScriptP.txt'];
TP=permute(T,[1 2 3 6 7 4 5]);
PDP=permute(PD,[1 2 3 6 7 4 5]);
SensCSMapP=permute(SensCSMap,[1 2 3 6 7 4 5]);
mask_samplePP=permute(mask_sampleP,[1 2 3 6 7 4 5]);
yP=permute(y,[1 2 3 6 7 4 5]);
mP=permute(m,[1 2 3 6 7 4 5]);
ImSizeP=size(mP);
ImSize16P=FillOnesTo16(ImSizeP);

QQ=bart(['linopScript ' ScriptFN],ImSize16P,mP,TP,PDP,SensCSMapP,mask_samplePP);
QQ=bart(['linopScript -A ' ScriptFN],ImSize16P,yP,TP,PDP,SensCSMapP,mask_samplePP);

RecM=bart(['picsS -R W:3:0:0.05 -m ' ScriptFN],ImSize16P,yP,TP,PDP,SensCSMapP,mask_samplePP);

%%
m1=MPD.m;
p1=MPD.p;
d1=MPD.d;
%%
niter = 100;
ninneriter = [10 10 10];
doplot = 1;
dohogwild = 1;
W=[];
[m, p,d] = mpdrecon_AopN(y, Aop, M, P, D, Pm, Pp, Pd, m, p, d, W, niter, ninneriter, dohogwild, doplot);

MmF=M*m;
expPpF = exp(P * p);
expDdF = exp(D * d);

EstF=MmF.*expPpF.*expDdF;
EstF2=sum(EstF.*CSMask,4);
grmss(CurSliCombinedr-squeeze(EstF2))
figure, imshow3(m)
figure, imshow3(p)
figure, imshow3(d)