addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/pattern_test')

cd('/autofs/space/daisy_002/users/Gilad/phase_cyclingD')
setPath
%%
clc
clear
close all
%%
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';

% SliIdx=16;
% Nres=128;
% % Nres=96;
% Sz=[Nres Nres];


% CurSliSz=imresizeBySlices(CurSli192,Sz);
% kdata_full=permute(fft2cg(CurSliSz),[3 1 2 4]); % 50   128   128    12
kdata_full=permute(FullSig,[3 1 2 4]); % 50   128   128    12
%
% clear Sens SensS
% load([BaseP 'Sli' num2str(SliIdx) '_Sens.mat']);
% Sens=double(SensS);
% Sens=permute43(Sens);
% Sensr=imresizeBySlices(Sens,Sz);
Sensr=Sens;

% CurSliCombinedr=sum(CurSliSz.*conj(Sensr),4);
CurSliCombinedr=X.*SensMsk;
nChannels=size(Sens,4);

im_full_ref=permute(sos(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);
im_full_phase=permute(sop(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);
disp('Read data');
%%
% Block_size_y=8;
% Block_size_z=8;
% mask_sample = EPTI_sampling_mask_Regular_YZ_V83(kdata_full,Block_size_y,Block_size_z);
% mask_sample2 = circshift(mask_sample,[0,-2,-1,0]);
% mask_sample2 = mask_sample+mask_sample2;
% mask_sample = mask_sample2;
mask_sample=repmat(perm31(Msk),[1 1 1 nChannels]);

mask_sampleP=permute(mask_sample,[2 3 4 5 1]);

mask_sampleP=mask_sampleP(:,:,1,:,:);

kdata_fullP=permute(kdata_full,[2 3 4 5 1]);
kdata=mask_sample.*kdata_full;

y=mask_sampleP.*kdata_fullP;
disp('Masked');
%%
% nEchos=nt;
% TEs_ms=TEs_GRE;
TEs_ms=EchoTimes_ms.';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);
DataWithEchos=CurSliCombinedr;
% PhaseChangePerEcho=2*pi*P_dB*dt;
ShowAbsAngle(DataWithEchos)

nCS=1;
disp('ok');
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
if(nCS>1)
% FatArea=SmoothBySlices(s_vals(:,:,2)>0.5,[5 5],1)>0.1;
FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.01,[5 5],1)>0.1;

FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.1,[5 5],1)>0.1;
% %%
% FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>10.01,[5 5],1)>0.1;
CSMask=cat(4,ones(Sz),FatArea);
fgmontage(CSMask);
SensCSMap=sens_map.*CSMask;
else
    SensCSMap=sens_map;
end
disp('ok CS');
%% Create linear operators

CS_Dim=4;
Ch_Dim=3;
TS_Dim=5;

F_for=@(x) fft2cg(x);
F_adj=@(x) ifft2cg(x);
M_for=@(x) x.*mask_sampleP;
M_adj=@(x) x.*mask_sampleP;
S_for=@(x) sum(x.*SensCSMap,CS_Dim);
S_adj=@(x) sum(x.*conj(SensCSMap),Ch_Dim);

A_for=@(x) M_for(F_for(S_for(x)));
A_adj=@(x) S_adj(F_adj(M_adj(x)));
Aop=OpFromFunc(A_for,A_adj);

DataSize=size(y);
DataSize(Ch_Dim)=nChannels;
ImSize=DataSize;
ImSize(Ch_Dim)=1;
ImSize(CS_Dim)=nCS;

OperatorTest2(Aop,ImSize,DataSize);
%%
x = randn(ImSize) + 1j*randn(ImSize);
yy = randn(DataSize) + 1j*randn(DataSize);
Ax = Aop*x;
Aty = Aop'*yy;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));

LS_ScriptFN=[pwd filesep 'Cart2CS.txt'];
% SensCSMap is 0
% mask_sampleP is 1
WriteLinopToFile(LS_ScriptFN,{'fmac 0 8','fftc 3','fmac 1 0'});

ImSz16=FillOnesTo16(ImSize);
Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sampleP);
grmss(Ax-Ax_LS)
%%
% M = Repmat([1, 1, 1, 1, nEchos, 1]);
T=repmat(permute(eye(nCS),[3 4 1 2]),[1 1 1 1 nEchos 1]);
TD=permute(eye(nCS),[3 4 1 2]).*permute(NTEs,[5 4 3 2 1]);
TDx=permute(eye(nCS),[3 4 1 2]).*permute((EchoTimes_ms.'),[5 4 3 2 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
% B = fmac(NTEs(1:30), 1i*TD,3,4:5);
B = fmac(NTEs(1:nEchos)/1000, -1i*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
% D = fmac(NTEs, TD,3,4:5);
D = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%% Get initial solution and phase wraps
m0=abs(X_Estr(:,:,1));
p0=angle(X_Estr(:,:,1));

b0=B0_Hz_Estr;
d0=R2S_Estr;

if(nCS>1)
    m0=cat(3,m0,m0*0);
    p0=cat(3,p0(:,:,1),p0(:,:,1)*0);
    b0=cat(3,b0(:,:,1),b0(:,:,1));
    d0=cat(3,d0(:,:,1),d0(:,:,1)*0);
end

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expDd0 = exp(D * d0);

Est0=Mm0.*expPp0.*expBb0.*expDd0;
Est0X=squeeze(sum(Est0,4));

MPBD=RefValueMPBD;
MPBD.m=m0;
MPBD.p=p0;
MPBD.b=b0;
MPBD.d=d0;

disp('ok initials');

W=[];
%% Create proximal operators
lambdam = repmat(0.003,[1 nCS])*5;
lambdap = repmat(0.05,[1 nCS]);
lambdab = repmat(0.05,[1 nCS]);
lambdad = repmat(0.003,[1 nCS]);

WaveletLevel=3;
Pm = wave_threshMS_real('db4', WaveletLevel, lambdam);
Pp = wave_threshMS_real('db6', WaveletLevel, lambdap);
Pb = wave_threshMS_real('db6', WaveletLevel, lambdab);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
disp('ok proximals');

dohogwild = 1;
%%
BaseSP='/autofs/space/daisy_002/users/Gilad/';

ScriptFN=[BaseSP 'MPBDScriptP.txt'];

mP=permute(m0,[1 2 3 6 7 4 5]);
ImSizeP=size(mP);
ImSize16P=FillOnesTo16(ImSizeP);

yP=permute(y,[1 2 3 6 7 4 5]);

SensCSMapP=permute(SensCSMap,[1 2 3 6 7 4 5]);
TP=permute(T,[1 2 3 6 7 4 5]);
mask_samplePP=permute(mask_sampleP,[1 2 3 6 7 4 5]);

BARTS=struct();
BARTS.cmd=['picsS -w 1 -R W:3:0:0.001 -m ' ScriptFN];
BARTS.ImSize16=ImSize16P;
BARTS.yP=yP;
BARTS.TP=TP;
BARTS.Others={SensCSMapP mask_samplePP};

BaseFP='/tmp/';

Flds=setdiff(fieldnames(BARTS),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTS.(Flds{i})))
        for j=1:numel(BARTS.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTS.(Flds{i}){j});
            BARTS.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTS.(Flds{i}));
        BARTS.(Flds{i})=tmpFN;
    end
end

disp('Prepared for BART');
%
niter = 500;
ninneriter = [10 10 10 10];
W=[];
%

BARTS_Aop=struct();
BARTS_Aop.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aop.ImSz16=ImSz16;
BARTS_Aop.Others={SensCSMap mask_sampleP};

BARTStruct=BARTS_Aop;
Flds=setdiff(fieldnames(BARTStruct),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTStruct.(Flds{i})))
        for j=1:numel(BARTStruct.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTStruct.(Flds{i}){j});
            BARTStruct.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTStruct.(Flds{i}));
        BARTStruct.(Flds{i})=tmpFN;
    end
end
BARTS_Aop=BARTStruct;

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,y,BARTS_Aop.Others{:});

disp('Prepared for BART_Aop');
%%
AlphaFacMPBD=struct('m',1,'p',1,'b',1000,'d',1);
MPBD = MPBDrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr,BARTS,10,AlphaFacMPBD);
MPBD = MPBDrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr,BARTS,Inf,AlphaFacMPBD);
MPBD = MPBDrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr,BARTS,1,AlphaFacMPBD);
%%
figure;
subplot(2,2,1);gmontage(MPBD.m,[0 1000]);title('Cur MPBD');
subplot(2,2,2);gmontage(MPBD.p,[-pi pi]);
subplot(2,2,3);gmontage(MPBD.b,[-300 300]);
subplot(2,2,4);gmontage(MPBD.d,[0 0.1]);

figure;
subplot(2,2,1);gmontage(m0,[0 1000]);title('MPBD_0');
subplot(2,2,2);gmontage(p0,[-pi pi]);
subplot(2,2,3);gmontage(b0,[-300 300]);
subplot(2,2,4);gmontage(d0,[0 0.1]);

figure;
subplot(2,2,1);gmontage(abs(X(:,:,1)),[0 1000]);title('MPBD GT');
subplot(2,2,2);gmontage(angle(X(:,:,1)),[-pi pi]);
subplot(2,2,3);gmontage(B0_Hzr,[-300 300]);
subplot(2,2,4);gmontage(R2Sr,[0 0.1]);

fgmontage(b0-MPBD.b,[-5 5])
fgmontage(SB0Var_Hzr,[-5 5])

fgmontage(MPBD.b-B0_Hzr,[-5 5])
fgmontage(b0-B0_Hzr,[-5 5])
%%
Mm = M * MPBD.m;
expPp = exp(P * MPBD.p);
expBb = exp(B * MPBD.b);
expDd = exp(D * MPBD.d);
    
CurEst=Mm.*expPp.*expBb.*expDd;
CurEstX=squeeze(sum(CurEst,4));
%%









%%
% ninneriter = [10 10 10 10]/2;
% MPBD = MPBDrecon_AopNR_BART(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr,BARTS,10);
% %%
% m1=MPBD.m;
% p1=MPBD.p;
% b1=MPBD.b;
% d1=MPBD.d;
% %%
% Mm = M * MPBD.m;
% expPp = exp(P * MPBD.p);
% expBb = exp(B * MPBD.b);
% expDd = exp(D * MPBD.d);
%     
% CurEst=Mm.*expPp.*expBb.*expDd;
% CurEstX=squeeze(sum(CurEst,4));
% Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sampleP);
%%
% ksp_adj1=Aop'*y;

% bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,y,BARTS_Aop.Others{:});
% 
% Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sampleP);
% Ax_LS=bart(['linopScript  ' LS_ScriptFN],BARTS_Aop.ImSz16,x,BARTS_Aop.Others{:});
% AHAx_LS=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,Ax_LS,BARTS_Aop.Others{:});


%% Prepare BART for m,p
BaseSP='/autofs/space/daisy_002/users/Gilad/';

ImSize=size(m0);
ImSize16=FillOnesTo16(ImSize);

b=MPBD.b;
d=MPBD.d;

expDdF = exp(D * d);
expBbF = exp(B * b);

BD=expDdF.*expBbF;

PDP=permute(BD,[1 2 3 6 7 4 5]);
%%














%%
niter = 500;
ninneriter = [10 10 10 10];
% ninneriter = [10 10 10 10]/2;
W=[];
MPBD = MPBDrecon_AopNR(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,CurSliCombinedr);
%%






%%
P = fmac(NTEs(1:15), 1i*cat(3,T*PFac,TD),3,4:5);

P = fmac(NTEs(1:30), 1i*T,3,4:5);

%%
MPBD = MPBDrecon_AopNR(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, 50, [1 1 1], dohogwild, 1,CurSliCombinedr);
%%
p=MPBD.p;
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
ninneriter = [0 0 7 7];

MPBD = MPBDrecon_AopNR(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 1);
%%
niter = 1;
ninneriter = [0 10 5 5];
doplot = 1;
W=[];
MPBD = MPBDrecon_AopNR(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 1);
%%
niter = 1;
ninneriter = [0 0 0 7];

MPBD = MPBDrecon_AopNR(y, Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 1);
%%
m1=m;
p1=p;
b1=b;
d1=d;
%%
m=MPBD.m;
p=MPBD.p;
b=MPBD.b;
d=MPBD.d;

MmF=M*m;
expPpF = exp(P * p);
expBbF = exp(B * b);
expDdF = exp(D * d);

EstF=MmF.*expPpF.*expBbF.*expDdF;
EstF2=sum(EstF.*CSMask,4);
ErrGT=grmss(CurSliCombinedr-squeeze(EstF2));
fgmontage(EstF2(:,:,[1 10 24 45]),[0 3.5]);title(ErrGT)
%%
CurDiff=abs(CurSliCombinedr-EstF2);
fgmontage(CurDiff(:,:,[1 10 24 45]),[-3 3]);title(ErrGT)
%%
MPBD.p=p;
MPBD.m=m;
%%
CurDiff=CurSliCombinedr-squeeze(EstF2);
ShowAbsAngle(cat(4,EstF2(:,:,[1 10 24 45]),CurSliCombinedr(:,:,[1 10 24 45]),CurDiff(:,:,[1 10 24 45])));
%%
fgmontage(CurSliCombinedr(:,:,[1 10 24 45]),[0 3.5]);title('GT')
%%
b(:,:,2)=b(:,:,1);
%%
fgmontage(m)
fgmontage(p)
fgmontage(b)
fgmontage(d,[-0.03 0])
%%
[~, s_valsEst] = llr_thresh_OnPreparedBlocks(H_for(squeeze(EstF2)), 0);

fgmontage(s_valsEst)

%% Call BART given P,D
% y 128   128    12     1    50
% sens_map      128x128x12
% SensCSMap      128x128x12x2  
% PD 128   128     1     2    50
%  1     1     2     2    50
ImSize=size(m);
ImSize16=FillOnesTo16(ImSize);

expDdF = exp(D * d);
expBbF = exp(B * b);

BD=expDdF.*expBbF;

EstF=sum(m.*T,3).*expPpF.*expDdF;

yEst=Aop*EstF;
S_for=@(x) sum(x.*SensCSMap,4);
F_for=@(x) fft2cg(x);
M_for=@(x) x.*mask_sampleP;

BaseSP='/autofs/space/daisy_002/users/Gilad/';

SensCSMapFM=bart('fftmod 3',SensCSMap);
yFM=bart('fftmod 3',y);
ScriptFN=[BaseSP 'MPBDScript.txt'];
QQ=bart(['linopScript ' ScriptFN],ImSize16,m,T,PD,SensCSMap,mask_sampleP);

RecM=bart(['picsS -R W:3:0:0.05 -m ' ScriptFN],ImSize16,y,T,PD,SensCSMap,mask_sampleP);
%%
RecM=bart(['picsS -R W:3:0:0.05 -m ' ScriptFN],ImSize16P,yP,TP,PDP,SensCSMapP,mask_samplePP);
%%
RecM=bart(BARTS.cmd,BARTS.ImSize16,BARTS.yP,BARTS.TP,PDP,BARTS.Others{:});
%%
m=abs(RecM);
p=angle(RecM);

MmF=M*m;
expPpF = exp(P * p);

EstF=MmF.*expPpF.*expBbF.*expDdF;
EstF2=sum(EstF.*CSMask,4);
ErrGT=grmss(CurSliCombinedr-squeeze(EstF2))

%%
QQ=bart(['linopScript ' ScriptFN],ImSize16P,mP,TP,PDP,SensCSMapP,mask_samplePP);
QQ=bart(['linopScript -A ' ScriptFN],ImSize16P,yP,TP,PDP,SensCSMapP,mask_samplePP);


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