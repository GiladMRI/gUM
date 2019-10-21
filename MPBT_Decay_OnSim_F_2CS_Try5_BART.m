%% Sim
if(~exist('Ma7TMPBD','var'))
    Ma7TMPBD=load('/autofs/cluster/kawin/Gilad/All_Orientation-0x.mat');
end

CurSli=permute(Ma7TMPBD.CurSetAll(:,:,:,134),[2 3 1 4]);

C=CurSli(:,:,1).*exp(1i*CurSli(:,:,2));
B0_Hz=CurSli(:,:,4);
B0_Hz(~isfinite(B0_Hz))=0;
% D=CurSli(:,:,3);
T2S_ms=CurSli(:,:,3);
T2S_ms=max(4,min(T2S_ms,300));

nEchos=30;
TimeBetweenEchos_ms=2;
EchoTimes_ms=TimeBetweenEchos_ms*(0:nEchos-1);
EchoTimes_ms3=permute32(EchoTimes_ms);
R2S=1./T2S_ms;
disp('OK loaded values for sim');
%%
N=64;
Sz=[N N];
%%
B0_Hzr=imresizeBySlices(B0_Hz,Sz);
T2S_msr=imresizeBySlices(T2S_ms,Sz);
T2S_msr=max(4,min(T2S_msr,300));

R2Sr=1./T2S_msr;
R2Sr=max(1e-4,min(0.1,R2Sr));
Cr=imresizeBySlices(C,Sz);
%%
if(~exist('Ma7TSensCollection','var'))
    Ma7TSensCollection=load('/autofs/cluster/kawin/Gilad/CCSensMaps.mat');
end
Sens=perm43(imresizeBySlices(perm31(Ma7TSensCollection.SensCC(:,:,:,15)),Sz));
SensMsk=grmss(Sens,3:4)>0.05;
disp('Loaded sens');
%%
TotalOSRatio=1.5;
TotalOSRatio=4; % Try easy scenario
AccPerEcho=nEchos/TotalOSRatio;
SigCenterSize=8;

sAccPerEcho=sqrt(AccPerEcho);

clear Msk
for i=1:nEchos
    disp(i);
    Msk(:,:,i)=squeeze(bart(['poisson -Y ' num2str(N) ' -Z ' num2str(N) ' -y ' num2str(sAccPerEcho) ' -z ' num2str(sAccPerEcho) ' -C ' num2str(SigCenterSize) ' -s ' num2str(rand*100000)]));
end
disp('ok Msk');
%%
addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/pattern_test')

cd('/autofs/space/daisy_002/users/Gilad/phase_cyclingD')
setPath
%%
clc
clear
close all
%%
% [X Y 1 Sens CS Echos]
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);

BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';

BaseSP='/autofs/space/daisy_002/users/Gilad/';

% kdata_full=gpermute(FullSig,[1 2 5 4 6 3]);
% XP=permute(X,[1 2 TS_Dim:-1:3]);
Sensr=Sens;
mask_sample=gpermute(Msk,[TS_Dim 3]);

% DataWithEchos=XP.*SensMsk;

nChannels=size(Sens,Ch_Dim);

% sig=mask_sample.*kdata_full;
disp('Masked');
%%
TEs_ms=EchoTimes_ms.';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

nCS=1;
disp('ok');
%% Check that the data is low rank after (maximal) temporal Hankelization
% HankelTemporalLen=3;
% [~, ~, ~,H]=ghankel(nEchos,HankelTemporalLen,Sz);
% [~, s_vals] = llr_thresh_OnPreparedBlocks(H*squeeze(DataWithEchos), 0);
% fgmontage(s_vals)
%%
if(nCS>1)
    % FatArea=SmoothBySlices(s_vals(:,:,2)>0.5,[5 5],1)>0.1;
    FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.01,[5 5],1)>0.1;
    
    FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.1,[5 5],1)>0.1;
    % %%
    % FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>10.01,[5 5],1)>0.1;
    CSMask=cat(4,ones(Sz),FatArea);
    fgmontage(CSMask);
    SensCSMap=Sensr.*CSMask;
else
    SensCSMap=Sensr;
end
disp('ok CS');
%% Create linear operators
F_for=@(x) fft2cg(x);
F_adj=@(x) ifft2cg(x);
M_for=@(x) x.*mask_sample;
M_adj=@(x) x.*mask_sample;
S_for=@(x) sum(x.*SensCSMap,CS_Dim);
S_adj=@(x) sum(x.*conj(SensCSMap),Ch_Dim);

A_for=@(x) M_for(F_for(S_for(x)));
A_adj=@(x) S_adj(F_adj(M_adj(x)));
Aop=OpFromFunc(A_for,A_adj);

DataSize=FillOnesTo16(Sz);
DataSize(Ch_Dim)=nChannels;
DataSize(TS_Dim)=nEchos;
ImSize=FillOnesTo16(Sz);
ImSize(TS_Dim)=nEchos;
ImSize(CS_Dim)=nCS;

OperatorTest2(Aop,ImSize,DataSize);
%%
x = randn(ImSize) + 1j*randn(ImSize);
yy = randn(DataSize) + 1j*randn(DataSize);
Ax = Aop*x;
Aty = Aop'*yy;
Out=abs(x(:)'*Aty(:) - conj(yy(:)'*Ax(:)));

LS_ScriptFN=[pwd filesep 'Cart2CS_ITS.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% Copy with sens/CS mask, sum over CS, FFT and sample mask
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

ImSz16=FillOnesTo16(ImSize);
Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample);
grmss(Ax-Ax_LS)
%%
% M = Repmat([1, 1, 1, 1, nEchos, 1]);
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
% B = fmac(NTEs(1:30), 1i*TD,3,4:5);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
% D = fmac(NTEs, TD,3,4:5);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%% Get initial solution and phase wraps
mGT=abs(Cr(:,:,1));
pGT=angle(Cr(:,:,1));

bGT=B0_Hzr;
tGT=T2S_msr;

MmGT=M*mGT;
expPpGT = exp(P * pGT);
expBbGT = exp(B * bGT);
expTtGT = exp(TT * (1./tGT));

GT=MmGT.*expPpGT.*expBbGT.*expTtGT;
GTX=squeeze(sum(GT,CS_Dim));

sig=bart(['linopScript ' LS_ScriptFN],ImSz16,GT,SensCSMap,mask_sample);

disp('Produced GT and sig');
%%
m0=mGT;
p0=pGT;

b0=bGT;
t0=tGT;
%% That's it for SplitProx_MPBT_OnSim_Cart_1CS
%%
b0=SmoothBySlices(bGT,[20 20],3);
t0=SmoothBySlices(tGT,[20 20],3);

if(nCS>1)
    m0=cat(CS_Dim,m0,m0*0);
    p0=cat(CS_Dim,p0,p0*0);
    b0=cat(CS_Dim,b0,b0);
    t0=cat(CS_Dim,t0,t0*0);
end

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Est0=Mm0.*expPp0.*expBb0.*expTt0;
Est0X=squeeze(sum(Est0,CS_Dim));

MPBT=RefValueMPBT;
MPBT.m=m0;
MPBT.p=p0;
MPBT.b=b0;
MPBT.t=t0;

disp('ok initials');
%% Create proximal operators
lambdam = repmat(0.003,[1 nCS])*5;
lambdap = repmat(0.05,[1 nCS]);
lambdab = repmat(0.0,[1 nCS]);
lambdat = repmat(0.003,[1 nCS]);

WaveletLevel=3;
Pm = wave_threshMS_real('db4', WaveletLevel, lambdam);
Pp = wave_threshMS_real('db6', WaveletLevel, lambdap);
Pb = wave_threshMS_real('db6', WaveletLevel, lambdab);
Pb = wave_threshMS_real('db4', WaveletLevel, lambdab);
Pt = wave_threshMS_real('db4', WaveletLevel, lambdat);
disp('ok proximals');

dohogwild = 1;
%%
% ScriptFN=[BaseSP 'MPBDScriptP.txt'];
ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
WriteLinopToFile(ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;
BARTS=struct();
BARTS.cmd=['picsS -w 1 -R W:3:0:0.001 -m ' ScriptFN];
BARTS.ImSz16=MPImSize16;
BARTS.sig=sig;
BARTS.Others={SensCSMap mask_sample};

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

BARTS_Aop=struct();
BARTS_Aop.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aop.ImSz16=ImSz16;
BARTS_Aop.Others={SensCSMap mask_sample};

BARTStruct=BARTS_Aop;
Flds=setdiff(fieldnames(BARTStruct),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTStruct.(Flds{i})))
        for j=1:numel(BARTStruct.(Fllds{i}))
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

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

disp('Prepared for BART_Aop');

niter = 500;
ninneriter = [10 10 10 10];
W=[];
%%
% ninneriter = [0 0 0 10];
ninneriter = [0 0 10 0];
AlphaFacMPBT=struct('m',1,'p',1,'b',10,'t',1e7);
MPBT = MPBTrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, TT, Pm, Pp, Pb, Pt, MPBT, W, niter, ninneriter, dohogwild, 0,GT,BARTS,Inf,AlphaFacMPBT);
%%
% MPBD = MPBDrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,DataWithEchos,BARTS,10,AlphaFacMPBD);
% %%
% MPBD = MPBDrecon_AopNR_BART_LS(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,DataWithEchos,BARTS,1,AlphaFacMPBD);
%%
B0Ss=cat(3,bGT,MPBT.b,b0);
B0Ss=B0Ss-bGT;
fgmontagex(B0Ss,[-30 30],'Size',[1 3]);title('GT                              T-Decay cycling                            Start point','FontSize',16);
%%
T2Ss=cat(3,tGT,MPBT.t,t0);
fgmontagex(T2Ss,[0 100],'Size',[1 3]);title('GT                              T-Decay cycling                            Start point','FontSize',16);
%%
ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
ElemsWS={m0,p0,b0,t0};
for i=1:numel(ElemsWS)
    writecfl([ToBARTP 'ElemsWS_' num2str(i)],ElemsWS{i});
end
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
for i=1:numel(ElemsL)
    writecfl([ToBARTP 'ElemsL_' num2str(i)],ElemsL{i});
end
Mm=M*m0;
expPp0 = exp(P * p);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Othersp=Mm.*expBb0.*expTt0;
Othersb=Mm.*expPp0.*expTt0;
Otherst=Mm.*expPp0.*expBb0;
ElemsAlphas=[AlphaFacMPBT.m/ lipschitz(M),AlphaFacMPBT.p/ lipschitz(P) / (max(abs(Othersp(:)))^2 + eps),...
    AlphaFacMPBT.b / lipschitz(B) / (max(abs(Othersb(:)))^2 + eps),...
    AlphaFacMPBT.t / lipschitz(TT) / (max(abs(Otherst(:)))^2 + eps)];

% alphatbase=1e11;
alphatbase=1e7;
ElemsAlphas(4)=alphatbase / lipschitz(TT) / (max(abs(Otherst(:)))^2 + eps);
writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas);
writecfl([ToBARTP 'ninneriter'],ninneriter);

ElemsLambda=[0.015,0.05,0.05,0.003];
% ElemsLambda=[0.015,0.05,0.05,300];

lambdat = repmat(ElemsLambda(4),[1 nCS]);
Pt = wave_threshMS_real('db4', WaveletLevel, lambdat);

writecfl([ToBARTP 'ElemsLambda'],ElemsLambda);

writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[0 0 0 10];
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');
%%
m=m0;
p=p0;
b=b0;
t=t0;

Mm = M * m;
expPp = exp(P * p);
expBb = exp(B * b);
expTt = exp(TT * (1./t));
Others=Mm.*expPp.*expBb;
% Now T-decay
%     Mm = M * m;
%     Mm=Mm.*expPp.*expBb;
%     alphat = 1.0*AlphaFacMPBT.t / lipschitz(T) / (max(abs(Mm(:)))^2 + eps) * h;
%     disp('t iters');
%     for itinner = 1:ninneriter(4)
%         expTt = exp(T * (1./t));
%         CurEst=Mm .* expTt;
CurEst=Others .* expTt;
AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
r=ksp_adj- AHACurEst;
rConjOthers=r.*conj(Others);
cexpTt=conj(expTt);
ToAdjLinop=(rConjOthers .* conj(expTt));
AfterAdjLinop=(TT' * ToAdjLinop);
dut=-((1./t).^2);
h=1;

% alphat = 1.0*AlphaFacMPBT.t / lipschitz(TT) / (max(abs(Mm(:)))^2 + eps) * h;
alphat = 1.0*alphatbase / lipschitz(TT) / (max(abs(Mm(:)))^2 + eps) * h;
% ForProx=t + alphat * real((T' * (conj(Mm) .* conj(expTt) .* r)) .* (-((1./t).^2)) );
ForReal=AfterAdjLinop .* dut;
AfterReal=real(ForReal);
ForProx=t + alphat * AfterReal;
AfterProx=Pt(ForProx,alphat);
disp('ok');
%         t = Pt(t + alphat * real(-  (T' * (conj(Mm) .* conj(expTt) .* r)) .* ((1./t).^2)    ), alphat);
%         t = Pt(t + alphat * real((T' * (conj(Mm) .* conj(expTt) .* r)) .* (-((1./t).^2)))    ), alphat);
%%
QQ=bart(['splitProx -i 200 ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
XDimsOut=readcfl([ToBARTP 'XDimsOut']);
mDimsOut=readcfl([ToBARTP 'mDimsOut']);
%
Maps=cell(1,4);
for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
%%
% ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

% alpham = 1.0 *AlphaFacMPBT.m/ lipschitz(M);
% alphap = 1.0 *AlphaFacMPBT.p/ lipschitz(P) / (max(abs(Mm(:)))^2 + eps);
% alphab = 1.0*AlphaFacMPBT.b / lipschitz(B) / (max(abs(Mm(:)))^2 + eps);
% alphat = 1.0*AlphaFacMPBT.t / lipschitz(T) / (max(abs(Mm(:)))^2 + eps);
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
subplot(2,2,1);gmontage(mGT,[0 1000]);title('MPBD GT');
subplot(2,2,2);gmontage(pGT,[-pi pi]);
subplot(2,2,3);gmontage(bGT,[-300 300]);
subplot(2,2,4);gmontage(dGT,[0 0.1]);
%%
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
CurEstX=sum(CurEst,4);

grmss(CurEstX-DataWithEchos)/grmss(DataWithEchos)
%%
MPBD1={MPBD.m,MPBD.p,MPBD.b,MPBD.d};








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