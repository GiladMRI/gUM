%% Real data
sTwixQ=load('/autofs/space/daisy_001/users/data/Gilad/gep_CL/meas_MID01928_FID43869_gre_4echo_24_26_G2/sTwixX.mat');
QQ1=sTwixQ.sTwixX.image(:,:,:,:,18,:,:,:,:,:,:,:,:,:,:,:,:);
sTwixQ2=load('/autofs/space/daisy_001/users/data/Gilad/gep_CL/meas_MID01929_FID43870_gre_4echo_37_26_G2/sTwixX.mat');
QQ2=sTwixQ2.sTwixX.image(:,:,:,:,18,:,:,:,:,:,:,:,:,:,:,:,:);

QSens=load('/autofs/space/daisy_001/users/data/Gilad/gep_CL/meas_MID01928_FID43869_gre_4echo_24_26_G2/Sens.mat');
QSens=QSens.SensB(:,:,:,18,1);

QQ1R=sTwixQ.sTwixX.refscan(:,:,:,:,18,:,:,:,:,:,:,:,:,:,:,:,:);
QQ2R=sTwixQ2.sTwixX.refscan(:,:,:,:,18,:,:,:,:,:,:,:,:,:,:,:,:);

save('ForSplitProxCart.mat','QQ1','QQ2','QSens','QQ1R','QQ2R');
%%
load('ForSplitProxCart.mat','QQ1','QQ2','QSens','QQ1R','QQ2R');
[squeeze(QQ1(15,15,93+(-2:2),1)) squeeze(QQ1R(15,15,12+(-2:2),1))]
[squeeze(QQ2(15,15,93+(-2:2),1)) squeeze(QQ2R(15,15,12+(-2:2),1))]
QQ1(:,:,93+(-11:12),1,1,1,1,:)=QQ1R(:,:,12+(-11:12),1,1,1,1,:);
QQ2(:,:,93+(-11:12),1,1,1,1,:)=QQ2R(:,:,12+(-11:12),1,1,1,1,:);
DataBoth=cat(4,squeeze(QQ1),squeeze(QQ2));
DataBoth=perm32(DataBoth(:,:,:,[1 5 2 6 3 7 4 8]));

Rec1=bart('pics -m -S -R W:3:0:0.000001 ',perm43(perm84(DataBoth)),perm43(QSens));
% 2-undersampled on dim 2 starting from index 2
FirstEchoTime=2.4;
TimeBetweenEchos_ms=1.3;
nEchos=8;
EchoTimes_ms=FirstEchoTime+TimeBetweenEchos_ms*(0:nEchos-1);
EchoTimes_ms3=permute32(EchoTimes_ms);

Sz=gsize(Rec1,1:2);

Sens=perm43(QSens);
SensMsk=grmss(Sens,3:4)>0.05;

%%
AccPerEcho=2;
clear Msk
for i=1:nEchos
    disp(i);
    Msk(:,:,i)=squeeze(bart(['poisson -Y ' num2str(Sz(2)/2) ' -Z ' num2str(1) ' -y ' num2str(AccPerEcho) ' -z ' num2str(1) ' -s ' num2str(rand*100000)]));
end
disp('ok Msk');
MskR=Msk;
Msk=zeros([1 Sz(2) nEchos]);
Msk(1,2:2:end,:)=MskR;
Msk(:,93+(-11:12),:)=1;
Msk=repmat(Msk,[Sz(1) 1 1]);
BaseSP='/autofs/cluster/kawin/Gilad/RealPois_';

sig=perm43(perm74(DataBoth)).*perm73(Msk);
SensCSMap=Sens;
mask_sample=perm73(Msk);
RealData=true;
%% Sim
if(~exist('Ma7TMPBD','var'))
    Ma7TMPBD=load('/autofs/cluster/kawin/Gilad/All_Orientation-0x.mat');
end
RealData=false;
CurSli=permute(Ma7TMPBD.CurSetAll(:,:,:,134),[2 3 1 4]);

C=CurSli(:,:,1).*exp(1i*CurSli(:,:,2));
B0_Hz=CurSli(:,:,4);
B0_Hz(~isfinite(B0_Hz))=0;
% D=CurSli(:,:,3);
T2S_ms=CurSli(:,:,3);
T2S_ms=max(4,min(T2S_ms,300));

nEchos=4;
% FirstEchoTime=4.6;
FirstEchoTime=0;
TimeBetweenEchos_ms=4.74;
EchoTimes_ms=FirstEchoTime+TimeBetweenEchos_ms*(0:nEchos-1);
EchoTimes_ms3=permute32(EchoTimes_ms);
R2S=1./T2S_ms;
disp('OK loaded values for sim');
%%
% N=64;
% Sz=[N N];
Sz=gsize(CurSli,1:2);
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
%% Mask: Here we plan poisson sampling masks

% TotalOSRatio=1.5;
BetterThanCurrent43Ration=0.25;
TotalOSRatio=(nEchos/3)*BetterThanCurrent43Ration; % Try easy scenario
AccPerEcho=nEchos/TotalOSRatio;
SigCenterSize=16;

sAccPerEcho=sqrt(AccPerEcho);

clear Msk
for i=1:nEchos
    disp(i);
    Msk(:,:,i)=squeeze(bart(['poisson -Y ' num2str(Sz(1)) ' -Z ' num2str(Sz(2)) ' -y ' num2str(sAccPerEcho) ' -z ' num2str(sAccPerEcho) ' -C ' num2str(SigCenterSize) ' -s ' num2str(rand*100000)]));
end
disp('ok Msk');
MskA=Msk;
SigCenterSize2=16;
Msk(Sz(1)/2+(-SigCenterSize2:SigCenterSize2),Sz(2)/2+(-SigCenterSize2:SigCenterSize2),:)=1;
BaseSP='/autofs/cluster/kawin/Gilad/Pois_';
Msk=squeeze(readcfl([BaseSP 'Msk']));
%% Mask: Here is simply 6x-accelerate on PE with CAIPI over echos 
% comment or skip this part to keep the poisson masks.
Msk=zeros([Sz nEchos]);
Msk(:,3:6:end,1:2:end)=1;
Msk(:,6:6:end,2:2:end)=1;
Msk(Sz(1)/2+(-SigCenterSize2:SigCenterSize2),Sz(2)/2+(-SigCenterSize2:SigCenterSize2),:)=1;
disp('ok Msk 3/6');
BaseSP='/autofs/cluster/kawin/Gilad/';
%%
% [X Y 1 Sens CS Echos]
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);

nChannels=size(Sens,Ch_Dim);

Cr=Cr.*SensMsk;

disp('Masked');
%%
TEs_ms=EchoTimes_ms.';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

nCS=1;
disp('ok');
%
SensCSMap=Sens;
mask_sample=perm73(Msk);
disp('ok CS');
%% Create linear operators
DataSize=FillOnesTo16(Sz);
DataSize(Ch_Dim)=nChannels;
DataSize(TS_Dim)=nEchos;
ImSize=FillOnesTo16(Sz);
ImSize(TS_Dim)=nEchos;
ImSize(CS_Dim)=nCS;
%
LS_ScriptFN=[pwd filesep 'Cart2CS_ITS.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% Copy with sens/CS mask, sum over CS, FFT and sample mask
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

ImSz16=FillOnesTo16(ImSize);
% Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample);
%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
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

sig=bart(['linopScript -d 5 ' LS_ScriptFN],ImSz16,GT,SensCSMap,mask_sample);
%%
sigFN=[BaseSP 'sig'];
SensFN=[BaseSP 'Sens'];
MskFN=[BaseSP 'Msk'];
ksp_adjFN=[BaseSP 'ksp_adj'];

writecfl(sigFN,sig);
writecfl(SensFN,SensCSMap);
writecfl(MskFN,mask_sample);

if(RealData)
    ksp_adj=bart(['linopScript -A ' LS_ScriptFN],ImSz16,sig,SensFN,MskFN);
else
    ksp_adj=bart(['linopScript -N ' LS_ScriptFN],ImSz16,GT,SensFN,MskFN);
end
disp('Produced GT and sig');

writecfl(ksp_adjFN,ksp_adj);
%%
if(RealData)
    RecA=bart(['picsS -m -S -g -R W:3:64:0.00000001 ' LS_ScriptFN],ImSz16,sigFN,SensFN,MskFN);
else
    RecA=bart(['picsS -m -S -g -R W:3:64:0.0001 ' LS_ScriptFN],ImSz16,sigFN,SensFN,MskFN);
end
RecA=squeeze(RecA).*SensMsk;
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(squeeze(RecA),1:nEchos,TimeBetweenEchos_ms,FirstEchoTime);
UpdatedT2SMap_ms=min(300,max(4,abs(UpdatedT2SMap_ms)));
disp('ok');
if(RealData)
    PDBase0A=PDBase0;
    PDBase0=abs(PDBase).*exp(FirstEchoTime./UpdatedT2SMap_ms).*exp(1i*angle(PDBase0A));
%     
%     PDBase0=PDBase0A.*(grmss(s_vals,3)>1e-4);
%     PDBase0=exp(1i*angle(PDBase0A)).*min(abs(PDBase0),median(abs(PDBase0(abs(PDBase0)>0)))*20);
    UpdatedT2SMap_ms=UpdatedT2SMap_ms*0+50;
end
%%
DA=RecA-GTX;
DAF=Fitted0-GTX;

ErrShowFac=10;
StepARes=cat(4,GTX,RecA,DA*ErrShowFac,DAF*ErrShowFac);

fgmontagex(StepARes,[0 1e3])
%%
Elem0{1}=abs(PDBase0);
Elem0{2}=angle(PDBase0);
Elem0{3}=UpdatedB0Map_Hz;
Elem0{4}=UpdatedT2SMap_ms;% Cs=SmoothBySlices(mGT.*exp(1i*pGT),[20 20],10);
% Elem0{1}=abs(Cs);
% Elem0{2}=angle(Cs);
% Elem0{3}=SmoothBySlices(bGT,[20 20],10);
% Elem0{4}=tGT*0+50;
%%
ninneriterBART=[2 2 2 2];
% ElemsAlphas=[1e-3,1e-8,1e-6,1e-4];
% ElemsLambda=[1e-6,1e-8,1e-8,1e-6];
ElemsAlphas=[1e-3,1e-8,1e-6,1e-4];
ElemsLambda=[1e-6,1e-8,1e-8,1e-6];
ElemsAlphas=[1e-3,1e-8,1e-6,1e-4];
ElemsLambda=[1e-6,1e-3,1e-8,1e-5];
if(RealData)
    ElemsAlphas=[1e-4,1e+2,1e+2,1e+7];
    ElemsLambda=[1e-6,1e-7,0,0];
%     ElemsLambda=[0,0,0,0];
    ninneriterBART=[2 2 0 2];
end
for i=1:numel(Elem0)
    writecfl([BaseSP 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
    writecfl([BaseSP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[Sz 1 1 1 1 1]));
end

writecfl([BaseSP 'ElemsAlpha'],ElemsAlphas.');
writecfl([BaseSP 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([BaseSP 'ElementTypes'],ElementTypes.');

writecfl([BaseSP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([BaseSP 'Elem' num2str(i-1) '.hdr']);
    delete([BaseSP 'Elem' num2str(i-1) '.cfl']);
end

ContinueRun=false;
if(ContinueRun)
    for i=1:numel(Elem0)
        writecfl([BaseSP 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end
%
QQ=bart(['splitProx -i 10000 -s 30 -d 2 -S ' ksp_adjFN ' -g -f -F ' BaseSP ' ' LS_ScriptFN],ImSz16,1,SensFN,MskFN);
%%
ErrVec=readcfl([BaseSP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<1,1))-1);
% figure;plot(ErrVec)
%
Maps=cell(1,4);
for i=1:4
    Maps{i}=readcfl([BaseSP 'Elem' num2str(i-1)]);
end

Mm=M*Maps{1};
expPp = exp(P * Maps{2});
expBb = exp(B * Maps{3});
expTt = exp(TT * (1./Maps{4}));

% Maps10k=Maps;

Rec=Mm.*expPp.*expBb.*expTt;
RecX=squeeze(sum(Rec,CS_Dim));
fgmontagex(RecX)
%%
DA=abs(RecA)-abs(squeeze(Rec1));
DX=abs(RecX)-abs(squeeze(Rec1));

ErrShowFac=10;
QQ=cat(4,squeeze(Rec1),squeeze(RecA),squeeze(RecX),abs(DA)*ErrShowFac,abs(DX)*ErrShowFac);

%%
DX=RecX-GTX;

StepBRes=cat(4,GTX,RecA,RecX,DA*ErrShowFac,DAF*ErrShowFac,DX*ErrShowFac);
fgmontagex(StepBRes,[0 1e3]);title(num2str([grmss(DA) grmss(DAF) grmss(DX)],' %.2f '));

% ErrShowFac=5;
StepBRes=cat(4,RecA,RecX,DA*ErrShowFac,DX*ErrShowFac);
fgmontagex(mean(abs(StepBRes),3),[0 1e3],'Size',[2 2]);

DBase=GTX-gmean(GTX);
Denom=norm(DBase(:));
NRMSEA=norm(DA(:))./Denom;
NRMSEX=norm(DX(:))./Denom;
%%
fgmontagex(StepBRes(:,:,3,:),[0 1e3]);title(num2str([grmss(DA) grmss(DAF) grmss(DX)],' %.2f '));