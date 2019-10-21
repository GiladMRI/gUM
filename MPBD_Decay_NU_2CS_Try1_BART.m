addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/pattern_test')

cd('/autofs/space/daisy_002/users/Gilad/phase_cyclingD')
setPath
%% For MLN?
% save('ForMLN.mat','TrajM','Sensr','sig','UpdatedB0Map_RS','UpdatedT2SMap_ms_RS');
% PDBase0_RS
%% For fixing P
% fgmontage(angle(exp(1i*(MPBD.p-angle(PDBase0_RS(:,:,1))))))
%%
CurSlizLoc=CurzLocs(Ord(SlicesToRead));
[MinRefD, CorrespondingRefSliIdx]=min(abs(RefzLocs-CurSlizLoc));

RefPDBase0=rot90(padarray(RefMaps.PDBase0(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
RefT2SM=rot90(padarray(RefMaps.UpdatedT2SMap_ms(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
RefPDBase0=min(abs(RefPDBase0),grmss(RefPDBase0)*6).*exp(1i*angle(RefPDBase0));
RefMidEcho=RefPDBase0.*exp(-EchoTimes_ms(CurEchoHere)./RefT2SM);
%%
RepSetIToShow=2;
CurEchoHere=floor(nEchos/2);
SeveralRecs=cat(3,THLRMultiShot_RS(:,:,7,RepSetIToShow),Fitted0_RS(:,:,7,RepSetIToShow),Rec_CompgB0_RS_MX(:,:,10,RepSetIToShow),...
    Fitted0_RS_LLR(:,:,10,RepSetIToShow),CurEst0bXs(:,:,CurEchoHere),CurEstXs(:,:,CurEchoHere));
fgmontagex(SeveralRecs);
QQ=get(gca,'Children');
WSz=size(QQ.CData);
WSzX=WSz./[2 3];

StartX=10;
StartY=30;
text(StartX,StartY,'THLR','FontSize',16,'Color','red')
text(StartX+WSzX(1),StartY,'THLR Fitted','FontSize',16,'Color','red')
text(StartX+WSzX(1)*2,StartY,'LLR','FontSize',16,'Color','red')
text(StartX,StartY+WSzX(2),'LLR Fitted','FontSize',16,'Color','red')
text(StartX+WSzX(1),StartY+WSzX(2),'After MP opt','FontSize',16,'Color','red')
text(StartX+WSzX(1)*2,StartY+WSzX(2),'Cycled','FontSize',16,'Color','red')

fgmontagex(RefMidEcho);title('From GRE-ME');
%%
load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms/meas_MID03499_FID20468_gSpi2d/For_NU_MPBD3_S3.mat');
% BaseEstX=Fitted0_RS_LLR(:,:,:,2);
BaseEstX=Fitted0_RS_LLRS(:,:,:,10);

% UpdatedT2SMap_ms1=UpdatedT2SMap_ms_RS_LLR(:,:,2);
% UpdatedB0Map1=UpdatedB0Map_RS_LLR(:,:,2);
UpdatedT2SMap_ms1=UpdatedT2SMap_ms_RS_LLRS(:,:,2,SliI);
UpdatedB0Map1=UpdatedB0Map_RS_LLRS(:,:,2,SliI);

% PDBase=PDBase0_RS_LLR(:,:,2);
PDBase=PDBase0_RS_LLRS(:,:,SliI);
%%
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
CurTraj=STraj3MMed(:,:,CurReps);
CurSig=QQ.CurSig(1,:,CurRep,:);
%%

Sensr=CurSens;
mask_sample=CurTraj;
CurTSC=1;
sig=CurSig;

nChannels=size(Sensr,Ch_Dim);
Sz=gsize(Sensr,1:2);
%%
Mx=grmss(BaseEstX)*4;
ShowAbsAngle(BaseEstX,1,[0 Mx])

SensMsk=grmss(Sensr,3:4)>0.05;

[Out B BN]=CalcSlicesSNR(grmss(BaseEstX(:,:,1),3:30),true,7);
Msk=imfillholesBySlices(~BN).*SensMsk;
se = strel('disk', 3);
DMsk=imdilate(Msk,se,'same');
%%
nEchos=80;
%%
nTraj=size(CurTraj,2);

TSBMed=GetTSCoeffsByLinear(nTraj,nEchos);
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

TrajPartMed=1:nTraj;
Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];
%%
clear STraj3MMed
for CurRep=1:size(sig,3)
    disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
    SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
    KernsByRepMed{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRep,TSBMed);
    KernsPMedC{CurRep}=permute(KernsByRepMed{CurRep},[1 2 7 6 5 4 3]);
end
KernsPMMed=cat(3,KernsPMedC{:});
clear KernsPMedC
disp('Prepared Med');
%%
% [X Y 1 Sens CS Echos]
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';

BaseSP='/autofs/space/daisy_002/users/Gilad/';
%%
EchoTimes_ms=linspace(0,AcqDwellTime_us*nTraj/1000,nEchos);
%%
InnerTSDiff2_ms=TotalAcqTime_ms/(nEchos-1);
EchoTimes2_ms=TE0_ms+linspace(0,TotalAcqTime_ms,nEchos);
EchoTimes_ms=EchoTimes2_ms;
EchoTimes2_ms3=perm32(EchoTimes2_ms);
R1TSCx2=exp(-EchoTimes2_ms3./UpdatedT2SMap_ms1).*exp(-1i*2*pi*UpdatedB0Map1.*EchoTimes2_ms3/1e3);
X_Estr2=R1TSCx2.*PDBase.*DMsk;
X_Estr2(isnan(X_Estr2))=0;
X_Estr=X_Estr2;
disp('X_Estr estimated');
TEs_ms=EchoTimes_ms.';

% TEs_ms=linspace(0,AcqDwellTime_us*nTraj/1000,nEchos).';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

nCS=1;
disp('ok');
%%
if(nCS>1)
    % FatArea=SmoothBySlices(s_vals(:,:,2)>0.5,[5 5],1)>0.1;
    FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.01,[5 5],1)>0.1;
    
    FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>0.1,[5 5],1)>0.1;
    % %%
    % FatArea=SmoothBySlices((s_vals(:,:,2)./s_vals(:,:,1))>10.01,[5 5],1)>0.1;
    CSMask=cat(CS_Dim,ones(Sz),FatArea);
    fgmontage(CSMask);
    SensCSMap=Sensr.*CSMask;
else
    SensCSMap=Sensr;
end
disp('ok CS');
%% Create linear operators
DataSize=FillOnesTo16(Sz);
DataSize(Ch_Dim)=nChannels;
DataSize(TS_Dim)=nEchos;
ImSize=FillOnesTo16(Sz);
ImSize(TS_Dim)=nEchos;
ImSize(CS_Dim)=nCS;
%
% LS_ScriptFN=[pwd filesep 'Cart2CS_ITS.txt'];
LS_ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/nuftAllTSC_N.txt';

CurTSB=TSBPMed;
CurKerns=sum(KernsPMMed,3);

ImSz16=FillOnesTo16(ImSize);

x=rand(ImSz16)+1i*rand(ImSz16);

disp('a');
% Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample,CurTSB,CurTSC,CurKerns);
disp('b');
AHAx_LS=bart(['linopScript -N ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample,CurTSB,CurTSC,CurKerns);
disp('c');
%%
T=grepmat(1,nEchos,TS_Dim);
TDx=gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,[],[TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,[],[TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,[],[TS_Dim]);
D = fmac(NTEs, -TDx,[],[TS_Dim]);
disp('ok operators');
%% Get initial solution and phase wraps
% m0=abs(X_Estr(:,:,1));
% p0=angle(X_Estr(:,:,1));

m0=abs(Fitted0_MLN(:,:,1)).*grmss(X_Estr(:,:,1))./grmss(Fitted0_MLN(:,:,1));
p0=angle(Fitted0_MLN(:,:,1));

% b0=B0_Hz_Estr;
% d0=R2S_Estr;
% 
% b0=UpdatedB0Map1;
% d0=1./UpdatedT2SMap_ms1;

b0=RefB0Mr;
d0=1./UpdatedT2SMap_ms_MLN;

d0=min(d0,0.1);
d0=max(d0,1e-4);

if(nCS>1)
    m0=cat(CS_Dim,m0,m0*0);
    p0=cat(CS_Dim,p0,p0*0);
    b0=cat(CS_Dim,b0,b0);
    d0=cat(CS_Dim,d0,d0*0+mean(d0(:)));
end

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expDd0 = exp(D * d0);

Est0=Mm0.*expPp0.*expBb0.*expDd0;
Est0X=sum(Est0,CS_Dim);

MPBD=RefValueMPBD;
MPBD.m=m0;
MPBD.p=p0;
MPBD.b=b0;
MPBD.d=d0;

disp('ok initials');
%% Create proximal operators
lambdam = repmat(15e-4,[1 nCS]);
lambdap = repmat(1e-6,[1 nCS]);
lambdab = repmat(1e-6,[1 nCS]);
lambdad = repmat(1e-5,[1 nCS]);

WaveletLevel=3;
Pm = wave_threshMS_real('db4', WaveletLevel, lambdam);
Pp = wave_threshMS_real('db6', WaveletLevel, lambdap);
Pb = wave_threshMS_real('db6', WaveletLevel, lambdab);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
disp('ok proximals');

dohogwild = 1;
%%
LS_ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/nuftAllTSC_N.txt';

MP_ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/nuftMP_N.txt';

% RepsToUse=1:3;
RepsToUse=1;
CurKerns=sum(KernsPMMed(:,:,RepsToUse,:,:,:,:,:),3);

% ScriptFN=[BaseSP 'MPBDScriptP.txt'];
% ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
% MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
% WriteLinopToFile(ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;
BARTS=struct();
BARTS.cmd=['picsS -w 1 -d 5 -R W:3:0:0.001 -m ' MP_ScriptFN];
BARTS.ImSz16=MPImSize16;
BARTS.sig=sig(:,:,RepsToUse,:);
BARTS.Others={SensCSMap mask_sample(:,:,RepsToUse) CurTSB CurTSC CurKerns};

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
disp('Prepared for BART MP');
%%
BARTS_Aop=struct();
BARTS_Aop.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aop.ImSz16=ImSz16;
BARTS_Aop.Others={SensCSMap mask_sample(:,:,RepsToUse) CurTSB CurTSC CurKerns};

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

disp('ksp_adj:');
ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig(:,:,RepsToUse,:),BARTS_Aop.Others{:});
disp('ok ksp_adj:');
%%
disp('AHAEst0:');
AHAEst0= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,Est0,BARTS_Aop.Others{:});
disp('Prepared for BART_Aop');
%%
MFacX=grmss(ksp_adj)/grmss(AHAEst0);
m0=m0*MFacX;
MPBD.m=m0;

disp('Applied m fac');
%
niter = 500;
% ninneriter = [10 10 10 10];
ninneriter = [10 10 10 10]/2;
W=[];
MPBD0={MPBD.m,MPBD.p,MPBD.b,MPBD.d};
%%
AlphaFacMPBD=struct('m',1e-3,'p',1e+0,'b',1e-1,'d',1e-3);
%%
clear Relr
Relr=0;
global Relr
Relr=0;
%% BART for MP: 200-400 sec
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, 1, [0 0 0 0], dohogwild, 0,[],BARTS,1,AlphaFacMPBD);
MPBD0b={MPBD.m,MPBD.p,MPBD.b,MPBD.d};
%%
AlphaFacMPBD=struct('m',1e-2,'p',1e+0,'b',1e-2,'d',1e-1);
% lambdad = repmat(1e-4,[1 nCS]);
% lambdad = repmat(3e-5,[1 nCS]);
% lambdad = repmat(1e-5,[1 nCS]);
lambdad = repmat(1e-6,[1 nCS]);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
% lambdam = repmat(3e-4,[1 nCS]);
lambdam = repmat(3e-5,[1 nCS]);
Pm = wave_threshMS_real('db4', WaveletLevel, lambdam);
lambdap = repmat(1e-5,[1 nCS]);
Pp = wave_threshMS_real('db6', WaveletLevel, lambdap);
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,Inf,AlphaFacMPBD);
%%
niter=500;

for gIter=1:10
    NW=400;
    pvals=rand(1,NW)*2*pi-pi;
    ivals=randi(3,1,NW)-2;
    for ww=1:NW
        W{ww}=(angle(exp(1i*(MPBD.p+pvals(ww))))>0)*2*pi*ivals(ww);
    end
    %     W=gmat2cell(linspace(-pi,pi,300));
    [MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,Inf,AlphaFacMPBD);
end
%%
AlphaFacMPBD=struct('m',1e-3,'p',1e+0,'b',1e-1,'d',1e-1);
ninneriter = [0 0 0 10]/2;
lambdad = repmat(0.0001,[1 nCS]);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,Inf,AlphaFacMPBD);
%%
AlphaFacMPBD=struct('m',1e-3,'p',1e+0,'b',1e-1,'d',1e-1);
% ninneriter = [0 0 0 10]/2;
ninneriter = [10 10 0 10]/2;
ninneriter = [2 2 0 5];
lambdad = repmat(0.00003,[1 nCS]);
Pd = wave_threshMS_real('db4', WaveletLevel, lambdad);
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,Inf,AlphaFacMPBD);
%%
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,10,AlphaFacMPBD);
%%
[MPBD, Relr] = MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, 0,[],BARTS,1,AlphaFacMPBD);
%%
figure;plot(Relr);
%%
figure;
subplot(2,2,1);gmontage(MPBD.m,[0 1e-3]);title('Cur MPBD');
subplot(2,2,2);gmontage(angle(exp(1i*(MPBD.p-pi*1.5))),[-pi pi]);
subplot(2,2,3);gmontage(MPBD.b,[-300 300]);
subplot(2,2,4);gmontage(MPBD.d,[0 0.1]);
subplot(2,2,4);gmontage(Msk./MPBD.d,[0 100]);% colormap hot

Mm = M * MPBD.m;
expPp = exp(P * MPBD.p);
expBb = exp(B * MPBD.b);
expDd = exp(D * MPBD.d);
    
CurEst=Mm.*expPp.*expBb.*expDd;
CurEstX=sum(CurEst,CS_Dim);
CurEstXs=squeeze(CurEstX);
%%
figure;
subplot(2,2,1);gmontage(MPBD0b{1},[0 1e-3]);title('MPBD_0 after MP opt');
subplot(2,2,2);gmontage(MPBD0b{2},[-pi pi]);
subplot(2,2,3);gmontage(MPBD0b{3},[-300 300]);
% subplot(2,2,4);gmontage(MPBD0b{4},[0 0.1]);
subplot(2,2,4);gmontage(1./MPBD0b{4},[0 100]);

figure;
subplot(2,2,1);gmontage(m0,[0 1e-3]);title('MPBD_0');
subplot(2,2,2);gmontage(p0,[-pi pi]);
subplot(2,2,3);gmontage(b0,[-300 300]);
subplot(2,2,4);gmontage(d0,[0 0.1]);
subplot(2,2,4);gmontage(1./d0,[0 100]);

%%
MPBD.m=MPBD.m.*DMsk;
MPBD.p=MPBD.p.*DMsk;
MPBD.b=MPBD.b.*DMsk;
MPBD.d=MPBD.d.*DMsk;
%%
figure;
subplot(2,2,1);gmontage(abs(X(:,:,1)),[0 1000]);title('MPBD GT');
subplot(2,2,2);gmontage(angle(X(:,:,1)),[-pi pi]);
subplot(2,2,3);gmontage(B0_Hzr,[-300 300]);
subplot(2,2,4);gmontage(R2Sr,[0 0.1]);
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
CurEstXs=squeeze(CurEstX);

grmss(CurEstX-DataWithEchos)/grmss(DataWithEchos)
%%
CurEst0b=(M*MPBD0b{1}).*(exp(P * MPBD0b{2})).*(exp(B * MPBD0b{3})).*(exp(D * MPBD0b{4}));
CurEst0bX=sum(CurEst0b,4);
CurEst0bXs=squeeze(CurEst0bX);
%%
MPBD1={MPBD.m,MPBD.p,MPBD.b,MPBD.d};
%%
MPBD.m=MPBD0b{1};
MPBD.p=MPBD0b{2};
MPBD.b=MPBD0b{3};
MPBD.d=MPBD0b{4};
%%



OutGifFN=[mainP filesep 'S_' num2str(SlicesToRead) '.gif'];
delete(OutGifFN);
Mx=grmss(CurEstXs)*6;
figure;pause(1)
imagesc(abs(CurEstXs(:,:,1)),[0 Mx]);colormap gray;removeTicks;pause(1)
gif(OutGifFN)
for i=2:size(CurEstXs,3)
    imagesc(abs(CurEstXs(:,:,i)),[0 Mx]);colormap gray;removeTicks;
    gif
end
disp(OutGifFN);





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