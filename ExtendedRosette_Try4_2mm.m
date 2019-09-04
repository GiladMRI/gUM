addpath('/autofs/space/daisy_002/users/Gilad/gUM/tOptGrad_V0.2/minTimeGradient');
addpath('/autofs/space/daisy_002/users/Gilad/gUM/tOptGrad_V0.2/minTimeGradient/mex interface');
%%
% A  2piKmin=nPetals*2*dk*alpha
% B  (Kmax-Kmin)=dk*nPetals
% A Kmin=nPetals*dk*alpha/pi
% B Kmin=Kmax-dk*nPetals
% nPetals*dk*alpha/pi=Kmax-dk*nPetals
% nPetals*dk*(alpha/pi+1)=Kmax
% nPetals=Kmax/(dk*(alpha/pi+1))
96/(2*(1/pi+1))
%%
Gcm2mTm=10;
GCmms2Tms=10;
radm2cm=1/(2*pi*100);

%%
FOV=192;

res_mm=2;
Nres=FOV/res_mm;
kMax_radm=pi*1000/res_mm;

Kmax=FOV/(2*res_mm);

Fac=(kMax_radm/Kmax)*radm2cm;

dk=2;
alpha=0.5;
dkinner=dk*alpha;

% Kmin=40;
nPetalsT=Kmax/(dk*(alpha/pi+1))
KminT=Kmax-dk*nPetalsT

% nPetals=36;
% nPetals=36;
% Kmin=24;
%
% nPetals=52;
% Kmin=23;

nPetals=ceil(nPetalsT);
Kmin=ceil(KminT);
nPetals=32;
% nPetals=37;

RPK=randperm(nPetals);
RP=randperm(nPetals-1)+1;
BP=(rand(1,nPetals-1)>0.5)*nPetals;

RP=mod(11.*(1:nPetals-1),nPetals-1)+2;
RP=mod(3.*(1:nPetals-1),nPetals-1)+2;
BP=BP*0;
BP=mod(1:nPetals-1,2)*nPetals;
BP=nPetals-BP;
% RP=2:(nPetals);
% RPK=[1:2:nPetals 2:2:nPetals];
% [~,RPK]=sort(RPK);
BP=BP*0;
RP1=RP+BP;
RP2=RP+nPetals-BP;
PetalOrderBase=[1 RP1; RP2 1+nPetals];
for i=2:nPetals
    A=mod(PetalOrderBase(1,i)-PetalOrderBase(1,i-1)+2*nPetals,2*nPetals);
    B=mod(PetalOrderBase(1,i)-PetalOrderBase(1,i-1)+3*nPetals,2*nPetals);
    BP(i-1)=(A>B)*nPetals;
        BP(i-1)=(B>A)*nPetals;
    
    RP1=RP+BP;
    RP2=RP+nPetals-BP;
    PetalOrderBase=[1 RP1; RP2 1+nPetals];
end
%%
KRadAlpha=0.7;
PhiAdd=pi/9;

dd=2;

Phis=linspace(0,pi,nPetals+1);
Phis=Phis(1:end-1);

PhisE=linspace(0,2*pi,nPetals*2+1);
PhisE=PhisE(1:end-1);

RP1=RP+BP;
RP2=RP+nPetals-BP;
PetalOrder=[1 RP1; RP2 1+nPetals];
setdiff(1:(nPetals*2),unique(PetalOrder))
%
Ks=linspace(Kmin,Kmax,nPetals);
KsOrdered=Ks(RPK);
%
% Build petals
% PetalI=3;

% Fac=1;

nRs=20;
nArcs=50;
nConnects=20;

figure;
Traj=[];

subplot(2,2,1);
for PetalI=1:nPetals
    CurK=KsOrdered(PetalI);
    % Ke=KsOrdered(PetalOrder(2,PetalI));
    Ps=PhisE(PetalOrder(1,PetalI));
    Pe=PhisE(PetalOrder(2,PetalI));
    Rs=linspace(0,CurK*KRadAlpha,nRs)*exp(1i*Ps);
    Re=linspace(CurK*KRadAlpha,0,nRs)*exp(1i*Pe);
    Rs=complex(Rs);
    Re=complex(Re);
    ArcPhiStart=Ps+PhiAdd;
    ArcPhiEnd=Pe-PhiAdd;
    if(ArcPhiEnd<ArcPhiStart)
        ArcPhiEnd=ArcPhiEnd+2*pi;
    end
    if(ArcPhiEnd<ArcPhiStart)
        ArcPhiEnd=ArcPhiEnd+2*pi;
    end
    %     ArcPhiEnd=ArcPhiEnd+2*pi;
    %     disp([ArcPhiStart ArcPhiEnd]);
    Arc=CurK*exp(1i*linspace(ArcPhiStart,ArcPhiEnd,nArcs));
    
    
    Arc=CurK*exp(1i*linspace(ArcPhiStart,ArcPhiEnd,ceil(sum(abs(diff(Arc)))/dd)  ));
    dPhiArc=ArcPhiEnd-ArcPhiStart;
    WaveL=pi/10;
    nWave=floor(dPhiArc/WaveL);
    AWave=0*dk*0.25;    
    if(CurK<30)
        Arc=(abs(Arc)+AWave*sin(linspace(0,pi*nWave,numel(Arc)))).*exp(1i*angle(Arc));
    end
    
    
    
    Rs=linspace(0,CurK*KRadAlpha,ceil(sum(abs(diff(Rs)))/dd))*exp(1i*Ps);
    Re=linspace(CurK*KRadAlpha,0,ceil(sum(abs(diff(Re)))/dd))*exp(1i*Pe);
    
    RAWave=0*dk*0.25;
    RWaveL=5;
    RsnWave=floor(CurK*KRadAlpha/RWaveL);
    RenWave=floor(CurK*KRadAlpha/RWaveL);
    Rs=Rs+RAWave*(rand*2-1)*sin(linspace(0,pi*RsnWave,numel(Rs)))*exp(1i*(Ps+pi/4));
    Re=Re+RAWave*(rand*2-1)*sin(linspace(0,pi*RenWave,numel(Re)))*exp(1i*(Pe+pi/4));
    
    Cnct1=conj(PetalConnectf(CurK*KRadAlpha,CurK,PhiAdd,dd))*exp(1i*2*pi/4)*exp(1i*Ps);
    
    
    Cnct2=-(PetalConnectf(CurK*KRadAlpha,CurK,PhiAdd,dd)*exp(1i*2*pi/4)*exp(1i*Pe));
    Cnct2=fliplr(Cnct2);
    if(isempty(Arc))
        Cnct2=Cnct2(2:end);
    end
    
    Rs=Rs(1:end-1);
    if(PetalI>1)
        Rs=Rs(2:end);
    end
    Re=Re(2:end);
    Arc=Arc(2:end-1);
    nArc(PetalI)=numel(Arc);
    Traj=[Traj Rs Cnct1 Arc Cnct2 Re];
    
    plot(real(Rs),imag(Rs),'k.-');
    hold on;
    plot(Re,'k.-');
    plot(Arc,'k.-');
    
    plot(Cnct1,'g.-');
    plot(Cnct2,'g.-');
end
axis square
axis equal
%


Gmax_mTm=38;
Smax_Tms=155;

Gmax_GCcm=Gmax_mTm/Gcm2mTm;
Smax_GCmms=Smax_Tms/GCmms2Tms;

Traj1cm=Traj*Fac;
DwellTimeGrad_ms=10e-3;

Traj2D=[real(Traj1cm);imag(Traj1cm)].';
%
[C,time,g,s,k, phi, sta, stb] = minTimeGradient(Traj2D,0, 0, 0, Gmax_GCcm, Smax_GCmms,DwellTimeGrad_ms);

TrajC=C(:,1)+1i*C(:,2);
%
% figure;plot(k(:,1),k(:,2))
sC=s(:,1)+1i*s(:,2);
gC=g(:,1)+1i*g(:,2);

kC=k(:,1)+1i*k(:,2);
kC=kC/Fac;

SOverflowIdx=find(abs(sC)>Smax_GCmms*1.1);
GOverflowIdx=find(abs(gC)>Gmax_GCcm*1.1);
hold on;
% plot(k(SOverflowIdx,1),k(SOverflowIdx,2),'ro','LineWidth',4)
% plot(k(GOverflowIdx,1),k(GOverflowIdx,2),'ko','LineWidth',4)

plot(kC(SOverflowIdx),'ro','LineWidth',4)
plot(kC(GOverflowIdx),'ko','LineWidth',4)

title(time)
xlabel(['#Petals: ' num2str(nPetals) ', Kmin=' num2str(Kmin) ', dkinner=' num2str(dkinner,'%.1f')]);
axis([-1 1 -1 1]*Kmax*1.1);
subplot(2,2,3);
plot(gC);
xlabel(GOverflowIdx);
axis([-1 1 -1 1]*Gmax_GCcm);
axis square
axis equal
subplot(2,2,4);
plot(sC);
xlabel(SOverflowIdx);
axis([-1 1 -1 1]*Smax_GCmms);
axis square
axis equal
%%
% if false
    dk=2;
    nLoops=Kmax/dk;
    VD=1;
    NPnts=2000;
    Traj=(linspace(0,1,NPnts).^VD)*Kmax.*exp(1i*linspace(0,2*pi*nLoops,NPnts));
    figure;
    subplot(2,2,1);
    plot(Traj);
    Traj1cm=Traj*Fac;
    DwellTimeGrad_ms=10e-3;
    
    Traj2D=[real(Traj1cm);imag(Traj1cm)].';
    %
    [C,time,g,s,k, phi, sta, stb] = minTimeGradient(Traj2D,0, 0, 0, Gmax_GCcm, Smax_GCmms,DwellTimeGrad_ms);
    
    TrajC=C(:,1)+1i*C(:,2);
    %
    % figure;plot(k(:,1),k(:,2))
    sC=s(:,1)+1i*s(:,2);
    gC=g(:,1)+1i*g(:,2);
    
    kC=k(:,1)+1i*k(:,2);
    kC=kC/Fac;
    
    SOverflowIdx=find(abs(sC)>Smax_GCmms*1.1);
    GOverflowIdx=find(abs(gC)>Gmax_GCcm*1.1);
    hold on;
    % plot(k(SOverflowIdx,1),k(SOverflowIdx,2),'ro','LineWidth',4)
    % plot(k(GOverflowIdx,1),k(GOverflowIdx,2),'ko','LineWidth',4)
    
    plot(kC(SOverflowIdx),'ro','LineWidth',4)
    plot(kC(GOverflowIdx),'ko','LineWidth',4)
    
    title(time)
    xlabel(['#loops: ' num2str(nLoops) ', VD=' num2str(VD,'%.1f')]);
    axis([-1 1 -1 1]*Kmax*1.1);
    subplot(2,2,3);
    plot(gC);
    xlabel(GOverflowIdx);
    axis([-1 1 -1 1]*Gmax_GCcm);
    axis square
    axis equal
    subplot(2,2,4);
    plot(sC);
    xlabel(SOverflowIdx);
    axis([-1 1 -1 1]*Smax_GCmms);
    axis square
    axis equal
% end
%% Take 50-echo data
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';
AcqParams=load([BaseP 'meas_prot_Fully_Sample.mat']);
EchoTime=0.93; % ms
nEchos=50;
TotalAcqTime=nEchos*0.93; % ms
nTrajTrg=TotalAcqTime*1000/2.5;
%%
SliIdx=26;
Sz=[Nres Nres];

load([BaseP 'Sli' num2str(SliIdx) '.mat']);
CurSli192=padarray(padarray(CurSliS(:,1:184,:,:),[4,8,0,0],'pre'),[17,0,0,0],'post');

CurSli128=imresizeBySlices(CurSli192,Sz);
%%
clear Sens SensS
load([BaseP 'Sli' num2str(SliIdx) '_Sens.mat']);
Sens=double(SensS);
Sens=permute43(Sens);
Sensr=imresizeBySlices(Sens,Sz);

nChannels=size(Sens,4);
%%
% CurSliCombined=sum(CurSli192.*conj(Sens),4);
CurSliCombinedr=sum(CurSli128.*conj(Sensr),4);
GT=CurSliCombinedr(:,:,1);
%%
Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];
%%
nufftStructR = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(kC.'),Sz), Sz, [6 6], Sz*2, Sz/2); % st.om
SN=nufftStructR.sn;

PaddedWithSensSN=padarray(CurSli128.*SN,[Sz 0 0],'post');
FPaddedWithSens=fft2(PaddedWithSensSN);
RFPaddedWithSens=double(reshape(FPaddedWithSens,[prod(Sz)*4 50 12]));

PaddedWithSensSN_MagOnly=padarray(abs(CurSliCombinedr).*Sensr.*SN,[Sz 0 0],'post');
RFPaddedWithSens_MagOnly=double(reshape(fft2(PaddedWithSensSN_MagOnly),[prod(Sz)*4 50 12]));
%%
kC=interp1(1:numel(kC),kC,1:0.25:numel(kC)-0.001).';
nTraj=numel(kC);
HnTraj=floor(nTraj/2);

nufftStructR = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(kC.'),Sz), Sz, [6 6], Sz*2, Sz/2); % st.om

DataFullTrajAllEchos=reshape(nufftStructR.p*reshape(RFPaddedWithSens,[prod(Sz)*4 nEchos*nChannels]),[],nEchos,nChannels);
TSB=GetTSCoeffsByLinear(nTraj,nEchos);
DataWithDecay=sum(DataFullTrajAllEchos.*TSB,2);

DataFullTrajAllEchos_MagOnly=reshape(nufftStructR.p*reshape(RFPaddedWithSens_MagOnly,[prod(Sz)*4 nEchos*nChannels]),[],nEchos,nChannels);
DataWithDecay_MagOnly=sum(DataFullTrajAllEchos_MagOnly.*TSB,2);
%% Rec "first echo"
DataFirstEcho=squeeze(DataFullTrajAllEchos(:,1,:));
DataFirstEcho=permute(DataFirstEcho,[4 1 3 2])/sqrt(prod(Sz));

RecrFirstEchoh1=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC(1:HnTraj,1).'),DataFirstEcho(1,1:HnTraj,1,:),Sensr);
RecrFirstEchoh2=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC(HnTraj:end,1).'),DataFirstEcho(1,HnTraj:end,1,:),Sensr);
RecrFirstEchoF=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC(1:end,1).'),DataFirstEcho(1,1:end,1,:),Sensr);

fgmontage(cat(3,RecrFirstEchoh1,RecrFirstEchoh2,RecrFirstEchoF));
RecrN=Recr*grmss(GT)/grmss(Recr);
ssimScr=ssim(abs(RecrN),abs(GT));
ShowAbsAngle(Recr);subplot(1,2,1);title(ssimScr);
figure;imshowpair(abs(RecrN),abs(GT),'checkerboard');title(ssimScr);
%% Rec B0-agnostic
DataWithDecayP=permute(DataWithDecay,[4 1 2 3])/sqrt(prod(Sz));

RecrB0Agn=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC(1:end,1).'),DataWithDecayP,Sensr);

DataWithDecay_MagOnlyP=permute(DataWithDecay_MagOnly,[4 1 2 3])/sqrt(prod(Sz));
RecrB0Agn_MagOnly=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC(1:end,1).'),DataWithDecay_MagOnlyP,Sensr);
%% Rec LLR on Magnitue
Combined2D=reshape(CurSliCombinedr,[],50);
[U,S,V]=svd(abs(Combined2D),'econ');
dS=diag(S);


nTS=4;
TSB_Echos=V(:,1:nTS);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);
clear TSBk
for i=1:4
    TSBk(:,i)=interp1(1:nEchos,V(:,i),linspace(1,nEchos,nTraj));
end
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
%% Rec LLR on Mag only using components
ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_Script.txt';
RecrDecayMagLLRC=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayMagLLRCT=sum(RecrDecayMagLLRC.*TSB_EchosP,7);
%% with Teoplitz - Not ready yet. For general components (as here), need #components^2 thing
ScriptFN_Toep='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_ScriptN.txt';

TSKern=NUFFT_to_Toep_2blocks(nufftStructR, TSBk);
TSKernP=permute(TSKern,[1 2 7 6 5 4 3]);

RecrDecayMagLLRC_Toep=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN_Toep],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP,TSKernP);
RecrDecayMagLLRCT_Toep=sum(RecrDecayMagLLRC_Toep.*TSB_EchosP,7);
%%
AA=bart(['linopScript ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],rand([Sz 1 1 1 1 nTS ones(1,9)]),Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
AAAdj=bart(['linopScript -A ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],AA,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
AAN=bart(['linopScript -N ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],rand([Sz 1 1 1 1 nTS ones(1,9)]),Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);

AA2=bart(['linopScript ' ScriptFN_Toep],[Sz 1 1 1 1 nTS ones(1,9)],rand([Sz 1 1 1 1 nTS ones(1,9)]),Sensr,CTo3Rows(kC(1:end,1).'),TSBkP,TSKernP);
%% Rec LLR on Mag only using TS
ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_Script.txt';
nTS=8; % Good result
% nTS=50;
% TSBk=GetTSCoeffsByLinear(nTraj,nTS);
TSBk=GetTSCoeffsByLinearWide(nTraj,nTS,nTraj);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

RecrDecayMagLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayMagLLRT=sum(RecrDecayMagLLR.*TSB_EchosP,7);
%% Rec LLR using TS with B0
nTS=15;
RadFac=10;
% TSBk=GetTSCoeffsByLinear(nTraj,nTS);
TSBk=GetTSCoeffsByLinearWide(nTraj,nTS,nTraj/RadFac);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
% TSB_Echos=GetTSCoeffsByLinear(nEchos,nTS);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos/RadFac);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

RecrDecayLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayLLRT=sum(RecrDecayLLR.*TSB_EchosP,7);
%% With warmstart
WarmstartFN=['/autofs/space/daisy_002/users/Gilad/gUM/' 'Warmstart'];
WarmstartFN2=['/autofs/space/daisy_002/users/Gilad/gUM/' 'Warmstart2'];

Warmstart=RecrDecayLLR;
writecfl(WarmstartFN,Warmstart);

RecrDecayLLR2=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ' -W ' WarmstartFN ' ' ScriptFN ],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);

Warmstart2=RecrDecayLLR2;
writecfl(WarmstartFN2,Warmstart2);

TSBk=GetTSCoeffsByLinear(nTraj,nTS);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinear(nEchos,nTS);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

RecrDecayLLRx=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ' -W ' WarmstartFN2 ' ' ScriptFN ],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayLLRT=sum(RecrDecayLLR.*TSB_EchosP,7);

%% Rec using LR on TH
RecrDecayLRH=bart(['picsS -m -u 0.010 -R K:67:67:0.01:3:1:2:6 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayLRHT=sum(RecrDecayLRH.*TSB_EchosP,7);

Warmstart=RecrDecayLRH;
writecfl(WarmstartFN,Warmstart);

RecrDecayLRH2=bart(['picsS -m -u 0.010 -R K:67:67:0.01:3:1:2:6 ' ' -W ' WarmstartFN ' ' ScriptFN ],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayLRHT2=sum(RecrDecayLRH.*TSB_EchosP,7);

%%
RecrDecayLRH=bart(['picsS -m -u 0.010 -R K:67:67:0.01:3:1:1:6 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayLRHT2=sum(RecrDecayLRH2.*TSB_EchosP,7);
%% Open first
ScriptFNOpenFirst='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_Scriptx.txt';

TSBk=GetTSCoeffsByLinearWide(nTraj,nTS,nTraj);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

RecrDecayMagLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);

TSBk=GetTSCoeffsByLinearWide(nTraj,nTS,nTraj);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

% TSB_Echos=V(:,1:nTS);
TSB_EchosP2=permute(TSB_Echos,[7 3 2 4 5 6 1]);
TSB_EchosP3=permute(TSB_Echos,[3 4 2 1]);

TSBk=GetTSCoeffsByLinear(nTraj,nEchos);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);

% Sensr_PEst=Sensr.*exp(1i*angle(permute(CurSliCombinedr(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));
Sensr_PEst=Sensr.*exp(1i*angle(permute(RecrDecayLRHT(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));

RecrDecayOPOPLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFNOpenFirst],[Sz nTS 1 1 1 1 ones(1,9)],DataWithDecayP,Sensr_PEst,CTo3Rows(kC(1:end,1).'),TSBkP,TSB_EchosP2);
RecrDecayOPOPLLRT=squeeze(sum(RecrDecayOPOPLLR.*TSB_EchosP3,3));
%%
RecrDecayOPOPLLRT_GivenRealPhase=RecrDecayOPOPLLRT;
ShowAbsAngle(RecrDecayOPOPLLRT_GivenRealPhase(:,:,1:15:end));title('RecrDecayOPOPLLRT_GivenRealPhase');
%%
RecrDecayOPOPLLRT_GivenLRTHPhase=RecrDecayOPOPLLRT;
ShowAbsAngle(RecrDecayOPOPLLRT_GivenLRTHPhase(:,:,1:15:end));title('RecrDecayOPOPLLRT_GivenLRTHPhase');
%% Wavelet denoising
IdentScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/IdentS.txt';
QQ=bart(['picsS -m -R W:3:0:0.03 ' IdentScriptFN],[Sz ones(1,14)],angle(RecrDecayLRHT(:,:,20)));
%%
RecrDecayLRH2end=RecrDecayLRHT2(:,:,2:end);
RecrDecayLRH1endm1=RecrDecayLRHT2(:,:,1:end-1);

RecrDecayLRH2end=RecrDecayLRHT2(:,:,:,:,:,:,2:end);
RecrDecayLRH1endm1=RecrDecayLRHT2(:,:,:,:,:,:,1:end-1);

FmacScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/fmacS.txt';
QQ=bart(['picsS -m -R W:3:0:0.01 ' FmacScriptFN],[Sz ones(1,14)],RecrDecayLRH2end,RecrDecayLRH1endm1);
QQ=bart(['picsS -m -R T:3:0:9 ' FmacScriptFN],[Sz ones(1,14)],RecrDecayLRH2end,RecrDecayLRH1endm1);
QQ=bart(['picsS -m -R Q:9 ' FmacScriptFN],[Sz ones(1,14)],RecrDecayLRH2end,RecrDecayLRH1endm1);
JJJ=(QQ.^permute(0:49,[1 3 2]));

PRemoved=CurSliCombinedr.*exp(-1i*angle(JJJ));
PRemoved2=PRemoved./PRemoved(:,:,1);
fgmontage(angle(PRemoved2(:,:,1:9:end)));title('Residual phase');
%% Open first
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);


% TSB_Echos=V(:,1:nTS);
% TSB_Echos=GetTSCoeffsByLinear(nEchos,nTS);
TSB_EchosP2=permute(TSB_Echos,[7 3 2 4 5 6 1]);
TSB_EchosP3=permute(TSB_Echos,[3 4 2 1]);

TSBk=GetTSCoeffsByLinear(nTraj,nEchos);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);

% Sensr_PEst=Sensr.*exp(1i*angle(permute(CurSliCombinedr(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));
Sensr_PEst=Sensr.*exp(1i*angle(permute(JJJ(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));

RecrDecayOPOPLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFNOpenFirst],[Sz nTS 1 1 1 1 ones(1,9)],DataWithDecayP,Sensr_PEst,CTo3Rows(kC(1:end,1).'),TSBkP,TSB_EchosP2);
RecrDecayOPOPLLRT=squeeze(sum(RecrDecayOPOPLLR.*TSB_EchosP3,3));

RecrDecayOPOPLLRTd1=RecrDecayOPOPLLRT./RecrDecayOPOPLLRT(:,:,1);

fgmontage(angle(RecrDecayOPOPLLRTd1(:,:,1:9:end)));title('Residual phase est');
%%
ForP=RecrDecayOPOPLLRT.*exp(1i*angle(JJJ));
RecrDecay2end=ForP(:,:,2:end);
RecrDecay1endm1=ForP(:,:,1:end-1);

FmacScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/fmacS.txt';
QQ=bart(['picsS -m -u 1000000 -R W:3:0:1 ' FmacScriptFN],[Sz ones(1,14)],RecrDecay2end,RecrDecay1endm1);
% QQ=bart(['picsS -m  -u 10000000 -R T:3:3:1000 ' FmacScriptFN],[Sz ones(1,14)],RecrDecay2end,RecrDecay1endm1);
% QQ=bart(['picsS -m -R T:3:0:9 ' FmacScriptFN],[Sz ones(1,14)],RecrDecayLRH2end,RecrDecayLRH1endm1);
% QQ=bart(['picsS -m -R Q:9 ' FmacScriptFN],[Sz ones(1,14)],RecrDecayLRH2end,RecrDecayLRH1endm1);
JJJ2=(QQ.^permute(0:49,[1 3 2]));

PRemoved=CurSliCombinedr.*exp(-1i*angle(JJJ2));
PRemoved2=PRemoved./PRemoved(:,:,1);
fgmontage(angle(PRemoved2(:,:,1:9:end)));title('Residual phase2');
%% Open first
TSB_EchosBase=GetTSCoeffsByLinear(nEchos,nTS);
TSB_EchosBaseP=permute(TSB_EchosBase,[7 3 1 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);
TSB_EchosP2=permute(TSB_Echos,[7 3 2 4 5 6 1]);
TSB_EchosP3=permute(TSB_Echos,[3 4 2 1]);

TSBk=GetTSCoeffsByLinear(nTraj,nEchos);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);

% Sensr_PEst=Sensr.*exp(1i*angle(permute(CurSliCombinedr(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));
Sensr_PEst=Sensr.*exp(1i*angle(permute(JJJ2(:,:,linspace(1,nEchos,nEchos)),[1 2 7 6 5 4 3])));

ITSB_EchosBase=pinv(TSB_EchosBase);
ITSB_EchosBaseP=permute(ITSB_EchosBase,[7 3 2 4 5 6 1]);


Warmstart=squeeze(sum(RecrDecayOPOPLLRT.*exp(1i*angle(JJJ)).*exp(-1i*angle(JJJ2)).*ITSB_EchosBaseP,3));

Warmstart=squeeze(sum(RecrDecayLRHT2.*exp(-1i*angle(JJJ2)).*ITSB_EchosP,3));
writecfl(WarmstartFN,Warmstart);

RecrDecayOPOPLLR2=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 -W ' WarmstartFN ' ' ScriptFNOpenFirst],[Sz nTS 1 1 1 1 ones(1,9)],DataWithDecayP,Sensr_PEst,CTo3Rows(kC(1:end,1).'),TSBkP,TSB_EchosP2);
RecrDecayOPOPLLRT2=squeeze(sum(RecrDecayOPOPLLR2.*TSB_EchosP3,3));

RecrDecayOPOPLLRT2d1=RecrDecayOPOPLLRT2./RecrDecayOPOPLLRT2(:,:,1);

fgmontage(angle(RecrDecayOPOPLLRTd1(:,:,1:9:end)));title('Residual phase est');
%%















%%

SRecrDecayLRH=squeeze(RecrDecayLRH);

TT=CurSliCombinedr(:,:,4:3:end)./CurSliCombinedr(:,:,1:3:end-3);
TTT=sum(abs(CurSliCombinedr(:,:,4:3:end)).*exp(1i*angle(TT)),3);

SRecrDecayLRHT=squeeze(RecrDecayLRHT);
% HRecrDecayLRHT=
RecH=H_for(squeeze(RecrDecayLRHT));
%%
Sensr_PEst=Sensr.*exp(1i*angle(RecrDecayLRH));
Sensr_PEst=Sensr.*exp(1i*angle(permute(CurSliCombinedr(:,:,linspace(1,nEchos,nTS)),[1 2 7 6 5 4 3])));
% Sensr_PEst=Sensr.*exp(1i*angle(permute(CurSliCombinedr(:,:,linspace(1,50,50)),[1 2 7 6 5 4 3])));
RadFac=8;
TSBk=GetTSCoeffsByLinearWide(nTraj,nTS,nTraj/RadFac);
TSBkP=permute(TSBk,[7 1 3 4 5 6 2]);
TSB_Echos=GetTSCoeffsByLinearWide(nEchos,nTS,nEchos/RadFac);
TSB_EchosP=permute(TSB_Echos,[7 3 1 4 5 6 2]);

RecrDecayPEstLLR=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecayP,Sensr_PEst,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayPEstLLRT=sum(RecrDecayPEstLLR.*TSB_EchosP,7);
%%
ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_Scriptx.txt';

I0=zeros([Sz 1 1 1 1 nEchos]);
ZZ=bart(['linopScript ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],I0,Sensr_PEst,CTo3Rows(kC(1:end,1).'),TSBkP);

%% Hankelize
WhichEchosToUse=1:nEchos;
WhichEchosToUse=15:45;
WhichEchosToUse=5:3:45;
nEchosToUse=numel(WhichEchosToUse);
nEchosToUse=nTS;
HankelTemporalLen=2;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(nEchosToUse,HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[Sz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[Sz, numel(HankelMat)]).*DeHankelizingMatP,3));

tmp=SmoothBySlices(RecrDecayLRHT(:,:,15:40),[20 20],3);
tmp=RecrDecayLRHT(:,:,WhichEchosToUse);
tmp=exp(1i.*angle(tmp));
[ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(tmp));
xx=V_LLR(:,:,2,1)./V_LLR(:,:,1,1);

[ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(squeeze(RecrDecayLRHT(:,:,WhichEchosToUse))));
xx=V_LLR(:,:,2,1)./V_LLR(:,:,1,1);

[AA, s_vals] = llr_thresh_OnPreparedBlocks(H_for(squeeze(RecrDecayLLRT)), 10000);
BB=H_inv(AA);
fgmontage(s_vals,'Size',[1 2])


%%
I0=zeros([Sz 1 1 1 1 nTS]);
ZZ=bart(['linopScript ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],I0,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);

%%
% FMSensr=bart('fftmod 7 ',Sensr);
RecrB0AgnS=bart(['picsS -m -R W:3:0:0.00001 ' ScriptFN],[Sz 1 1 1 1 1 ones(1,9)],DataWithDecayP,Sensr,CTo3Rows(kC(1:end,1).'),ones(1,nTraj));

nTraj=100^2*5;
SpiC=linspace(0,80,nTraj).*exp(1i*linspace(0,2*pi*80,nTraj));
BTraj=[real(SpiC); imag(SpiC); imag(SpiC)*0];


Data=bart('nufft ',BTraj,AAA);
ScriptFN=[BaseFolder filesep 'Snuft_v1.txt'];

% QQ=bart(['picsS ' ScriptFN],[100 100 ones(1,14)],Data,SAAA,BTraj);
% ShowAbsAngle(QQ)

QQ=bart(['picsS -m -R T:3:3:100 ' ScriptFN],[100 100 ones(1,14)],Data,SAAA,BTraj);

%% Now per timepoint
nTraj=numel(kC);
PerEcho=floor(nTraj/nEchos);
nTrajM=PerEcho*nEchos;
tmp=linspace(1,nTrajM+1,nEchos+1);
StartLocs=floor(tmp(1:end-1));
EndLocs=[StartLocs(2:end)-1 nTrajM];

PerEcho=1000;
StartLocs=randi(nTraj-PerEcho,1,nEchos);
EndLocs=StartLocs+PerEcho-1;
DataT=zeros(PerEcho,12,nEchos);
for i=1:nEchos
    disp(i);
    nufftStructRT{i} = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(kC(StartLocs(i):EndLocs(i),1).'),Sz), Sz, [6 6], Sz*2, Sz/2); % st.om
    DataT(:,:,i)=nufftStructRT{i}.p*squeeze(RFPaddedWithSens(:,i,:));
end
clear TrajT
for i=1:nEchos
    TrajT(:,:,i)=CTo3Rows(kC(StartLocs(i):EndLocs(i),1).');
end

Data=permute(DataT,[4 1 3 2])/sqrt(prod(Sz));
%%
kC=Rows2Complex(Traj).';
%%
nEchosToUse=4;
SelectedEchos=floor(linspace(1,nEchos,nEchosToUse));
FullDataSelectedEchos=RFPaddedWithSens(:,SelectedEchos,:);

PerEcho=4500;
PerEcho=floor(numel(kC)/nEchosToUse);
StartLocs=randi(nTraj-PerEcho,1,nEchosToUse);
StartLocs=floor(linspace(1,nTraj-PerEcho-1,nEchosToUse));
EndLocs=StartLocs+PerEcho-1;
DataT=zeros(PerEcho,12,nEchosToUse);
for i=1:nEchosToUse
    disp(i);
    nufftStructRT{i} = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(kC(StartLocs(i):EndLocs(i),1).'),Sz), Sz, [6 6], Sz*2, Sz/2); % st.om
    DataT(:,:,i)=nufftStructRT{i}.p*squeeze(FullDataSelectedEchos(:,i,:));
end
clear TrajT
for i=1:nEchosToUse
    TrajT(:,:,i)=CTo3Rows(kC(StartLocs(i):EndLocs(i),1).');
end

Data=permute(DataT,[4 1 3 2])/sqrt(prod(Sz));
DataP=permute(Data,[1 2 5 4 8 7 6 3]);
TrajTP=permute(TrajT,[1 2 8 7 6 5 4 3]);
%%
SensWithP=Sensr.*PGT(:,:,1,1,1,1,1,SelectedEchos);
RecTr=bart('pics -m -u 0.1 -R K:131:131:0.001:3:1:2:7 -t ',TrajTP,DataP,SensWithP);

RecTrNoP=bart('pics -m -u 0.1 -R K:131:131:0.001:3:1:2:7 -t ',TrajTP,DataP,Sensr);

RecTrNoP3=bart('pics -m -u 0.03 -R K:131:131:0.001:3:1:2:7 -t ',TrajTP,DataP,Sensr);

RecTrNoP2=bart('pics -m -u 0.3 -R K:131:131:0.001:3:1:2:7 -t ',TrajTP,DataP,Sensr);

%%
% RecTr=bart('pics -m -R W:3:0:0.00001 -t ',CTo3Rows(kC.'),Data,Sensr);

RecTr1=bart('pics -m -R W:3:0:0.00001 -t ',TrajT,Data,Sensr);
ShowAbsAngle(RecTr);subplot(1,2,1);title('RecTr')


RecTrW=bart('pics -m -R W:3:0:0.00001 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);

RecTr=bart('pics -m -R L:3:3:0.01 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);

RecTr=bart('pics -m -u 1000 -R K:131:131:0.01:3:1:1:7 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);

RecTr=bart('pics -m -u 0.01 -R K:131:131:0.01:3:1:2:7 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);

RecTr=bart('pics -m -u 1 -R K:131:131:0.01:3:1:2:7 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);

RecTr=bart('pics -m -u 1 -R K:131:131:0.01:3:3:2:7 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr);
%%
RecTr=bart('pics -m -u 0.03 -R K:131:131:0.001:3:1:2:7 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr.*PGT);

RecTr=bart('pics -b 1 -m -R L:3:3:100 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr.*PGT);
RecTr=bart('pics -m -R T:128:128:100 -t ',permute(TrajT,[1 2 8 7 6 5 4 3]),permute(Data,[1 2 5 4 8 7 6 3]),Sensr.*PGT);
%%
PGT=permute(exp(1i*angle(CurSliCombinedr)),[1 2 8 7 6 5 4 3]);
%%
HankelTemporalLen=2;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(nEchosToUse,HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[Sz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[Sz, numel(HankelMat)]).*DeHankelizingMatP,3));

[AA, s_vals] = llr_thresh_OnPreparedBlocks(H_for(squeeze(RecTrNoP)), 10000);
BB=H_inv(AA);
fgmontage(s_vals,'Size',[1 3])

%%
% 'falsecolor' 'blend' | 'diff' | 'montage' 'checkerboard'
%% Apply NUFT
% I1=phantom(Nres);
% Data=bart('nufft ',CTo3Rows(kC.'),I1);

Data=bart('nufft ',CTo3Rows(kC.'),CurSli192(:,:,1,:));
%%

%% Call BART for recon
Rec=bart('pics -R Q:1 -t ',CTo3Rows(kC.'),Data,Sens);

Rec=bart('pics -R Q:1 -t ',CTo3Rows(kC.'),Data,ones(Sz));
%%
Rec=bart('pics -m -R W:3:0:0.001 -t ',CTo3Rows(kC.'),Data,Sens);
Recr=bart('pics -m -R W:3:0:0.001 -t ',CTo3Rows(kC.'),Data,Sensr);
% 
% Rec=bart('pics -R W:3:0:0.001 -t ',CTo3Rows(kC.'),Data,ones(Sz));
%%
% function [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C,rv, g0, gfin, gmax, smax,T, ds,show)

% function [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C,rv, g0, gfin, gmax, smax,T, ds,show)

%   [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C, RIV/RV, g0, gfin, gmax, smax,T, ds, show)
%   [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C, [], [], [], gmax, smax, T, [], [])
%
%   Given a k-space trajectory C = [x y z], gradient and slew constraints,
%   This function will return a new parameterization that will meet these constraints
%   while getting from one point to the other in minimum time.
%   The constraints can be either magnitude (|g| < Gmax and |s| < Smax) or as a set of
%   individual for each direction (|gx|, |gy|, |gz| < Gmax and |sx|, |sy|, |sz| < Smax).
%   It will thus either call the function minTimeGradientRIV.m or
%   minTimeGradientRV.m based on the specified solution type.
%
%   Input Values :
%
%   C       -   The Curve in k-space given in any parametrization [1/cm]
%               C should be inputed as an Nx2 or Nx3 trajectory.
%   RIV/RV   -  0 for rotationally invariant solution (magnitude constraints),
%                1 for rotationally variant solution (individual gradient/slew constraints).
%                Leave empty for default (0 for the rotationally invariant solution).
%   g0      -   Initial gradient amplitude (leave empty for g0 = 0)
%   gfin    -   Gradient value at the end of the trajectory. If not possible,
%               the result would be the largest possible ampltude.
%               (Leave empty if you don't care to get maximum gradient.)
%   gmax    -   Maximum gradient [G/cm] (3.9 default)
%   smax    -   Maximum slew [G/Cm/ms]  (14.5 default)
%   T       -   Sampling time interval [ms] (4e-3 default)
%   ds      -   step size for ODE integration, leave empty to use default value
%   show    -   Show plots while optimizing (Warning: This will make the
%               process considerably slower!)
%
% return values:
%   C       - reparametrized curve, sampled at T[ms]
%   time    - total time to get to the end
%   g       - gradiet waveform [G/cm]
%   s       - slew rate [G/cm/ms]
%   k       - exact k-space corresponding to gradient g (This function reparametrizes
%             C, then takes a derivative. Numerical errors in the derivative can lead to
%             deviation.
%   phi     - Geometry constraints on the amplitude vs. arclength
%   sta     - Solution for the forward ODE
%   stb     - Solution for the backward ODE
