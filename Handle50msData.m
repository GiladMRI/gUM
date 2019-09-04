PhaseCyclingBaseP='/autofs/space/daisy_002/users/Gilad/phase_cyclingD/';
addpath(PhaseCyclingBaseP);
addpath([PhaseCyclingBaseP 'utils'])
addpath([PhaseCyclingBaseP 'data']);
addpath([PhaseCyclingBaseP 'linops']);
addpath([PhaseCyclingBaseP 'prox']);
addpath([PhaseCyclingBaseP 'phase_cycling']);
%%
BaseP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/';
FN='meas_MID855_gBP_ASL_SMS_Spi_TI_VD1_ST10_FID51469';
FN='meas_MID857_gBP_ASL_SMS_Spi_TI_VD1_ST11_FID51471';
FN='meas_MID859_gBP_ASL_SMS_Spi_TI_VD1_ST12_FID51473'; % bad?
FN='meas_MID863_gBP_ASL_SMS_Spi_TI_VD1_ST12_FID51477';
FN='meas_MID861_gBP_ASL_SMS_Spi_TI_VD1_ST13_FID51475';
FN='meas_MID865_gBP_ASL_SMS_Spi_TI_VD1_ST14_FID51479';
RefFldMapP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/meas_MID853_BP_fieldmap_9echos_2mm_Full_FID51467/';

BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';

sTwix = mapVBVD([BaseP FN '.dat']);
Data=sTwix.image(:,:,:,:,:,3,:,:,:,:,:,:,:,:,:,:,:);
SData=squeeze(Data);
SData=CombineDims(SData,[4 1]);

sTwix.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal
mkdir([BaseP FN]);
%%
asSlice=sTwix.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[sTwix.hdr.Phoenix.sSliceArray.asSlice{:}];
end

nSlices=numel(asSlice);
for s=1:nSlices
    try
        SlbLoc(1,s)=asSlice(s).sPosition.dSag;
    catch
        SlbLoc(1,s)=0;
    end
    try
        SlbLoc(2,s)=asSlice(s).sPosition.dCor;
    catch
        SlbLoc(2,s)=0;
    end
    try
        SlbLoc(3,s)=asSlice(s).sPosition.dTra;
    catch
        SlbLoc(3,s)=0;
    end
end

RotMat = transpose(Quat2RotMat(sTwix.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;
disp('ok');
%% Load the B0
% sTwixB0=sTwix;
%%
Grads=load('GAll.mat');
CurGrad=Grads.GAll(:,5,1);
%%
g=CurGrad;
FOV_mm=200;
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_ms=10e-3;

k=cumsum([0;g;0])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

kK=k*FOV_mm/1000/2/pi;
figure;
subplot(2,3,1);
plot(kK);
axis square;
axis equal;
ylabel('k');
%%
nAcqPoints=size(SData,1);
SpBW=1e6;
Acq_dT_us=1e6/SpBW;
AcqPointsPerGrad=GradDwellTime_ms*1e9/SpBW;
EstGradDelay_us=9;
EstGradDelay_acqPoints=EstGradDelay_us/Acq_dT_us;

TrajQ=interp1(0:numel(kK)-1,kK,((0:nAcqPoints-1)+AcqPointsPerGrad-EstGradDelay_us)/AcqPointsPerGrad).';
sum(isnan(TrajQ(:)))
BARTTraj=[real(TrajQ) imag(TrajQ) imag(TrajQ)*0].';

disp('Interpolated traj');
figure;plot(grmss(SData,2));hold on;plot(abs(TrajQ)/10000);legend([gmat2cell(num2str((1:12).'),1);{'traj radius'}])
Crossings=quickFindPeaks(diff(abs(TrajQ),2,1).');
setXaxis(Crossings(1)+[-400 400]);
%%
SensB=load([RefFldMapP 'Sens.mat']);
SensB=SensB.SensB;
SnsSzB=gsize(SensB,1:2);

SliOffset=0;
SensX=permute(SensB(:,:,:,SliOffset+(1:nSlices),1),[1 2 3 5 4]);
SensX=gflip(SensX,1:2);
%
FirstEcho=load([RefFldMapP 'FirstEcho.mat']);
FirstEcho=FirstEcho.FirstEcho;

FirstEcho=gflip(FirstEcho(:,:,:,SliOffset+(1:nSlices)),1:2);
Mg=grmss(FirstEcho,3);

B0S=load([RefFldMapP 'B0S.mat']);
B0S=B0S.B0S;

B0Q2=B0S(:,:,SliOffset+(1:nSlices));
B0Q2=gflip(B0Q2,1:2);
disp('ok');
%%
Data2D=CombineDims(SData,[3 1]);

[U,S,sccmtxFull] = svd(Data2D,'econ');
ncc=7;
sccmtx=sccmtxFull(:,1:ncc).';
sccmtxP=permute(sccmtx,[3 2 4 5 6 7 1]);
sccmtxP2=permute(sccmtx,[4 3 2 5 6 7 1]);

DataCC=permute(sum(SData.*sccmtxP,2),[1 7 3:6 2]);

SensCC=permute(sum(SensX.*sccmtxP2,3),[1 2 7 4:6 3]);

disp('CCed');
%%
nReps=size( SData,3);
nChannels=size( sTwix,2);

dx=RotatedLocs(2,1)/sTwix.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;

kx=BARTTraj(1,:)*2*pi;
ky=BARTTraj(2,:)*2*pi;

sNormal=sTwix.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal;
dy=RotatedLocs(1,1)/sTwix.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
modx=double(exp(1i*(dx*kx+dy*ky))');

DataCCX=DataCC.*modx;
disp('ok mod')
%% just pics NUFFT single channel
Nres=96;
Sz=[Nres Nres];
Hnres=Nres/2;
OKPoints=abs(TrajQ)<Hnres;
OKPoints=OKPoints*0+1>0;
DataCCXB=DataCCX(OKPoints,:,:);
BARTTrajB=BARTTraj(:,OKPoints);
DataCCP=permute(DataCCXB,[4 1 5 2 6 3]);
SensM=ones(Nres);
Rec1=bart('pics -t ',BARTTrajB,DataCCP(1,:,1,1,1,5),SensM);

% fgmontage(Mg);MaximizeFig
fgmontage(Rec1);MaximizeFig;
%% just pics NUFFT
SensCCP=permute43(SensCC);
Rec1c=bart('pics -t ',BARTTrajB,DataCCP(1,:,1,:,1,5),SensCCP);
fgmontage(Rec1c);MaximizeFig;
% fgmontage(Mg);MaximizeFig
% fgmontage(B0Q2,[-100 100]);MaximizeFig
title('Fieldmap B0 Hz');removeTicks; colorbar
%%
Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];
%
STraj=TrajQ.';
STraj3=CTo3Rows(STraj);

Sz16=FillOnesTo16(Sz);

SnufftStruct = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
%%
nTS=20;
TSB=GetTSCoeffsByLinear(nAcqPoints,nTS);

Kerns=NUFFT_to_Toep_2blocks(SnufftStruct,TSB);

TimePoints_ms=linspace(0,nAcqPoints*1e3/SpBW,nTS);
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TSC=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3);
% [x y z Channels 1 TS]

ScriptFN_TS=[BaseSP 'nuftTSC_N.txt'];

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # Data is [1 #Traj 1 Ch]
TSBP=permute(TSB,[3 1 4 5 6 7 2]);
TSCP=permute(TSC,[1 2 7 6 5 4 3]);
KernsP=permute(Kerns,[1 2 7 6 5 4 3]);
%%
% NUFTB=bart(['linopScript ' ScriptFN_TS],Sz16,I1,SensCCP,STraj3,TSBP,TSCP,KernsP);
% NUFTB_A=bart(['linopScript -A ' ScriptFN_TS],Sz16,NUFTB,SensCCP,STraj3,TSBP,TSCP,KernsP);
% NUFTB_N=bart(['linopScript -N ' ScriptFN_TS],Sz16,I1,SensCCP,STraj3,TSBP,TSCP,KernsP);
% ShowAbsAngle(NUFTB_A);colorbar
% ShowAbsAngle(NUFTB_N);colorbar
% %%
% ImSize=Sz16;
% DataSize=size(NUFTB);
% x = randn(ImSize) + 1j*randn(ImSize);
% y = randn(DataSize) + 1j*randn(DataSize);
% Ax = bart(['linopScript ' ScriptFN_TS],Sz16,x,SensCCP,STraj3,TSBP,TSCP,KernsP);
% Aty = bart(['linopScript -A ' ScriptFN_TS],Sz16,y,SensCCP,STraj3,TSBP,TSCP,KernsP);
% Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
%% given B0
RecTS=bart(['picsS ' ScriptFN_TS],Sz16,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP,TSCP,KernsP);
fgmontage(RecTS);MaximizeFig;
%%
RecTS_W=bart(['picsS -m -R W:3:0:10.01 ' ScriptFN_TS],Sz16,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP,TSCP,KernsP);
fgmontage(RecTS_W);MaximizeFig;
%% Work on part of the data
CurIdxs=25000+(1:25000);

nTSp=20;
TSBp=GetTSCoeffsByLinear(numel(CurIdxs),nTSp);
STrajp=STraj(:,CurIdxs);
STraj3p=CTo3Rows(STrajp);

SnufftStructp = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STrajp),Sz), Sz, [6 6], Sz*2); % st.om

Kernsp=NUFFT_to_Toep_2blocks(SnufftStructp,TSBp);

TimePoints_msp=linspace(0,numel(CurIdxs)*1e3/SpBW,nTSp);
TimePoints_ms3p=permute(TimePoints_msp,[1 3 2]);
TSCp=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3p);


TSBPp=permute(TSBp,[3 1 4 5 6 7 2]);
TSCPp=permute(TSCp,[1 2 7 6 5 4 3]);
KernsPp=permute(Kernsp,[1 2 7 6 5 4 3]);


RecTS_Wp=bart(['picsS -m -R W:3:0:10.01 ' ScriptFN_TS],Sz16,DataCCP(:,CurIdxs,:,:,:,5),SensCCP,STraj3p,TSBPp,TSCPp,KernsPp);
fgmontage(RecTS_Wp);MaximizeFig;
%%
% RecTS_Wp1=RecTS_Wp;
% RecTS_Wp2=RecTS_Wp;
% %%
% R=RecTS_Wp1./RecTS_Wp2;
%% Now two levels
% L1, e.g. 5
% *                  *                 *                    *
% |           / \              /\                 /         |        with large and decreasing radius
% * * * * * * * * * * * * * * * * * * * * * * * * * * * *      to Level 2, e.g. 50
%                                                                    <- Apply estimated B0
% |               |           |             |               |        Normal TS to traj
% 
% ...........................................................      
nTS_L1=5;
Rad_L12=25; % 5,25 good result
nTS_L2=50;

TSB_L1=GetTSCoeffsByLinearWide(nTS_L2,nTS_L1,Rad_L12);
% figure;plot(TSB_L1)

TSB_L2=GetTSCoeffsByLinear(nAcqPoints,nTS_L2);

Kerns_L2=NUFFT_to_Toep_2blocks(SnufftStruct,TSB_L2);

TimePoints_ms_L2=linspace(0,nAcqPoints*1e3/SpBW,nTS_L2);
TimePoints_ms3_L2=permute(TimePoints_ms_L2,[1 3 2]);
TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
% [x y z Channels 1 TS]

ScriptFN_TS_2L=[BaseSP 'nuftTSC_N_2L.txt'];

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # file 5 is TSB_L1 [1 1 1 1 1 TS_L1 TS_L2] 
% # Data is [1 #Traj 1 Ch]
TSBP_L2=permute(TSB_L2,[3 1 4 5 6 7 2]);
TSCP_L2=permute(TSC_L2,[1 2 7 6 5 4 3]);
KernsP_L2=permute(Kerns_L2,[1 2 7 6 5 4 3]);
TSBP_L1=permute(TSB_L1,[3:7 2 1]);
Sz16_2L=Sz16;
Sz16_2L(6)=nTS_L1;
disp('Prepared 2L');
%%
% RecTS_2L=bart(['picsS -m -R W:3:0:0.1 ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
% fgmontage(RecTS_2L);MaximizeFig;
%% -R K:7:7:.03:HankelizationK:BlkSize:Option:Dim	Hankelized low-rank.
RegCmd='-R K:3:3:1:2:1:0:5';
RegCmd='-R K:3:3:0.001:2:1:0:5';
RecTS_2Lb=bart(['picsS -m ' RegCmd ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2Lb);MaximizeFig;title(RegCmd);
%% Open to L2
RecTS_2LbX=squeeze(sum(RecTS_2Lb.*TSBP_L1,6));
RecTS_2LbY=RecTS_2LbX.*TSC_L2;
fgmontage(RecTS_2LbY);MaximizeFig;title(RegCmd);

RecTS_2LbYx=RecTS_2LbY.*exp(-1i*angle(RecTS_2LbY(:,:,1)));
RecTS_2LbYxx=RecTS_2LbYx.*exp(-1i*angle(TSC_L2));
%%
FRecTS_2LbYxx=fft1cg(RecTS_2LbYxx,3);

figure;plot(grmss(FRecTS_2LbYxx,1:2)./max(grmss(FRecTS_2LbYxx,1:2)))
fgmontage(grmss(FRecTS_2LbYxx(:,:,[1:10 40:end]),3)./grmss(FRecTS_2LbYxx(:,:,[11:39]),3))
%% With warm start
TmpP='/tmp/';
ITSB_L1=permute(pinv(TSB_L1),[3 4 2 5 6 1]);
W1P=sum(RecTS_2LbX1.*ITSB_L1,3);
WarmStartFN=[TmpP 'WarmStart'];
writecfl(WarmStartFN,W1P);

RegCmd='-R K:3:3:1:2:1:0:5';
RecTS_2Lb_w1=bart(['picsS ' ' -W ' WarmStartFN ' -m ' RegCmd ' ' ScriptFN_TS_2L ],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2Lb_w1);MaximizeFig;title(['W1 ' RegCmd]);

% RecTS_2LbX=squeeze(sum(RecTS_2Lb.*TSBP_L1,6));
RecTS_2LbX_w1=squeeze(sum(RecTS_2Lb_w1.*TSBP_L1,6));
RecTS_2LbY_w1=RecTS_2LbX_w1.*TSC_L2;


%% Components:
Mag2D=reshape(abs(RecTS_2LbYx),[],50);
Mag2Dr=Mag2D./max(Mag2D,[],2);

[Um,Sm,Vm]=svd(Mag2D,'econ');

T2svalues_ms=linspace(5,300,40);
Decays=exp(-(1:50)./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

figure;
for i=1:4
    subplot(2,2,i);
    plot(Vm(:,i));hold on
    plot(Vd(:,i),'--');
end
% figure;plot(Vm(:,1:4));
% hold on;
% plot(Vd(:,1:4),'--');
%% Use decays
nComponents_L1=3;
nTS_L2=50;

TSB_L1=Vd(:,1:nComponents_L1);
% figure;plot(TSB_L1)

TSB_L2=GetTSCoeffsByLinear(nAcqPoints,nTS_L2);

Kerns_L2=NUFFT_to_Toep_2blocks(SnufftStruct,TSB_L2);

TimePoints_ms_L2=linspace(0,nAcqPoints*1e3/SpBW,nTS_L2);
TimePoints_ms3_L2=permute(TimePoints_ms_L2,[1 3 2]);
TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
% [x y z Channels 1 TS]

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # file 5 is TSB_L1 [1 1 1 1 1 TS_L1 TS_L2] 
% # Data is [1 #Traj 1 Ch]
TSBP_L2=permute(TSB_L2,[3 1 4 5 6 7 2]);
TSCP_L2=permute(TSC_L2,[1 2 7 6 5 4 3]);
KernsP_L2=permute(Kerns_L2,[1 2 7 6 5 4 3]);
TSBP_L1=permute(TSB_L1,[3:7 2 1]);
Sz16_2L=Sz16;
Sz16_2L(6)=nComponents_L1;
disp('Prepared 2L Components');
%% -R K:7:7:.03:HankelizationK:BlkSize:Option:Dim	Hankelized low-rank.
RegCmd='-R K:3:3:1:2:1:0:5';
RegCmd='-R K:3:3:0.001:2:1:0:5';
RecTS_2LCb=bart(['picsS -m ' RegCmd ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2LCb);MaximizeFig;title(RegCmd);

RegCmdLLR='-b 8 -R L:3:3:0.0001';
RecTS_2LCb_LLR=bart(['picsS -m ' RegCmdLLR ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2LCb_LLR);MaximizeFig;title(RegCmdLLR);

RecTS_2LCbX=squeeze(sum(RecTS_2LCb_LLR.*TSBP_L1,6));
RecTS_2LCbY=RecTS_2LCbX.*TSC_L2;

RecTS_2LCbYx=RecTS_2LCbY.*exp(-1i*angle(RecTS_2LCbY(:,:,1)));
RecTS_2LCbYxx=RecTS_2LCbYx.*exp(-1i*angle(TSC_L2));
FRecTS_2LCbYxx=fft1cg(RecTS_2LCbYxx,3);

figure;plot(grmss(FRecTS_2LCbYxx,1:2)./max(grmss(FRecTS_2LCbYxx,1:2)))
fgmontage(grmss(FRecTS_2LCbYxx(:,:,[1:10 40:end]),3)./grmss(FRecTS_2LCbYxx(:,:,[11:39]),3))
%% On F
ScriptFN_TS_F=[BaseSP 'nuftTSC_N_TS_F.txt'];

nTS_L2=50;

TSB_L2=GetTSCoeffsByLinear(nAcqPoints,nTS_L2);

Kerns_L2=NUFFT_to_Toep_2blocks(SnufftStruct,TSB_L2);

TimePoints_ms_L2=linspace(0,nAcqPoints*1e3/SpBW,nTS_L2);
TimePoints_ms3_L2=permute(TimePoints_ms_L2,[1 3 2]);
TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
% [x y z Channels 1 TS]

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # Data is [1 #Traj 1 Ch]
TSBP_L2=permute(TSB_L2,[3 1 4 5 6 7 2]);
TSCP_L2=permute(TSC_L2,[1 2 7 6 5 4 3]);
KernsP_L2=permute(Kerns_L2,[1 2 7 6 5 4 3]);
Sz16_2L=Sz16;
Sz16_2L(7)=nTS_L2;
MaskOnF=zeros(Sz16_2L);
MaskOnF(:,:,1,1,1,1,20:33)=1;
disp('Prepared TS F');
%%
RegCmd_IF='-R I:0:0.1';
RecTS_F=bart(['picsS -m ' RegCmd_IF ' ' ScriptFN_TS_F],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,MaskOnF);
fgmontage(RecTS_F);MaximizeFig;title(RegCmd_IF);

RecTS_FX=fft1cg(squeeze(RecTS_F),3);
RecTS_2LCbY=RecTS_2LCbX.*TSC_L2;





%%
nTS_L1=11;
Rad_L12=5;
nTS_L2=50;

TSB_L1=GetTSCoeffsByLinearWide(nTS_L2,nTS_L1,Rad_L12);
W1=squeeze(sum(RecTS_2LbX.*permute(TSB_L1,[3 4 1 2]),3));
% figure;plot(TSB_L1)

TSB_L2=GetTSCoeffsByLinear(nAcqPoints,nTS_L2);

Kerns_L2=NUFFT_to_Toep_2blocks(SnufftStruct,TSB_L2);

TimePoints_ms_L2=linspace(0,nAcqPoints*1e3/SpBW,nTS_L2);
TimePoints_ms3_L2=permute(TimePoints_ms_L2,[1 3 2]);
TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
% [x y z Channels 1 TS]

ScriptFN_TS_2L=[BaseSP 'nuftTSC_N_2L.txt'];

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # file 5 is TSB_L1 [1 1 1 1 1 TS_L1 TS_L2] 
% # Data is [1 #Traj 1 Ch]
TSBP_L2=permute(TSB_L2,[3 1 4 5 6 7 2]);
TSCP_L2=permute(TSC_L2,[1 2 7 6 5 4 3]);
KernsP_L2=permute(Kerns_L2,[1 2 7 6 5 4 3]);
TSBP_L1=permute(TSB_L1,[3:7 2 1]);
Sz16_2L=Sz16;
Sz16_2L(6)=nTS_L1;


disp('Prepared 2L');

W1P=permute(W1,[1 2 6 5 4 3]);
WarmStartFN=[BasePA 'WarmStart'];
writecfl(WarmStartFN,W1P);
%% -R K:7:7:.03:HankelizationK:BlkSize:Option:Dim	Hankelized low-rank.
RegCmd='-R K:3:3:1:2:1:0:5';
RecTS_2L_Stg2=bart(['picsS -m -W ' WarmStartFN ' ' RegCmd ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2L_Stg2);MaximizeFig;title(RegCmd);
%%


%% Check that the data is low rank after (maximal) temporal Hankelization
HankelTemporalLen=2;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(nTS_L2,HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[Sz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[Sz, numel(HankelMat)]).*DeHankelizingMatP,3));

[~, s_vals] = llr_thresh_OnPreparedBlocks(H_for(RecTS_2LbY), 0);

fgmontage(s_vals)
%%
[ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(RecTS_2LbY));

% For test: [ U_LLRa, s_LLRa, V_LLRa ] = batch_svd(H_for(repmat(permute(150*(0.8.^(1:50)),[1 3 2]),[Sz 1])));

R1=V_LLR(:,:,2,1)./V_LLR(:,:,1,1); % R1 is simply the decay
fgmontage(angle(R1));
% fgmontage(-log(abs(R1)),[0 0.2])
fgmontage(-1./log(abs(R1)),[0 100]) % TE is ~1ms, so ignore

%% Compare to estimated B0
[ U_LLRa, s_LLRa, V_LLRa ] = batch_svd(H_for(TSC_L2));

R1a=V_LLRa(:,:,2,1)./V_LLRa(:,:,1,1);
fgmontage(angle(R1a));
fgmontage(-log(abs(R1a)),[0 0.2]) % TE is ~1ms, so ignore
fgmontage(-1./log(abs(R1a))) % TE is ~1ms, so ignore

% R1=exp(-TE/T2*)
% T2*=1/(log(R1)/(-TE))=-TE/log(R1)
%% Build TSC_L2 from R1
MR1=abs(R1);
MR1(MR1==0)=1;
MR1=min(MR1,1);
R1c=MR1.*exp(1i*angle(R1));
TSC_L2x=R1c.^(-permute(0:nTS_L2-1,[1 3 2]));
TSC_L2x=exp(1i*angle(TSC_L2x));
% TSC_L2x=TSC_L2x*0+1;
%% Recon single image? or several? same as first, just different B0 est
nTS_L1=5;
Rad_L12=25;
nTS_L2=50;

TSB_L1=GetTSCoeffsByLinearWide(nTS_L2,nTS_L1,Rad_L12);
W1=squeeze(sum(RecTS_2LbX.*permute(TSB_L1,[3 4 1 2]),3));
% figure;plot(TSB_L1)

TSB_L2=GetTSCoeffsByLinear(nAcqPoints,nTS_L2);

Kerns_L2=NUFFT_to_Toep_2blocks(SnufftStruct,TSB_L2);

% TimePoints_ms_L2=linspace(0,nAcqPoints*1e3/SpBW,nTS_L2);
% TimePoints_ms3_L2=permute(TimePoints_ms_L2,[1 3 2]);
% TSC_L2=exp(1i.*2*pi*B0Q2*(1e-3).*TimePoints_ms3_L2);
% [x y z Channels 1 TS]

ScriptFN_TS_2L=[BaseSP 'nuftTSC_N_2L.txt'];

% # file 0 is sensitivity maps [x y z Ch Maps]
% # file 1 is sampling pattern/Trajectory [3 #Traj spokes]
% # file 2 is TSB [1 #traj 1 1 1 1 TS] 
% # file 3 is TSC [x y z 1 1 1 TS]
% # file 4 is Toeplitz kernel [2x 2y z 1 1 1 TS]
% # file 5 is TSB_L1 [1 1 1 1 1 TS_L1 TS_L2] 
% # Data is [1 #Traj 1 Ch]
TSBP_L2=permute(TSB_L2,[3 1 4 5 6 7 2]);
TSCP_L2=permute(TSC_L2x,[1 2 7 6 5 4 3]);
KernsP_L2=permute(Kerns_L2,[1 2 7 6 5 4 3]);
TSBP_L1=permute(TSB_L1,[3:7 2 1]);
Sz16_2L=Sz16;
Sz16_2L(6)=nTS_L1;

disp('Prepared 2L');

W1P=permute(W1,[1 2 6 5 4 3]);
WarmStartFN=[BasePA 'WarmStart'];
writecfl(WarmStartFN,W1P);
%% 
% RegCmd='-R K:3:3:10:2:1:0:5';
RegCmd='-R K:3:3:1:2:1:0:5';
% RegCmd='-R W:3:0:0.01';
RegCmd='-R K:3:3:1:2:1:0:5 -R W:3:0:1 ';
% RecTS_2L_Stg2x=bart(['picsS -m -W ' WarmStartFN ' ' RegCmd ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
RecTS_2L_Stg2x=bart(['picsS -m ' RegCmd ' ' ScriptFN_TS_2L],Sz16_2L,DataCCP(:,:,:,:,:,5),SensCCP,STraj3,TSBP_L2,TSCP_L2,KernsP_L2,TSBP_L1);
fgmontage(RecTS_2L_Stg2x);MaximizeFig;title(RegCmd);

RecTS_2L_Stg2xO=squeeze(sum(RecTS_2L_Stg2x.*TSBP_L1,6));
RecTS_2L_Stg2xO=RecTS_2L_Stg2xO.*TSC_L2x;
fgmontage(RecTS_2L_Stg2xO(:,:,1:6:end))
%%
[ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(RecTS_2L_Stg2xO));
R2=V_LLR(:,:,2,1)./V_LLR(:,:,1,1);
fgmontage(angle(R2));
fgmontage(-log(abs(R2)),[0 0.2])
%% Wavelet smooth B0? T2*
%%


%%
RecTS_2L_2_25=RecTS_2L;
%%
I1=phantom(Nres);
II=bart('linopScript /autofs/space/daisy_002/users/Gilad/gUM/dblszS.txt',FillOnesTo16(Sz),I1);


ScriptFN=[BaseSP 'nuft2D.txt'];
ScriptFN_N=[BaseSP 'nuft2DN.txt'];

Kern=NUFFT_to_Toep_2blocks(SnufftStruct,TSB(:,4));
NUFTB=bart(['linopScript ' ScriptFN],Sz16,I1,STraj3);
NUFTB=NUFTB.*(TSB(:,4).');
NUFTB_NA=bart(['linopScript -A ' ScriptFN],Sz16,NUFTB,STraj3);
ShowAbsAngle((NUFTB_NA))

NUFTB_N=bart(['linopScript -N ' ScriptFN],Sz16,I1,STraj3);
ShowAbsAngle((NUFTB_N))

NUFTB_N2=bart(['linopScript -N ' ScriptFN_N],Sz16,I1,STraj3,Kern);
ShowAbsAngle((NUFTB_N2))
%%