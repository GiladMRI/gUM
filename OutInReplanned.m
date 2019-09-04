NN=128*128;
FOVx=192;
dFOV=FOVx/1000;
paramLongROSamples=1024*18;
nInnerInterleaves=20;
% AccR=18.23;
AccR=NN/(paramLongROSamples/nInnerInterleaves);
% AccR=9.2;
% AccR=4.3;
% AccR=17.2; % 6
AccR=22;
GRAD_RASTER_TIME=10;

spBW=4e5;

TrajPointsPerGrad=GRAD_RASTER_TIME*spBW/1e6;
TimePerAcqPoint_us=1e6/spBW;
paramLongInterleaves=1;
VD=1;
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;

ExtraDelay=0;

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/nInnerInterleaves,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,nInnerInterleaves,1,ExtraDelay]);

dkToGFac=FOVx/1000/2/pi;
BARTTrajx=kTraj.'*FOVx/1000/2/pi;
nTrajG=size(BARTTrajx,2);
C=BARTTrajx(1,1:size(BARTTrajx,2)/2)+1i*BARTTrajx(2,1:size(BARTTrajx,2)/2);

MaxK=max(BARTTrajx(:));
Acc=ceil(MaxK*2).^2/nTrajG;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajxQ=kTrajQ.'*FOVx/1000/2/pi;
nTraj=size(BARTTrajxQ,2);
Traj=BARTTrajxQ;
MaxK
kTrajO=kTraj;
%% Now pass trajectory by minTimeGradient
Gcm2mTm=10;
GCmms2Tms=10;
radm2cm=1/(2*pi*100);

Gmax_mTm=38;
Smax_Tms=155;

Gmax_GCcm=Gmax_mTm/Gcm2mTm;
Smax_GCmms=Smax_Tms/GCmms2Tms;

res_mm=1.5;
Kmax=FOVx/(2*res_mm);

Nres=FOVx/res_mm;
kMax_radm=pi*1000/res_mm;


Fac=(kMax_radm/Kmax)*radm2cm;

Traj1cm=Rows2C(BARTTrajx)*Fac;
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
figure;
subplot(2,2,1);
plot(Rows2C(Traj))
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
kTrajO=C2Rows(kC.'/Fac).'; % Nx2
%%
kTraj=kTrajO.';
nTrajG=size(kTraj,2);
nInnerInterleaves=20;
TimePerAcqPoint_us=2.5;
MaxK=max(BARTTrajx(:));
Acc=ceil(MaxK*2).^2/nTrajG;
dkToGFac=FOVx/1000/2/pi;

gammaMHz=42.574; % Hz/T
TwoPiGammaMHz=gammaMHz*2*pi;

Rows2C=@(x) x(1,:) + 1i*x(2,:);
C2Rows=@(x) [real(x); imag(x)];
% kTraj is in rad/m
CkTraj_radTom=Rows2C(kTraj);
CkTraj_k=CkTraj_radTom*FOVx/1000/2/pi; % +- Half nres 

CGrad=[0 diff(CkTraj_radTom)]*1000/TwoPiGammaMHz/10;

CSlew=diff(CGrad)*100;

MaxGrad=max(abs(CGrad));
MaxSlew=max(abs(CSlew));

Ts=(0:(numel(CkTraj_k)-1))*10;
%%
OutInPlot;
CGradO=CGrad;
%%
addpath(genpath('/autofs/space/daisy_002/users/Gilad/mintgrad/'));
%%
mTmToGcm=0.1;
TmsToGcms=100;
radmTo1cm=1/100/2/pi;
Id=30;
Nest=Id;
g0=[0 0];
gf=[real(CGrad(Id)) imag(CGrad(Id))]*mTmToGcm;
% gf=[imag(CGrad(Id)) real(CGrad(Id))]*mTmToGcm;
deltakMoment=[real(CkTraj_radTom(Id)) imag(CkTraj_radTom(Id))]*radmTo1cm; % /cm, s/cm, s^2/cm
T		= 10/1e6;  % Sampling time (s).
gmax=4;
smax=15000;
%	    gmax	= 1x1;  Maximum gradient amplitude (G/cm). 4,15000,0,3);
%	    smax	= 1x1;	Maximum slew rate (G/cm/s).
t0	=0;
type=3;

[g,v] = mintimegrad(Nest,g0,gf,deltakMoment,T,gmax,smax,t0,type)
Cg=Rows2C(g.')/mTmToGcm;
%%
figure;
plot(CGrad(1:Id),'r'); hold on
plot(CGrad(Id:100),'k');
plot(Cg,'g')
%%
csCGrad=cumsum(CGrad);
csCG=cumsum(Cg);
figure;
plot(csCGrad(1:Id),'r'); hold on
plot(csCGrad(Id:100),'k');
plot(csCG,'g')
%%
kTraj=kTrajO(1:300,:).'; % 2xN
OutInPlot;
%%
CurGrad=[Cg CGradO(Id+0:end)];
% CGrad=[0 diff(CkTraj_radTom)]*1000/TwoPiGammaMHz/10;
% kTraj=C2Rows(cumsum(CGrad(1:260))*TwoPiGammaMHz*10/1000);
kTraj=C2Rows(cumsum(CurGrad(1:260))*TwoPiGammaMHz*10/1000);

% kTraj=kTrajO(1:300,:).';
OutInPlot;
