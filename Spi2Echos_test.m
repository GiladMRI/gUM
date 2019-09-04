FOVx=200;
dFOV=FOVx/1000;
paramLongROSamples=1024*4;
AccR=19.8;
nInnerInterleaves=6;

% paramLongROSamples=1024*15;
% AccR=2.8;

GRAD_RASTER_TIME=10;

spBW=4e5;

TrajPointsPerGrad=GRAD_RASTER_TIME*spBW/1e6;
paramLongInterleaves=1;
VD=1.3;
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;
% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);

% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,nInnerInterleaves,40]);

dkToGFac=FOVx/1000/2/pi;
BARTTrajx=kTraj.'*FOVx/1000/2/pi;
nTraj=size(BARTTrajx,2);
C=BARTTrajx(1,1:size(BARTTrajx,2)/2)+1i*BARTTrajx(2,1:size(BARTTrajx,2)/2);
%
% figure;plot(BARTTrajx(1,1:end),BARTTrajx(2,1:end),'.')
% figure;plot(BARTTrajx(1,1:end-40),BARTTrajx(2,1:end-40),'b.')
% hold on
% plot(BARTTrajx(1,end-39:end-19),BARTTrajx(2,end-39:end-19),'k.')
% plot(BARTTrajx(1,end-18:end),BARTTrajx(2,end-18:end),'r.')
%
MaxK=max(BARTTrajx(:));
Acc=ceil(MaxK*2).^2/nTraj;

figure;
subplot(2,2,1);
plot3(1:size(BARTTrajx,2),BARTTrajx(1,:),BARTTrajx(2,:))
%
HnTraj=nTraj/nInnerInterleaves;
CLRs='rgbcmkrg';
subplot(2,2,2);plot(BARTTrajx(1,1:HnTraj),BARTTrajx(2,1:HnTraj),[CLRs(1) '.'])
hold on
for i=2:nInnerInterleaves
plot(BARTTrajx(1,HnTraj*(i-1)+(1:HnTraj)),BARTTrajx(2,HnTraj*(i-1)+(1:HnTraj)),[CLRs(i) '.'])
end
axis equal
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
%
subplot(2,2,3);
plot(GradBuf*MaxGrad*1000);title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
setXaxis([1 nTraj]);

SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(SlewBuf*100);
MaxSlew=max(max(abs(SlewBuf(20:HnTraj,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
setXaxis([1 nTraj]);
%%
figure;
plot(real(C),imag(C),'b.');
hold on;
D=-C;
PhiEnd=angle(C(end));
dPhiEnd=angle(C(end))-angle(C(end-1));
n=ceil(pi/dPhiEnd);
dPhiBetter=(pi)/n;
Ext=abs(C(end)).*exp(1i*(PhiEnd+(1:n)*dPhiBetter));
CE=[C Ext];
plot(real(Ext),imag(Ext),'g.');
PhiEnd=angle(CE(end));
Ang=angle(CE);
dAng=diff(Ang);
dAngR=dAng(end:-1:1);
cdAngR=cumsum(dAngR);
% FC=CE(end:-1:1);
% D=abs(FCE(1:end-1)).*exp(1i*(cdAngR+pi+PhiEnd));
% D=-FC;
FC=C(end:-1:1);
% FCE=CE(end:-1:1);
D=-FC;
% D=abs(FC).*exp(1i*(-angle(FC)+PhiEnd*2*0));
% D=abs(FC).*exp(1i*(angle(FC)+pi));
% D=abs(FCE).*exp(1i*(angle(FCE)+pi));
% D=abs(FC).*exp(1i*(angle(FC)+pi*2/3));
% D=-C;
plot(real(D),imag(D),'r.');
% D=abs(FC).*exp(1i*(angle(FC)+pi*4/3));
% plot(real(D),imag(D),'g');
% D=abs(FC).*exp(1i*(angle(FC)+pi));
% plot(real(D),imag(D),'r');
M=max(abs(C))+5;
axis equal
axis([-M M -M M]);
Both=[CE D];

% figure;plot3(1:size(Both,2),real(Both),imag(Both))
%%
dGammaRad			= 267512897.638;
FacA=((GRAD_RASTER_TIME*1e-6)*dGammaRad)/1000;
CBase=C/dkToGFac;
% dTraj=diff(CBase)*1000/((10*1e-6)*267512897.638); % mT/m
dTraj=diff(CBase)/FacA; % mT/m
% dTraj=dTraj*1000/(FOVx/1000/2/pi);
% figure;plot(abs(dTraj))
dREnd=abs(dTraj(end)); % mT/m
usNeeded=dREnd*1000/paramLongSpSlewRate;
GradAccPerGradStep=10*paramLongSpSlewRate/1000;
TrajPointsNeededToStop=ceil(usNeeded/GRAD_RASTER_TIME);
dC=diff(CBase);
dCEnd=abs(dC(end));
FirstLoc=max(abs(CBase));
dY=linspace(abs(dCEnd),0,TrajPointsNeededToStop);
dYP=dY*1000/FacA; % mT/m
Y=cumsum(dY);
YEnd=Y(end);
YEndA=YEnd/GradAccPerGradStep/FacA;
FirstLocA=FirstLoc/GradAccPerGradStep/FacA;
nNeededToRecoverY=ceil(sqrt(YEndA*2));
nNeededForX=ceil(sqrt(FirstLocA*2));
NeededForRecovery=max(nNeededForX-TrajPointsNeededToStop,nNeededToRecoverY);
YRec=YEnd-(0:NeededForRecovery).^2*(YEnd/(NeededForRecovery^2));
nTotal=TrajPointsNeededToStop+NeededForRecovery;
X=FirstLoc-(0:nTotal).^2*(FirstLoc/(nTotal^2));
X=linspace(FirstLoc,0,nTotal+1);
E=X+1i*([Y YRec]);
% E=E/FacA;
E=E.*exp(1i*angle(CBase(end)));
% E=-E;
%
figure;plot(real(CBase),imag(CBase),'bo')

hold on;
plot(0,0,'r+')
plot(real(E),imag(E),'g.')
%%
dE=diff(E); % mT/m
figure;plot(real(dE),imag(dE),'-')

