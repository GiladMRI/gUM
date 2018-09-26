dFOV=0.192;
paramLongROSamples=5120;
spBW=400000;
AccR=1.9;
paramLongInterleaves=1;
VD=1.1;
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
figure;plot3(kTraj(:,1),kTraj(:,2),1:size(kTraj,1))
figure;plot3(1:size(kTraj,1),kTraj(:,1),kTraj(:,2))
%%
FOVx=dFOV*1000;
% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;

BARTTrajMS=BARTTrajx;
BARTTrajxC=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
for i=2:paramLongInterleaves
    CurRotTrajC=BARTTrajxC*exp(1i*2*pi*(i-1)/paramLongInterleaves);
    BARTTrajMS(1,:,i)=real(CurRotTrajC);
    BARTTrajMS(2,:,i)=imag(CurRotTrajC);
end
BARTTrajMS(3,end)=0;
disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])