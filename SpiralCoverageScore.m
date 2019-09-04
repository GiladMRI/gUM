FOVx=200;
dFOV=FOVx/1000;
paramLongROSamples=1024*6;
AccR=1.4;

paramLongROSamples=1024*3;
AccR=2.8;

spBW=4e5;
paramLongInterleaves=1;
VD=1.3;
paramLongSpGradAmp=35;
paramLongSpSlewRate=170;
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
        paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
    
kTrajM=kTraj;
kTrajM(:,1)=-kTrajM(:,1);
kTraj=[kTraj; kTrajM];
% kTraj=[kTraj; -kTraj];
clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;

EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;
HalfN=ceil(EffMaxRes/2);
% figure;plot(BARTTrajx(1,:),BARTTrajx(2,:),'.');
% grid on;
% set(gca,'XTick',-HalfN:HalfN);
% set(gca,'YTick',-HalfN:HalfN);
%
nTrajX=size(BARTTrajx,2);
C=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
M=zeros(2*HalfN+1);
for i=-HalfN:HalfN
    for j=-HalfN:HalfN
        Dist=C-i-1i*j;
        M(i+HalfN+1,j+HalfN+1)=min(abs(Dist));
    end
end
[X Y]=ndgrid(-HalfN:HalfN,-HalfN:HalfN);
R=sqrt(X.^2+Y.^2);
MB=M;
MB(R>HalfN-5)=0;
fgmontage(MB);
title(['Max dist = ' num2str(max(MB(:)))]); 
xlabel(HalfN);
ylabel(AccR);
%%
GFac=100;
NBig=HalfN*2*GFac;
HalfNBig=ceil(NBig/2);
MBig=zeros(NBig);
for i=1:nTrajX
    MBig(ceil(BARTTrajx(1,i)*GFac+HalfNBig),ceil(BARTTrajx(2,i)*GFac+HalfNBig))=MBig(ceil(BARTTrajx(1,i)*GFac+HalfNBig),ceil(BARTTrajx(2,i)*GFac+HalfNBig))+1;
end
Kernel=fspecial('Gaussian',GFac*[5 5],GFac*1);
MBigK=conv2(MBig,Kernel./max(Kernel(:)),'same');