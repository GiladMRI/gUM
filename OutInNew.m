% mex VDSpiralMex.cpp VDSpiral.cpp -DMATLAB

FOVx=200;
dFOV=FOVx/1000;
paramLongROSamples=1024*12;
AccR=1.4;
AccR=4;
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

ExtraDelay=5;
% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);

% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/nInnerInterleaves,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,6,1,ExtraDelay]);

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
CLRs='krbbcmkrg';
subplot(2,2,2);plot(BARTTrajx(1,1:HnTraj),BARTTrajx(2,1:HnTraj),[CLRs(1) '-'])
hold on
for i=2:nInnerInterleaves
    plot(BARTTrajx(1,HnTraj*(i-1)+(1:HnTraj)),BARTTrajx(2,HnTraj*(i-1)+(1:HnTraj)),[CLRs(i) '-'])
end
axis equal
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
%
subplot(2,2,3);
Ts=(0:(size(GradBuf,1)-1))*10;
plot(Ts,GradBuf*MaxGrad*1000,'.-');title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
% setXaxis([1 nTraj]);
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel('time (ms)');
ylabel('mT/m');
%
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(Ts(1:end-1),SlewBuf*100);
MaxSlew=max(max(abs(SlewBuf((20+ExtraDelay):(HnTraj-1),:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
% setXaxis([1 nTraj]);
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel('time (ms)');
ylabel('mT/m/s');
%%
CAIPISep_mm=36;
CAIPIDelay_us=200;
CAIPIPeriod_us=200;
OutInOut=true;
CAIPIVec=CAIPIBlips([paramLongSpGradAmp, paramLongSpSlewRate,CAIPISep_mm,CAIPIDelay_us,CAIPIPeriod_us,...
    2560*paramLongROSamples/1024]);

gamma=42.5774806;
cCAIPIVec=cumsum(CAIPIVec)*gamma*10*2*pi/1e6;

% if(OutInOut)
%     BaseLen=size(kTraj,1)/3;
%     TrajIdxs=(1:BaseLen)+BaseLen*(WhichE-1);
%     kTraj=kTraj(TrajIdxs,:);
%     GradBuf=GradBuf(TrajIdxs,:);
%     
%     CAIPIVec=CAIPIVec(TrajIdxs,:);
%     cCAIPIVec=cCAIPIVec(TrajIdxs,:);
% end

EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQRadius
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

% nTrajOrig=size(ADataIsPy,1);
nTraj=size(BARTTrajMS,2);
TimeInMs2=(0:nTraj-1)*2.5/1e3;
T2SCompStr='';

BARTTrajAct=BARTTrajMS(:,1:(nTraj));
% BARTTrajCorr=BARTTrajAct./[1.01; 1.01; 1.0];
BARTTrajCorr=RepDotMult(BARTTrajAct,1./[1.01; 1.01; 1.0]);

disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
% nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
XResmm=dFOV*1000/(2*MaxK);
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));

setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc) ]);
xlabel(['Res ' num2str(XResmm,'%.2f') 'mm']);
ylabel(['Nominal acc: ' num2str(AccR) ' VD=' num2str(VD)]);

subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['                Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
hold on;
plot(CAIPIVec,'k');
subplot(2,2,4);
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
subplot(2,2,3);
plot(cumsum(CAIPIVec));
title('CAIPI cumsum');