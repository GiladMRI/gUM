load('PSFs_analysis.mat','ToPlotT','AdjsT');
load('PSFs_analysisb.mat','ImT');

load('N2mktC.mat','N2mktC','N2mktCX');
%%
fgmontagex(perm43(squeeze(AdjsTM(:,:,1:12,12,[1 2 4 11 12]))),[0 30])
fgmontagex(perm43(squeeze(AdjsTM(:,:,1:12,12,[1 4 12]))),[0 30])

fgmontagex(perm43(squeeze(ImTM2(:,:,[1 2 4],[1 4 12]))),[0 1e-3])
%%
[az,el] = view
cmp=campos;
Rng=40;
figure;
surface(abs(ImTM2(61+(-Rng:Rng),61+(-Rng:Rng),2,4)));
view(15.60,31.433);removeTicks;set(gca,'ZTick',[]);
campos([87.1758 -215.0226    0.0710]);
% axis([-0.3550   24.6450   -2.5885   18.4115   -0.0019    0.0081]);
%%
figure;plot(abs(CTraj{12}),'k','LineWidth',2)
axis([-500 5300 -20 80]);
%%
fgmontagex(gflip(N2mktCM(:,:,12),1));colormap jet
%%
fgmontagex(perm43(squeeze(ImTM2(:,:,[1 2 4],[1 4 12]))),[0 1e-3])
%%
ImTM3=PartitionDim(ImTM2,3,2);
ImTM3=CombineDims(CombineDims(ImTM3,[5 1]),[3 2]);

ImTM4=log10(abs(ImTM3));
ImTM4(ImTM4<-4)=NaN;

figure;surf(ImTM4(:,:,1));view([az,el]);campos(cmp);removeTicks;set(gcf,'Color','w');grid off;%set(gca,'Visible','off')
figure;surf(ImTM4(:,:,4));view([az,el]);campos(cmp);removeTicks;set(gcf,'Color','w');grid off
figure;surf(ImTM4(:,:,12));view([az,el]);campos(cmp);removeTicks;set(gcf,'Color','w');grid off
%%
N2mktCM=cat(3,N2mktC{:});
ImTM=cat(5,ImT{:});
ImTM2=CombineDims(ImTM,[4 3]);

AdjsTM=cat(5,AdjsT{:});
%%
G2mm=load('Grads2mmb.mat');
% QQ=load('GAll1p9mmVD1PAT3Pause.mat');
QQ=load('GAll1p9mmVD1PAT3Pauseb.mat');
AllGrads=cat(1,G2mm.Grads{:}).';
AllGrads=cat(2,AllGrads,QQ.GAll(:,7:8));
G11mm=load('Grads11mmb.mat');
AllGrads11=cat(1,G11mm.Grads{:}).';
AllGrads=cat(2,AllGrads,AllGrads11);
%%
for TrajIdx=1:12
GTraj=AllGrads(:,TrajIdx);
% nEchosData=size(idatax,3);
% TotalAcqTime_ms=ES_ms*nEchosData;
k=cumsum([0; GTraj])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kK=k*FOV_mm/1000/2/pi;
CTraj{TrajIdx}=kK;
end
%%
TrajIdx=12;
figure;
subplot(4,1,1);
plot(abs(CTraj{2}))
subplot(4,1,2);

%%
TrajIdx=1;
for TrajIdx=1:5
GTraj=AllGrads(:,TrajIdx);
nEchosData=size(idatax,3);
TotalAcqTime_ms=ES_ms*nEchosData;
k=cumsum([0; GTraj])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kK=k*FOV_mm/1000/2/pi;

figure;
subplot(1,2,1);
plot(abs(kK));setXaxis([0 4800]);setYaxis([0 60]);

subplot(2,2,2);
gmontage(AdjsT{TrajIdx}(:,:,12,1:12),[0 60],'Size',[1 12]);

subplot(2,2,3);
gmontage(ImT{TrajIdx}(:,:,1:3),[0 60],'Size',[1 3]);

subplot(2,2,4);
imagesc(flip(N2mktC{TrajIdx},1),[0 1]);colormap jet
end
%%
figure;
subplot(3,4,1);
plot(Traj(1,:),Traj(2,:));
title(CurTtl);
xlabel(['#Traj=' num2str(nTraj)]);
setXaxis([-1 1]*N/2);
setYaxis([-1 1]*N/2);
removeTicks;
subplot(3,4,2:4);
plot((angle(CTrajx)+pi)/(2*pi),'r');hold on;plot(abs(CTrajx)/96,'k','LineWidth',2)
setXaxis([1 nTraj]);
subplot(3,4,[5 6 9 10]);
% gmontage(abs(Im),[0 (max(MxVecXX(2:end))+max(MxVecXX(1:end)))/2]);
gmontage(abs(Im),[0 max(MxVecXX(2:end))*1.2]);
% subplot(3,4,[2 3 4 6 7 8 10 11 12]);
subplot(3,4,[7 8 11 12]);
plot(MxVecXX,'b-','LineWidth',2);hold on
plot(MxVecYY,'k-','LineWidth',2);
plot(MxVecXY,'r-','LineWidth',2);
plot(MxVecYX,'g-','LineWidth',2);
setXaxis([0 N/2]);
legend({'1-1','2-2','1-2','2-1'});

PlotRange=N/2;
CTrajx=Traj(1,:)+1i*Traj(2,:);
figure;
subplot(2,4,1);
plot(Traj(1,:),Traj(2,:));
xlabel(['#Traj=' num2str(nTraj)]);
setXaxis([-PlotRange PlotRange]);
setYaxis([-PlotRange PlotRange]);
removeTicks;
subplot(2,4,2:4);
plot((angle(CTrajx)+pi)/(2*pi),'r');hold on;plot(abs(CTrajx)/PlotRange,'k','LineWidth',2)
setXaxis([1 nTraj]);
removeTicks;
title(CurTtl);
xlabel(num2str(ShotPerm))
subplot(2,4,5:6);
plot(ToPlot(:))
subplot(2,4,7:8);
gmontage(Adjs(:,:,(nB0vals_Hz+1)/2,1:(nB0vals_Hz+1)/2),[0 30])
title('-110:10:0 Hz');

% close all;
% 
% ToPlotT{TrajIdx}=ToPlot;
% AdjsT{TrajIdx}=Adjs;
% end
