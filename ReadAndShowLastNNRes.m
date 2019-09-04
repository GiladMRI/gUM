BB=load('/autofs/cluster/kawin/Gilad/TF/srezN/out/OnRealData01.mat');
C=squeeze(BB.x(:,:,:,1)+1i*BB.x(:,:,:,2));

CLast=C(end-127:end,:,:,:);
CLastx=PartitionDim(CLast,2,6);

ShowAbsAngle(CLastx);
subplot(1,2,1);colorbar