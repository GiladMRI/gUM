FN='/media/a/DATA/11Jul18/RL/ASL_EPI.nii';
A=gflip(permute(loadniidata(FN),[2 1 3 4]),[1 3]);
AE=A(:,:,:,2:2:end);
AO=A(:,:,:,1:2:end);
D=AE-AO;
Perf=mean(D(:,:,:,3:35),4);
M=mean(A,4);
%%
Fore=Perf;
Back=M;
PStr='EPI';
ClimF=[3 20];
ClimB=[10 2500];

WhichSlicesShow=1:12;

Fore(Back<0.1)=0;
A=(min(max(Fore,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
clear RGB
for s=1:numel(WhichSlicesShow)
    RGB(:,:,:,s)=ind2rgb(round(squeeze(A(:,:,WhichSlicesShow(s)))*255)+1,parula(256));
end
Msk=permute(Fore(:,:,WhichSlicesShow)>ClimF(1),[1 2 4 3]);

B=(min(max(Back,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+permute(B(:,:,WhichSlicesShow),[1 2 4 3]).*(1-Msk);
% Y=permute(X,[1 2 4 3]);
Y=PartitionDim(X,4,3);
Z=CombineDims(Y,[4 2]);
Z=CombineDims(Z,[4 1]);
figure;imshow(Z)
title(PStr)
xlabel(FN);