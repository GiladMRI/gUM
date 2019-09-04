BaseP='/autofs/space/daisy_002/users/Gilad/';
DatFN='meas_MID90_GRE_2D_G3_ASPIRE_FID41991.dat';
FN=[BaseP DatFN];
AData = mapVBVD(FN);
% ADataI=AData.image();
ADataIx=AData.image(:,:,:,:,:,:,:,:,:,:,:,:,:);
refscan=AData.refscan();
refscan=squeeze(refscan);
refscanP=permute(refscan,[1 3 4 2 5]);
Frefscan=ifft2cg(refscanP);
nSlices=32;
Ord=[1:2:nSlices 2:2:nSlices];
[~, ROrd]=sort(Ord);
Frefscan=Frefscan(:,:,ROrd,:,:);
ADataIx=squeeze(ADataIx);
% 32 slices
% 3rd dim is 3 undersampled
ADataIxP=permute(ADataIx,[1 3 4 2 5]);
F=fft2cg(ADataIxP);
FS=grmss(F,4:5);
grmss(ADataIx,1:4)
grmss(ADataIx,[1 2 4 5])
%% Crop both on RO
%%
SliI=15;
Sens=RunESPIRiTForSensMaps(squeeze(Frefscan(:,:,SliI,:,1)),29);
SzTF=[128 128];
SensX=imresizeBySlices(Sens(:,:,1:16),SzTF);
CurSliData=squeeze(ADataIxP(:,:,SliI,:,:));
CurSliDataX=padarray(CurSliData(43:end-42,:,1:16,1:6),[0 11 0 0],'both');
CurSliDataY=permute(CurSliDataX,[2 1 3 4]);
CurSliDataY(3:6:end,:,:,1:2:end)=0;
CurSliDataY(6:6:end,:,:,2:2:end)=0;
IF=ifft2cg(CurSliDataY); %  # batch_size,H,W,nTSC,MB
IFCombined=sum(IF.*conj(SensX),3);
In=permute(IFCombined,[5 1 2 4 3])*100;
save('/autofs/cluster/kawin/Gilad/TF/srezN/in/AA.mat','SensX','In');
% Msk3[3::6,:,::2,:,:,:]=1
% Msk3[::6,:,1::2,:,:,:]=1