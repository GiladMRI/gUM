nSlices=16;
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

BaseTFP3='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns/MLN/';
MLNPrefix3='A3S';
[MLNRXS3, SlicesFound3]=CollectMLNResultsfx(BaseTFP3,MLNPrefix3,nSlices);
MLNRXS3=gflip(MLNRXS3,2);

BaseTFP1='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00872_FID32107_ep2d_ge_EPTI_1p9_1shot_4dyns/MLN/';
MLNPrefix1='S';
[MLNRXS1, SlicesFound1]=CollectMLNResultsfx(BaseTFP1,MLNPrefix1,nSlices);

BaseTFPS='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1/MLN/';
MLNPrefixS3='gB03S';
MLNPrefixS1='gB0S';
[SMLNRXS3, SSlicesFound3]=CollectMLNResultsfx(BaseTFPS,MLNPrefixS3,nSlices);
[SMLNRXS1, SSlicesFound1]=CollectMLNResultsfx(BaseTFPS,MLNPrefixS1,nSlices);
SMLNRXS1=SMLNRXS1(:,:,:,ROrd);
SMLNRXS3=SMLNRXS3(:,:,:,ROrd);
SxMLNRXS3=CollectMLNResultsf('MID00860','ME_3S',nSlices);

% ESAllRec=double(cat(5,MLNRXS1,MLNRXS3,padarray(SMLNRXS1,[3 3],'Both'),padarray(SxMLNRXS3,[3 3],'Both')));
ESAllRec=double(cat(5,MLNRXS1,MLNRXS3,padarray(SMLNRXS1,[3 3],'Both'),padarray(SMLNRXS3,[3 3],'Both')));
ESAllRec=ESAllRec./grms(ESAllRec,1:3);

fgmontagex((squeeze(ESAllRec(:,:,5,[4 11],:))))
XSz=[2 4];
Ttls={'EPTI 1shot','EPTI 3shot','SKEPTIC 1shot','SKEPTIC 3shot'};
AddTopTtls(Ttls,XSz)

fgmontagex((squeeze(ESAllRec(:,:,[3 7],4,:))));caxis([0 3]);
AddTopTtls(Ttls,XSz)

fgmontagex((squeeze(ESAllRec(:,:,[3 7],11,:))));caxis([0 3]);
AddTopTtls(Ttls,XSz)
%%
EPTI1ShotP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00872_FID32107_ep2d_ge_EPTI_1p9_1shot_4dyns/';
load([EPTI1ShotP 'resGRAPPASquashed.mat'],'resS','resbS');
%%
directory_rawdata_fully = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
FNBaseFully='meas_MID00868_FID32103_ep2d_ge_sms1_EPTI_1p9_fully';
load([directory_rawdata_fully FNBaseFully filesep 'Recs411.mat'],'idata4','idata11');
%%
SG4=grmss(ifft2cg(resS(:,:,:,4)),3);
SG4b=grmss(ifft2cg(resbS(:,:,:,4)),3);
F4=grmss(idata4(:,1:120,15,:),4);
F4b=grmss(idata4(:,1:120,45,:),4);
%%
fgmontagex(SG4b);colorbar
fgmontagex(F4b);colorbar
fgmontagex(F4);colorbar
fgmontagex(circshift(padarray(MLNRXS3(:,:,7,4)/256,[0 0],'Both'),1,2));colorbar;title('MLNRXS3')
fgmontagex(circshift(padarray(SxMLNRXS3(:,:,7,4)/256,[2 2],'Both'),5,2));colorbar;title('SxMLNRXS3');caxis([0 1.2]);
fgmontagex(circshift(padarray(MLNRXS1(:,:,7,4)/256,[2 2],'Both'),5,2));colorbar;title('MLNRXS1')
fgmontagex(circshift(padarray(SMLNRXS1(:,:,7,4)/256,[2 2],'Both'),5,2));colorbar;title('SMLNRXS1')