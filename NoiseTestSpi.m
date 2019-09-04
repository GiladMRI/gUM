FN='meas_MID38_gBP_ep2d_bold_multiecho_Noise_FID33419.dat';

AData = mapVBVD(FN);
N=squeeze(AData.image(1:100,:,:,:,1:12,3,:,:,:,:,:,:,:,:,:));

N=squeeze(AData.image(:,:,:,:,1:12,3,:,:,11,:,:,:,:,:,:));

N=N(:,:,1:12,:,:);
% ADataI=AData.image();
% ADataIx=AData.image(:,:,:,:,:,3,:,:,:,:,:,:,:);
