info.GradFile = 'sstraj_45S1B1_15_0_2R1' ;
[kx,ky,kz,adr,h,gxr,gyr,gzr]=readgrads(info,1);

LData=mapVBVD('meas_MID48_sp_ep2d_diff_multiband_91_150_FID817.dat');
LDataI=LData.image();
SLDataI=squeeze(LDataI);
%%
figure
for i=1:8
    subplot(2,4,i);
    plot(grmss(SLDataI(:,:,i),2));
end