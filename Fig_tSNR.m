BaseXP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/';
MOEt=rot90(loadniidata([BaseXP 'MLN_scaled_odd-even_tSNR.nii']));
Mt=rot90(loadniidata([BaseXP 'MLN_scaled_tSNR.nii']));
EOEt=rot90(loadniidata([BaseXP 'resL1ESPIRiTCCS1A_scaled_odd-even_tSNR.nii']));
Et=rot90(loadniidata([BaseXP 'resL1ESPIRiTCCS1A_scaled_tSNR.nii']));
figure;
ha = tight_subplot(1, 2, 0.01, 0.01, 0.01);
axes(ha(1));
gmontage(Et(7:120,20:108,:),[0 200]);title('ESPIRIT tSNR');
axes(ha(2));
gmontage(Mt(7:120,20:108,:),[0 200]);title('MLN tSNR');

figure;
ha = tight_subplot(1, 2, 0.01, 0.01, 0.01);
axes(ha(1));
gmontage(EOEt(7:120,20:108,:),[0 2]);title('ESPIRIT Perfusion tSNR');
axes(ha(2));
gmontage(MOEt(7:120,20:108,:),[0 2]);title('MLN Perfusion tSNR');