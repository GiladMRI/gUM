M=loadniidata('/media/a/DATA1/FromRoy/S021/niftis/fieldmap/RoyHaa_270715_C001S021_20150727_001_006_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E00_M.nii.gz');
P=loadniidata('/media/a/DATA1/FromRoy/S021/niftis/fieldmap/RoyHaa_270715_C001S021_20150727_001_008_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E00_P.nii.gz');
C=double(M.*exp(1i*(P*2*pi/4096)));

save('Brain128x128x38x32ch.mat','C');
% Focus on sli 10
s=10;
CurSliD=squeeze(C(:,:,s,:));
%%
Senss1=RunESPIRiTForSensMaps(CurSliD,38,[70 70]);
ShowAbsAngle(Senss1)
%%
s2=22;
CurSliD2=squeeze(C(:,:,s2,:));
%%
Senss2=RunESPIRiTForSensMaps(CurSliD2,38,[70 70]);
ShowAbsAngle(Senss2)
%% SVD
DataBoth=CurSliD+CurSliD2;
[sccmtx, Energy]=SCCfBySlices(DataBoth);
ncc=8;
Senss1C=ApplySCCBySlices(Senss1,sccmtx,ncc);
Senss2C=ApplySCCBySlices(Senss2,sccmtx,ncc);
SensBoth=cat(4,Senss1C,Senss2C);
save([BaseRoyP 'SensBoth.mat'],'SensBoth');