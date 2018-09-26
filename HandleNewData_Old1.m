FN='/home/a/SpiralData/27Apr18/meas_MID80_gSpiral_BOLD_1Sli_1Rep_1Shot_VD10_6ADCs_25_125_FID15656.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID78_gSpiral_BOLD_1Sli_1Rep_1Shot_VD13_6ADCs_25_125_FID15654.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID85_gSpiral_BOLD_1Sli_1Rep_2Shot_VD10_4ADCs_25_125_FID15660.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID87_gSpiral_BOLD_1Sli_4Rep_2Shot_VD10_4ADCs_25_125_FID15662.dat';

FN='/media/a/DATA1/13May18/Phantom/meas_MID364_gBP_VD11_U19_FID17753.dat';
%% Read raw
AData = mapVBVD(FN);
ADataI=AData.image();
ADataIs=squeeze(ADataI);

FOVx=AData.hdr.Meas.ReadFOV;
dFOV=FOVx/1000;


paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{4};
spBW =AData.hdr.MeasYaps.sWiPMemBlock.alFree{9};
AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.alFree{5};
VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};
paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.alFree{7};
paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.alFree{6};

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1/4:(size(kTraj,1)-0.1));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1/4:(size(kTraj,1)-0.1));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;
ADataIsP=double(permute(ADataIs,[1 3 2]));
ADataIsPx=reshape(ADataIsP,[paramLongROSamples,3]);
ADataIsPy=ADataIsPx(1:size(BARTTrajx,2),:);
%%
clear recoAD
Sz3=[88 88];
recoAD(:,:,c) = bart('pics -r:0.1 -R T:7:0:0.0010 -t ',BARTTrajx, ADataIsPy(:,c).', ones(Sz3));
ShowAbsAngle(recoAD(:,:,c))
%%
clear recoAD

Sz3=[88 88];
recoAD(:,:,c) = bart('pics -r:0.1 -R T:7:0:0.0010 -t ',BARTTrajx(:,1:end-2), ADataIsPy(3:end,c).', ones(Sz3));
ShowAbsAngle(recoAD(:,:,c))