% RES::::::::: 128
% 01 TEndOut 0.08245, NCur 32.4148824
% 02 TEndOut 0.00224, NCur 16.7074412
% 03 TEndOut 0.00466, NCur 8.8537206
% 04 TEndOut 0.00907, NCur 4.9268603
% 05 TEndOut 0.01652, NCur 6.8902905
% 06 TEndOut 0.01173, NCur 5.9085754
% 07 TEndOut 0.01373, NCur 6.3994329
% 08 TEndOut 0.01265, NCur 6.1540041
% 09 TEndOut 0.01317, NCur 6.2767185
% 10 TEndOut 0.01290, NCur 6.3380757
% 11 TEndOut 0.01278, NCur 6.3073971
% 12 TEndOut 0.01284, NCur 6.3227364
% 13 TEndOut 0.01281, NCur 6.3304061
% 14 TEndOut 0.01279, NCur 6.3265713
% 15 TEndOut 0.01280, NCur 6.3284887
% VDSpiral
% --------
% FOV 1.920000e-01
% res 1.503977e-03
% AccR 6.328489e+00
% nInterleaves 1.000000e+00
% alpha 1.700000e+00
% MaxGradAmp 5.570423e-03
% MaxSlewRate 2.466902e+01
% VecLen 1280
% Maximum gradient Actual: 3.472300e-02 T/m
% 
% 
% 
% RES::::::::: 130
% 01 TEndOut 0.05840, NCur 33.0469316
% 02 TEndOut 0.00155, NCur 17.0234658
% 03 TEndOut 0.00324, NCur 9.0117329
% 04 TEndOut 0.00631, NCur 5.0058664
% 05 TEndOut 0.01152, NCur 3.0029332
% 06 TEndOut 0.01932, NCur 4.0043998
% 07 TEndOut 0.01444, NCur 4.5051331
% 08 TEndOut 0.01282, NCur 4.7554998
% 09 TEndOut 0.01213, NCur 4.6303165
% 10 TEndOut 0.01246, NCur 4.5677248
% 11 TEndOut 0.01264, NCur 4.5364290
% 12 TEndOut 0.01273, NCur 4.5207811
% 13 TEndOut 0.01277, NCur 4.5129571
% 14 TEndOut 0.01279, NCur 4.5090451
% 15 TEndOut 0.01281, NCur 4.5110011
% 16 TEndOut 0.01280, NCur 4.5100231
% VDSpiral
% --------
% FOV 1.920000e-01
% res 1.474771e-03
% AccR 4.510023e+00
% nInterleaves 1.000000e+00
% alpha 1.700000e+00
% MaxGradAmp 5.570423e-03
% MaxSlewRate 2.466902e+01
% VecLen 1281
% Maximum gradient Actual: 2.954327e-02 T/m
%%
ACCsX=1:0.01:7;
for k=1:numel(ACCsX)
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,ACCsX(k),paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
    EMR(k)=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;
end
figure;plot(ACCsX,EMR)
%%
figure;
ACCsY=[2.56 2.57 4.8 5.5 7];
MM=140;
for y=1:numel(ACCsY)
    CurACC=ACCsY(y);
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,CurACC,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
    CEMR=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;
    SlewBuf=diff(GradBuf,1,1)*MaxGrad;
    kTrajX=kTraj*dFOV/(2*pi);
    gsubplot(numel(ACCsY),3,y,1);
    plot(kTrajX(:,1),kTrajX(:,2));
    setXaxis([-MM MM]/2);
    setYaxis([-MM MM]/2);
    title(CurACC)
    ylabel(CEMR);
    gsubplot(numel(ACCsY),3,y,2);
    plot(GradBuf*MaxGrad*1000);title(MaxGrad*1000);
    setYaxis([-paramLongSpGradAmp paramLongSpGradAmp]);
    setXaxis([0 1300]);
    gsubplot(numel(ACCsY),3,y,3);
    plot(SlewBuf*100000);
    hold on;
    plot([0 1280],[paramLongSpSlewRate paramLongSpSlewRate],'k');
    plot([0 1280],-[paramLongSpSlewRate paramLongSpSlewRate],'k');
    setXaxis([0 1300]);
    setYaxis([-1.3 1.3]*paramLongSpSlewRate);
end