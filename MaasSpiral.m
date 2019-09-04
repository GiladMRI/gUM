% FOVx=AData.hdr.Meas.ReadFOV;
% dFOV=FOVx/1000;
% 
% paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{20};
% spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
% AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{6};
% paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
% VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{5};
% paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
% paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
% MB=AData.hdr.MeasYaps.sWiPMemBlock.alFree{9};
% 
% CAIPISep_mm=AData.hdr.MeasYaps.sWiPMemBlock.adFree{7};
% CAIPIPeriod_us=AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
% CAIPIDelay_us=AData.hdr.MeasYaps.sWiPMemBlock.adFree{9};
% % MB
% nSlicesNoMB=nSlices/MB;
% 
% % if(MB>1)
% if(isempty(spBW))
%     spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{14};
%     paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};
%     paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};
%     paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
% end
% disp('Read data');
% %%
% SummaryStr.MB=MB;
% SummaryStr.paramLongSpGradAmp=paramLongSpGradAmp;
% SummaryStr.paramLongSpSlewRate=paramLongSpSlewRate;
% SummaryStr.paramLongROSamples =paramLongROSamples;
% SummaryStr.spBW=spBW;
% SummaryStr.AccR=AccR;
% SummaryStr.paramLongInterleaves=paramLongInterleaves;
% SummaryStr.VD=VD;
% 
% if(MB>1)
%     SummaryStr.CAIPISep_mm=CAIPISep_mm;
%     SummaryStr.CAIPIPeriod_us=CAIPIPeriod_us;
%     SummaryStr.CAIPIDelay_us=CAIPIDelay_us;
% end
% %%
% save([ScanP BaseFN filesep 'Data.mat'],'ADataIsL','AData');
% disp(['Saved ' [ScanP BaseFN filesep 'Data.mat']]);
% %%
% % save('StatusForTrajDesign1Band_2.mat');
% %% given everything but VD and Acc
% % OutInOut=false;
% OutInOut=paramLongROSamples>9000;
% if(OutInOut)
%     [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/3,spBW,AccR,...
%     paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);
% else
%     [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
% end
%%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;
%%
FOV_mm=200;
dFOV=FOV_mm/1000;
paramLongROSamples=1024*5;
spBW=4e5;
AccR=1.6;
paramLongInterleaves=1;
VD=1.3;
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
        paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
%%
kTrajC=kTraj(:,1)+1i*kTraj(:,2);
figure;plot(kTrajC);
%%
GradBufC=(GradBuf(:,1)+1i*GradBuf(:,2))*1000*MaxGrad;

gPlotTraj_radm(GradBufC,FOV_mm);
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(GradBufC,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 7]);

%%

figure;
QuickAnalyzeTrajJumpsf(GradBufC,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm)
%%
GAll68=load('GAll68.mat');
%%
g6=GAll68.GAll(:,1,1);
g6=g6(1:2400);
gPlotTraj_radm(g6,FOV_mm);
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(g6,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 7]);
%%
VD=1.3;
FOV=200;

FullROTime_ms=12;
Gmax_mTm=80;
Smax_Tms=170;
GradDwellTime_ms=10e-3;

prisma=safe_params_prisma();

FibPhi=(1+sqrt(5))/2;

GoldenAngle=180/FibPhi;
InterShotPhiRad=GoldenAngle*pi/180;

GoldenAnglex=360/(1+FibPhi);
InterShotPhiRad=GoldenAnglex*pi/180;
PNSThresh=95; % percent

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
FOV_mm=FOV;

ResV=[1.9 1.3];
nRes=numel(ResV);

disp('Base params set');
%
Gxms=cell(1,nRes);
Sxms=cell(1,nRes);
%%
InnerShotTime_ms=6;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRadxms=2*pi/(2*nInnerShots);

% UseMTG_forLowK=true;
UseMTG_forLowK=false;
% setting by PNS 95%
[Gxms{1}, Sxms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.8,UseMTG_forLowK,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% [Gxms{1}, Sxms{1}]=OutInSingleInnerShotf(FOV,ResV(1),1,2.8,UseMTG_forLowK,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% [Gxms{1}, Sxms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,3.2,UseMTG_forLowK,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% For 12 out in:
% [Gxms{1}, Sxms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,6.3,UseMTG_forLowK,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,0); 
% For 12 only out:
% [Gxms{1}, Sxms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,14.3,UseMTG_forLowK,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,0); 

% Only out: [GC{r,i} GSC{r,i}]=OutMultiInnerShotPNSf(FOV,ResV(r),VD,nLoopsVec(i),FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(Gxms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRadxms)).';
tmp=tmp(1:1200);
gPlotTraj_radm(tmp,FOV,[]);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 7]);
%%
tmp=tmp(1:1200);
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm)

figure;
QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm)