VD=1.3;
FOV=200;

FullROTime_ms=1;
Gmax_mTm=80;
Smax_Tms=170;
GradDwellTime_ms=10e-3;

prisma=safe_params_prisma();

FibPhi=(1+sqrt(5))/2;

GoldenAngle=180/FibPhi;

PNSThresh=95; % percent

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
FOV_mm=FOV;

disp('Base params set');
%%
clear GC GSC
%%
% ResV=[2 1.5 1];
ResV=15;
% nLoopsVec=[1.2]; % Out-In
nLoopsVec=[3.3]; % Out
for r=1:numel(ResV)
    for i=1:numel(nLoopsVec)
        disp([datestr(now) ' ' num2str([r i])]);
%         [GC{r,i} GSC{r,i}]=OutInMultiInnerShotPNSf(FOV,ResV(r),VD,nLoopsVec(i),FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
        [GC{r,i} GSC{r,i}]=OutMultiInnerShotPNSf(FOV,ResV(r),VD,nLoopsVec(i),FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
    end
end
disp('Finished');
%%
r=1;
for i=1:numel(nLoopsVec)
    gPlotTraj_radm(GC{r,i},FOV,prisma);MaximizeFig;% subplot(2,3,2);title(GSC{r,i}.nInner);
    subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(GC{r,i},GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);
end
%%
GNav=GC{r,i};
save('GNav1ms.mat','GNav');