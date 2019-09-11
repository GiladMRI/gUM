VD=1.3;
FOV=200;

FullROTime_ms=25;
Gmax_mTm=30;
Smax_Tms=170;
GradDwellTime_ms=10e-3;

ma7t=safe_params_ma7T();

FibPhi=(1+sqrt(5))/2;

GoldenAngle=180/FibPhi;
InterShotPhiRad=GoldenAngle*pi/180;

GoldenAnglex=360/(1+FibPhi);
InterShotPhiRad=GoldenAnglex*pi/180;
PNSThresh=95; % percent

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
FOV_mm=FOV;

% ResV=[2 1.25];
ResV=1.9;
nRes=numel(ResV);

disp('Base params set');
%
G5ms=cell(1,nRes);
S5ms=cell(1,nRes);
%%
InnerShotTime_ms=5;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad5ms=2*pi/(2*nInnerShots);

% setting by PNS 95%
% [G5ms{1}, S5ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.8,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,ma7t,PNSThresh,GoldenAngle); 
[G5ms{1}, S5ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.60,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,ma7t,PNSThresh,GoldenAngle); 
% [G5ms{2}, S5ms{2}]=OutInSingleInnerShotf(FOV,ResV(2),VD,2.07,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,ma7t,PNSThresh,GoldenAngle);

tmp=Row(G5ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad5ms));
gPlotTraj_radm(tmp,FOV,ma7t);MaximizeFig;
%
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
% That's it
Order5ms=[0 3 1 4 2];

G5msx=Row(G5ms{1}.*exp(1i*Order5ms*InterShotPhiRad5ms));
GC{1}=[G5msx.'];
%%
% for r=1:nRes
%     for i=1:1
%         gPlotTraj_radm(GC{r,i},FOV,ma7t);MaximizeFig;
%         subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(GC{r,i}.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);
%     end
% end
%% Show advance along time - every 10 ms
for r=1:numel(ResV)
    figure;
    for i=1:1
        GTrajaC=GC{r,i}.';
        g=GTrajaC.';
        k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
        kKA=k*FOV_mm/1000/2/pi;
        for k=1:4
            CurIdxs=(k-1)*500+(1:1000);
            subplot(2,4,k+(i-1)*4);
            plot(kKA(CurIdxs));
            axis([-1 1 -1 1]*FOV/ResV(r)/2);axis square;
            if(k==1)
                ylabel([num2str(ResV(r)) 'mm']);
            end
        end
    end
end
%% GAll: [5000×5×3 double]
for r=1:nRes
    GCM{r}=cat(1,GC{r,:}).';
end
GAll=cat(3,GCM{:}.');
%%
save('GAll5x5.mat','GAll');
disp('Saved');