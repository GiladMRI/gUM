VD=1.3;
FOV=200;

FullROTime_ms=50;
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

disp('Base params set');
%%
G5ms=cell(1,3);
S5ms=cell(1,3);
G10ms=cell(1,3);
S10ms=cell(1,3);
G5msx=cell(1,3);
G10msx=cell(1,3);
%%
InnerShotTime_ms=5;
InterShotPhiRad5ms=2*pi/20;

[G5ms{1}, S5ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.8,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[G5ms{2}, S5ms{2}]=OutInSingleInnerShotf(FOV,ResV(2),VD,2.2,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[G5ms{3}, S5ms{3}]=OutInSingleInnerShotf(FOV,ResV(3),VD,1.54,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

G5msx=Row(G5ms{3}.*exp(1i*(0:9)*InterShotPhiRad5ms));
gPlotTraj_radm(G5msx,FOV,prisma);MaximizeFig;
%%
InnerShotTime_ms=10;
InterShotPhiRad10ms=2*pi/10;

[G10ms{1}, S10ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,6.28,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[G10ms{2}, S10ms{2}]=OutInSingleInnerShotf(FOV,ResV(2),VD,5.1,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[G10ms{3}, S10ms{3}]=OutInSingleInnerShotf(FOV,ResV(3),VD,3.8,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

G10msx=Row(G10ms{3}.*exp(1i*(0:4)*InterShotPhiRad10ms));
gPlotTraj_radm(G10msx,FOV,prisma);MaximizeFig;
%%
Order5ms=[0 5 3 7 1 6 4 8 2 9];
Order10ms=[0 3 1 4 2];
for r=1:numel(ResV)
    G5msx{r}=Row(G5ms{r}.*exp(1i*Order5ms*InterShotPhiRad5ms));
    G10msx{r}=Row(G10ms{r}.*exp(1i*Order10ms*InterShotPhiRad10ms));
end
GC=[G5msx.' G10msx.'];
%%
for r=1:numel(ResV)
    for i=1:2
        gPlotTraj_radm(GC{r,i},FOV,prisma);MaximizeFig;
    end
end
%%
for m=1:numel(ResV)
    figure;
    for i=1:2
        GTrajaC=GC{m,i}.';
        g=GTrajaC;
        k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
        kK=k*FOV_mm/1000/2/pi;
        subplot(2,2,i*2-1);
        plot(kK)
        zIdxs=find((real(kK(1:end-1)).*real(kK(2:end)))<0);
        hold on;
        plot(kK(zIdxs),'ro');
        subplot(2,2,i*2);
        IkK=imag(kK(zIdxs));
        [SI, IkKOrd]=sort(IkK);
        dSI=diff(SI);
        plot(SI(1:end-1)+dSI/2,dSI)
        
        Crossings=find(  (real(kK(1:end-3)).*real(kK(4:end))<0) & (imag(kK(3:end-1)).*imag(kK(4:end))<0) );
        title(numel(Crossings));
    end
    MaximizeFig;
end
%% Show advance along time - every 10 ms
for r=1:numel(ResV)
    figure;
    for i=1:2
        GTrajaC=GC{r,i}.';
        for k=1:4
            CurIdxs=(k-1)*1000+(1:1000);
            subplot(2,4,k+(i-1)*4);
            g=GTrajaC(CurIdxs);
            k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
            kK=k*FOV_mm/1000/2/pi;
            plot(kK);
            axis([-1 1 -1 1]*FOV/ResV(r)/2);axis square;
        end
    end
end
%% Show advance along time - every 20 ms
for r=1:numel(ResV)
    figure;
    for i=1:2
        GTrajaC=GC{r,i}.';
        for k=1:4
            CurIdxs=(k-1)*1000+(1:2000);
            subplot(2,4,k+(i-1)*4);
            g=GTrajaC(CurIdxs);
            k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
            kK=k*FOV_mm/1000/2/pi;
            plot(kK);
            axis([-1 1 -1 1]*FOV/ResV(r)/2);axis square;
        end
    end
end
%% GAll: [5000×5×3 double]
for r=1:numel(ResV)
    GCM{r}=cat(1,GC{r,:}).';
end
GAll=cat(3,GCM{:});
%%
