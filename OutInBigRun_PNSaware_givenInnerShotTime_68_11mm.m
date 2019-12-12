VD=1;
FOV=240;

FullROTime_ms=48;
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

ResV=1.1;
nRes=numel(ResV);

disp('Base params set');
%
G6ms=cell(1,nRes);
S6ms=cell(1,nRes);
G8ms=cell(1,nRes);
S8ms=cell(1,nRes);

MaxKyShow=5;
%%
InnerShotTime_ms=8;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad8ms=2*pi/(2*nInnerShots);

[G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.67,true,InnerShotTime_ms,42.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(G8ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad8ms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);
%%
InnerShotTime_ms=12;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad12ms=2*pi/(2*nInnerShots);

[G12ms{1}, S12ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,4.3,true,InnerShotTime_ms,42.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(G12ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad12ms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);
%%
InnerShotTime_ms=16;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad16ms=2*pi/(2*nInnerShots);

[G16ms{1}, S16ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,6.0,true,InnerShotTime_ms,37.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(G16ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad16ms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);
%%
InnerShotTime_ms=24;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad24ms=2*pi/(2*nInnerShots);

[G24ms{1}, S24ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,9.3,true,InnerShotTime_ms,42,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(G24ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad24ms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);
%%
Grads11mm=[G8ms G12ms G16ms G24ms];
%%
nInnerShots=[6 4 3 2];
InterShotPhiRads=2*pi./(2*nInnerShots);
Order8ms=[0 4 2 5 3 1];
Order12ms=[0 2 1 3];
Order16ms=[0 1 2];
Order24ms=[0 1];
Orders={Order6ms Order8ms Order12ms Order16ms Order24ms};
clear Grads
for r=1:numel(Grads11mm)
    Grads{r}=Row(Grads2mm{r}.*exp(1i*Orders{r}*InterShotPhiRads(r)));
%     Grads{r+5}=Row(Grads2mm{r+5}.*exp(1i*Orders{r}*InterShotPhiRads(r)));
end
%%
Ttls={'8ms VD1  1.1mm','12ms VD1  1.1mm','16ms VD1  1.1mm','24ms VD1  1.1mm'};
BasePG='/autofs/cluster/kawin/Gilad/Grads/';
BaseIdx=12;
for i=1:numel(Grads)
    WriteGradFile([BasePG 'Traj' num2str(BaseIdx+i)],Grads{i}.',Ttls{i})
end
disp('Saved');
%%
save('Grads11mm.mat','Grads11mm','Grads');
%% That's it



%%
for r=1:nRes
    for i=1:2
        gPlotTraj_radm(GC{r,i},FOV,prisma);MaximizeFig;
        subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(GC{r,i}.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);
    end
end
%% Show advance along time - every 10 ms
for r=1:numel(ResV)
    figure;
    for i=1:2
        GTrajaC=GC{r,i}.';
        g=GTrajaC;
        k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
        kKA=k*FOV_mm/1000/2/pi;
        for k=1:4
            CurIdxs=(k-1)*1000+(1:1000);
            subplot(2,4,k+(i-1)*4);
            plot(kKA(CurIdxs));
            axis([-1 1 -1 1]*FOV/ResV(r)/2);axis square;
            if(k==1)
                ylabel([num2str(ResV(r)) 'mm']);
            end
        end
    end
end
%% Show advance along time - every 20 ms
for r=1:numel(ResV)
    figure;
    for i=1:2
        GTrajaC=GC{r,i}.';
        g=GTrajaC;
        k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
        kKA=k*FOV_mm/1000/2/pi;
        for k=1:4
            CurIdxsS=(k-1)*1000+1;
            CurIdxsE=min(CurIdxsS+1999,numel(kKA));
            CurIdxs=CurIdxsS:CurIdxsE;
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
GAll=cat(3,GCM{:});
%%
save('GAll68.mat','GAll');