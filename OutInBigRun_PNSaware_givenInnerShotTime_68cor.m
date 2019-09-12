VD=1.3;
FOV=200;

FullROTime_ms=48;
Gmax_mTm=28;
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

ResV=[1.9];
nRes=numel(ResV);

disp('Base params set');
%
G6ms=cell(1,nRes);
S6ms=cell(1,nRes);
G8ms=cell(1,nRes);
S8ms=cell(1,nRes);
%%
Gmax_mTm=30;
InnerShotTime_ms=6;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad6ms=2*pi/(2*nInnerShots);

% setting by PNS 95%
% [G6ms{1}, S6ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.8,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% [G6ms{1}, S6ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,3.0,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
[G6ms{1}, S6ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,3.31,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 

tmp=Row(G6ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad6ms));
gPlotTraj_radm(tmp,FOV,prisma,[1 3 2]);MaximizeFig;
% small test
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
set(findall(gcf,'-property','FontSize'),'FontSize',22)
%%
% gPlotTraj_radm(tmp,FOV,prisma,[1 2 3]);MaximizeFig;subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
% gPlotTraj_radm(tmp,FOV,prisma,[1 3 2]);MaximizeFig;subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
% gPlotTraj_radm(tmp,FOV,prisma,[3 2 1]);MaximizeFig;subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
% 
% skyra=safe_params_skyra();
% trio=safe_params_trio();
% ma7t=safe_params_ma7T();
% 
% gPlotTraj_radm(tmp,FOV,skyra);MaximizeFig;
% gPlotTraj_radm(tmp,FOV,trio);MaximizeFig;
% gPlotTraj_radm(tmp,FOV,ma7t);MaximizeFig;
% 
% [G6msHighPNS, S6msHighPNS]=OutInSingleInnerShotf(FOV,ResV(1),VD,3.5,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% tmpHighPNS=Row(G6msHighPNS.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad6ms));
% gPlotTraj_radm(tmpHighPNS,FOV,prisma);MaximizeFig;subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmpHighPNS.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
%%
Gmax_mTm=80;
Gmax_mTm=30;
InnerShotTime_ms=8;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRad8ms=2*pi/(2*nInnerShots);

% [G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,4.25,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,4.68,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(G8ms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRad8ms));
gPlotTraj_radm(tmp,FOV,prisma,[1 3 2]);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 3.5]);
set(findall(gcf,'-property','FontSize'),'FontSize',22)
%%
clear G6msx G8msx
nInnerShots6=8;
mod((nInnerShots6/FibPhi)*(0:(nInnerShots6-1)),nInnerShots6)
nInnerShots8=6;
mod((nInnerShots8/FibPhi)*(0:(nInnerShots8-1)),nInnerShots8)

Order6ms=[0 5 2 7 4 1 6 3];
Order8ms=[0 4 2 5 3 1];
for r=1:nRes
    G6msx{r}=Row(G6ms{r}.*exp(1i*Order6ms*InterShotPhiRad6ms));
    G8msx{r}=Row(G8ms{r}.*exp(1i*Order8ms*InterShotPhiRad8ms));
end
GC=[G6msx.' G8msx.'];
%%
for r=1:nRes
    for i=1:2
        gPlotTraj_radm(GC{r,i},FOV,prisma,[1 3 2]);MaximizeFig;
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
save('GAll68cor.mat','GAll');
disp('Saved');