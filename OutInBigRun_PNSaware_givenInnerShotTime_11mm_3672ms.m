VD=1;
FOV=240;

FullROTime_ms=36;
Gmax_mTm=80;
Smax_Tms=170;
GradDwellTime_ms=10e-3;

prisma=safe_params_prisma_bad();

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
MaxKyShow=6;
%%
InnerShotTime_ms=7.2;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRadms=2*pi/(2*nInnerShots);

% [G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.67,true,InnerShotTime_ms,42.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[Gms{1}, Sms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.15,true,InnerShotTime_ms,40.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(Gms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRadms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);

Gmsx{1}=Gms{1};
%%
InnerShotTime_ms=9;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRadms=2*pi/(2*nInnerShots);

% [G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.67,true,InnerShotTime_ms,42.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[Gms{1}, Sms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.9,true,InnerShotTime_ms,40.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(Gms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRadms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);

Gmsx{2}=Gms{1};
%%
InnerShotTime_ms=12;
nInnerShots=FullROTime_ms/InnerShotTime_ms;
InterShotPhiRadms=2*pi/(2*nInnerShots);

% [G8ms{1}, S8ms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,2.67,true,InnerShotTime_ms,42.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
[Gms{1}, Sms{1}]=OutInSingleInnerShotf(FOV,ResV(1),VD,4.2,true,InnerShotTime_ms,40.5,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);

tmp=Row(Gms{1}.*exp(1i*(0:(nInnerShots-1))*InterShotPhiRadms));
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp.',GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);setYaxis([0 MaxKyShow]);

Gmsx{3}=Gms{1};
%%
InnerShotTime=[7.2 9 12];
nInnerShots=FullROTime_ms./InnerShotTime;
InterShotPhiRads=2*pi./(2*nInnerShots);
Orders{1}=[0 3 1 4 2];
Orders{2}=[0 2 1 3];
Orders{3}=[0 1 2];
clear Grads
for r=1:numel(Gmsx)
    Grads{r}=Row(Gmsx{r}.*exp(1i*Orders{r}*InterShotPhiRads(r)));
%     Grads{r+5}=Row(Grads2mm{r+5}.*exp(1i*Orders{r}*InterShotPhiRads(r)));
end
%%
Ttls={'7.2ms VD1 1.1mm 36ms','9ms VD1 1.1mm 36ms','12ms VD1 1.1mm 36ms'};
BasePG='/autofs/cluster/kawin/Gilad/Grads/';
BaseIdx=17;
for i=1:numel(Grads)
    WriteGradFile([BasePG 'Traj' num2str(BaseIdx+i)],Grads{i}.',Ttls{i})
end
disp('Saved');
%%
save('Grads11mm36ms.mat','Gmsx','Grads');
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