VD=1.3;
FOV=200;

FullROTime_ms=50;
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
ResV=[2 1.5 1];
nLoopsVec=[1 1.4 2 3 4 6 8];
for r=1:numel(ResV)
    for i=1:numel(nLoopsVec)
        disp([datestr(now) ' ' num2str([r i])]);
        [GC{r,i} GSC{r,i}]=OutInMultiInnerShotPNSf(FOV,ResV(r),VD,nLoopsVec(i),FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle);
    end
end
%%
r=2;
for i=1:numel(nLoopsVec)
    gPlotTraj_radm(GC{r,i},FOV,prisma);MaximizeFig;subplot(2,3,1);title(GSC{r,i}.nInner);
end
%%
for m=1:numel(ResV)
    figure;
    for i=1:numel(nLoopsVec)
        GTrajaC=GC{m,i};
        g=GTrajaC;
        k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
        kK=k*FOV_mm/1000/2/pi;
        subplot(numel(nLoopsVec),2,i*2-1);
        plot(kK)
        zIdxs=find((real(kK(1:end-1)).*real(kK(2:end)))<0);
        hold on;
        plot(kK(zIdxs),'ro');
        subplot(numel(nLoopsVec),2,i*2);
        SI=sort(imag(kK(zIdxs)));
        dSI=diff(SI);
        plot(dSI)
        
        Crossings=find(  (real(kK(1:end-3)).*real(kK(4:end))<0) & (imag(kK(3:end-1)).*imag(kK(4:end))<0) );
        title(numel(Crossings));
    end
    MaximizeFig;
end
%%
G2mm6L=OutInMultiInnerShotf(FOV,res_mm,VD,6,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm4Le=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,30,70,GradDwellTime_ms);
G2mm4Ld=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,30,100,GradDwellTime_ms);
G2mm4Lc=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,30,120,GradDwellTime_ms);
G2mm4Lb=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm4L=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm3L=OutInMultiInnerShotf(FOV,res_mm,VD,3,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm2L=OutInMultiInnerShotf(FOV,res_mm,VD,2,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm1_4L=OutInMultiInnerShotf(FOV,res_mm,VD,1.4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G2mm1L=OutInMultiInnerShotf(FOV,res_mm,VD,1,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);

gPlotTraj_radm(G2mm8L,FOV)
gPlotTraj_radm(G2mm6L,FOV)
gPlotTraj_radm(G2mm4L,FOV)

gPlotTraj_radm(G2mm4L,FOV)
gPlotTraj_radm(G2mm4Le,FOV)
gPlotTraj_radm(G2mm4Ld,FOV)
gPlotTraj_radm(G2mm4Lc,FOV)
gPlotTraj_radm(G2mm4Lb,FOV)
gPlotTraj_radm(G2mm4L,FOV)

gPlotTraj_radm(G2mm4L,FOV);MaximizeFig
gPlotTraj_radm(G2mm3L,FOV);MaximizeFig
gPlotTraj_radm(G2mm2L,FOV);MaximizeFig
gPlotTraj_radm(G2mm1_4L,FOV);MaximizeFig
gPlotTraj_radm(G2mm1L,FOV);MaximizeFig
%%
g=G2mm3L;
g=G15mm2L;
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
kK=k*FOV/1000/2/pi;
AcqPointsPerGrad=10;
kPoints=interp1(1:numel(kK),kK,1:1/AcqPointsPerGrad:numel(kK));
%%
res_mm=1.5;

tic
G15mm4L=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
toc
G15mm3L=OutInMultiInnerShotf(FOV,res_mm,VD,3,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G15mm2L=OutInMultiInnerShotf(FOV,res_mm,VD,2,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G15mm1_4L=OutInMultiInnerShotf(FOV,res_mm,VD,1.4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G15mm1L=OutInMultiInnerShotf(FOV,res_mm,VD,1,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);

gPlotTraj_radm(G15mm4L,FOV);MaximizeFig
gPlotTraj_radm(G15mm3L,FOV);MaximizeFig
gPlotTraj_radm(G15mm2L,FOV);MaximizeFig
gPlotTraj_radm(G15mm1_4L,FOV);MaximizeFig
gPlotTraj_radm(G15mm1L,FOV);MaximizeFig
%%
res_mm=1;

tic
G1mm4L=OutInMultiInnerShotf(FOV,res_mm,VD,4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
toc
G1mm3L=OutInMultiInnerShotf(FOV,res_mm,VD,3,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G1mm2L=OutInMultiInnerShotf(FOV,res_mm,VD,2,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G1mm1_4L=OutInMultiInnerShotf(FOV,res_mm,VD,1.4,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
G1mm1L=OutInMultiInnerShotf(FOV,res_mm,VD,1,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);

gPlotTraj_radm(G1mm4L,FOV);MaximizeFig
gPlotTraj_radm(G1mm3L,FOV);MaximizeFig
gPlotTraj_radm(G1mm2L,FOV);MaximizeFig
gPlotTraj_radm(G1mm1_4L,FOV);MaximizeFig
gPlotTraj_radm(G1mm1L,FOV);MaximizeFig
%%

G2mm=cat(2,G2mm4L,G2mm3L,G2mm2L,G2mm1_4L,G2mm1L);
G15mm=cat(2,G15mm4L,G15mm3L,G15mm2L,G15mm1_4L,G15mm1L);
G1mm=cat(2,G1mm4L,G1mm3L,G1mm2L,G1mm1_4L,G1mm1L);
GAll=cat(3,G2mm,G15mm,G1mm);
%%
clear Trajs
%%
Ress=[1.5 2];
VDs=0.9:0.1:1.6;
Loops=0.8:0.1:6;

for r=1:numel(Ress)
    for v=1:numel(VDs)
        for l=1:numel(Loops)
            try
                Trajs{r,v,l}=OutInMultiInnerShotf(FOV,Ress(r),VDs(v),Loops(l),FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
            catch
            end
        end
    end
end

%%

gPlotTraj_radm(Trajs{2,2,1},FOV);


gPlotTraj_radm(Trajs{2,2,50},FOV);
%%
gPlotTraj_radm(Trajs{2,2,50},FOV);

tmp=Trajs{2,2,50};
F=find(abs(tmp)~=0,1,'last');
tmp2=interp1(1:F,tmp(1:F),linspace(1,F,5000)).'*F/5000;

gPlotTraj_radm(tmp2,FOV);






