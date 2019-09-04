VD=1.3;
FOV=200;

FullROTime_ms=50;
Gmax_mTm=35;
Smax_Tms=150;
GradDwellTime_ms=10e-3;
%%
res_mm=2;

G2mm8L=OutInMultiInnerShotf(FOV,res_mm,VD,8,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms);
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






