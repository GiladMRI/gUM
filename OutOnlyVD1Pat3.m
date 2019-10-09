VD=1;
FOV=200;

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

ResV=1.9;
nRes=numel(ResV);

disp('Base params set');
%
Gms=cell(1,nRes);
Sms=cell(1,nRes);
%%
% setting by PNS 95%
% [Gms, Sms]=OutInSingleInnerShotf(FOV,ResV(1),VD,7,false,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% [Gms, Sms]=OutSingleInnerShotf(FOV,ResV(1),VD,17,true,InnerShotTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,GoldenAngle); 
% [Gms, Sms]=OutMultiInnerShotPNSf(FOV,ResV(1),VD,17,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,prisma,PNSThresh,0);

% [g2Shots, OutS]=OutMultiInnerShotPNSf(FOV,res_mm,VD,nLoops,FullROTime_ms,Gmax_mTm,Smax_Tms,GradDwellTime_ms,PNSsystem,PNST,ShotPhi)
GmsA=OutSingleInnerShotPNSf(FOV,ResV(1),VD,18,true,Gmax_mTm,135,GradDwellTime_ms);

GmsA=OutSingleInnerShotPNSf(FOV,ResV(1),VD,17.7,true,Gmax_mTm,135,GradDwellTime_ms);


k=cumsum([0; GmsA]); % mT/m*ms * 2*pi*MHz/T = rad/m
kx=interp1(1:numel(k),k,linspace(1,numel(k),1601));
Gms=diff(kx).';

tmp=Row(Gms.*exp(1i*(0:(3-1))*0)).';
G1=tmp;
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);

tmp=Row(Gms.*exp(1i*(0:(3-1))*(2*pi/3))).';
G2=tmp;
gPlotTraj_radm(tmp,FOV,prisma);MaximizeFig;
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);

GPAT3=[G1 G2];

GAll68x=load('GAll68x.mat');
GAll68cor=load('GAll68cor.mat');

GAll=[GAll68x.GAll GAll68cor.GAll GPAT3];
save('GAll1p9mmVD1PAT3.mat','GAll');