VD=1;
FOV=200;
ResV=0.8;

FullROTime_ms=50;
Gmax_mTm=80;
Smax_Tms=170;
GradDwellTime_ms=10e-3;

ma7T=safe_params_ma7T();

FibPhi=(1+sqrt(5))/2;

GoldenAngle=180/FibPhi;
InterShotPhiRad=GoldenAngle*pi/180;

GoldenAnglex=360/(1+FibPhi);
InterShotPhiRad=GoldenAnglex*pi/180;
PNSThresh=95; % percent

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
FOV_mm=FOV;

nRes=numel(ResV);

disp('Base params set');
%
Gms=cell(1,nRes);
Sms=cell(1,nRes);
%%
nInnerShots=2;
TimeForOneSpiral_ms=25;

GradStepsp1=TimeForOneSpiral_ms/GradDwellTime_ms+1;

SlowAccurateTraj=true;
nLoopsPerSpiral=17;
GmsA=OutSingleInnerShotPNSf(FOV,ResV(1),VD,nLoopsPerSpiral,SlowAccurateTraj,Gmax_mTm,Smax_Tms,GradDwellTime_ms);

% GmsA=OutSingleInnerShotPNSf(FOV,ResV(1),VD,17.7,true,Gmax_mTm,135,GradDwellTime_ms);
if(SlowAccurateTraj)
    k=cumsum([0; GmsA]); % mT/m*ms * 2*pi*MHz/T = rad/m
else
    k=[cumsum([0; GmsA]); 0]; % mT/m*ms * 2*pi*MHz/T = rad/m
end
kx=interp1(1:numel(k),k,linspace(1,numel(k),GradStepsp1));
Gms=diff(kx).';

PhiBetweenInnerShots=2*pi/2;

tmp=Row(Gms.*exp(1i*(0:(nInnerShots-1))*PhiBetweenInnerShots)).';
G1=tmp;
gPlotTraj_radm(tmp,FOV,ma7T);MaximizeFig;
subplot(2,3,1);xlabel(nLoopsPerSpiral);
subplot(2,3,2);cla;QuickAnalyzeTrajJumpsf(tmp,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm);
%%
GAll=tmp;
save('GAll1p9mmVD1PAT3.mat','GAll');