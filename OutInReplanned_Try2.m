% addpath(genpath('/autofs/space/daisy_002/users/Gilad/mintgrad/'));
%%
% New out in
Gmax_mTm=38;
Smax_Tms=160;

FullROTime_ms=50;
FOV=200;
nForDesign=10000;
VD=1.3;
res_mm=2;
Kmax=FOV/(2*res_mm);
kMax_radm=pi*1000/res_mm;
nLoops=2.1;
R=(linspace(0,1,nForDesign).^VD)*kMax_radm;
Phi=linspace(0,2*pi*nLoops,nForDesign);
EPhi=exp(1i*Phi);
EPhi=EPhi.*exp(-1i*Phi(end))*1i;
Traj=R.*EPhi;

radm2cm=1/(2*pi*100);

ktoradm=(kMax_radm/Kmax)*radm2cm;

GradDwellTime_ms=10e-3;
[~,g,~]=minTimeGradient_radm(Traj,[], [], [], 38, 160,GradDwellTime_ms);

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;

k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
PhiEnd=angle(k(end));
g=g.*exp(-1i*PhiEnd);
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

% gPlotTraj_radm(g)
%% now connect
% gEnd=g(end);
% gStartNext=-conj(gEnd);
% kEnd=k(end);
% Nest=100;
% Cg=mintimegrad_radm(Nest,gEnd,gStartNext,-2*kEnd,DwellTimeGrad_ms,Gmax_mTm, Smax_Tms);
%
gEnd=g(end);
gStartNext=[];
kEnd=k(end);
Nest=50;
Cg=mintimegrad_radm(Nest,gEnd,gStartNext,-kEnd,DwellTimeGrad_ms,Gmax_mTm, Smax_Tms);
%%
% gBoth=[g;Cg];
% kBoth=cumsum([0; gBoth])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
% kAll=[kBoth; flipud(k)];
% figure;plot(kBoth);hold on;plot(flipud(-k))
%%
gAll=[g;Cg;flipud((Cg));flipud(((g)))];
% gPlotTraj_radm(gAll)
kAll=cumsum([0; gAll])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m

% [~,g2,~]=minTimeGradient_radm(interp1(1:numel(kAll),kAll,1:0.1:numel(kAll)),[], [], [], Gmax_mTm, Smax_Tms,GradDwellTime_ms/10);
kAll=[kAll;0];
[~,g2,~]=minTimeGradient_radm(kAll,[], [], 0, Gmax_mTm, Smax_Tms,GradDwellTime_ms);
kAll=cumsum([0; gAll])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
kAll=[kAll;0];
g2=diff(kAll)/GradDwellTime_ms/TwoPiGammaMHz;
% gPlotTraj_radm(g2) %% 6ms
g20=[0;g2;0];
TimePerInnerShot_ms=numel(g20)*GradDwellTime_ms;
%%
nInner=floor(FullROTime_ms/TimePerInnerShot_ms);
clear GCD
for i=1:nInner
    GCDOk(i)=gcd(i,nInner)==1;
end
GCDOk=double(GCDOk);
GCDOk(~GCDOk)=NaN;
GoldenAngle=111.25;
GoldenAngler=GoldenAngle/360;
DiffFromGolden=GCDOk.*(1:nInner)/nInner-GoldenAngler;
DiffFrommGolden=GCDOk.*(1:nInner)/nInner-(1-GoldenAngler);
DiffFrommBoth=min(abs(DiffFromGolden),abs(DiffFrommGolden));
[MinD,MinDI]=min(DiffFrommBoth);
SOrd=mod(MinDI*(0:(nInner-1)),nInner);
Phis=exp(1i*SOrd.*2*pi/nInner);
g2Phis=g20.*Phis;
g2Shots=g2Phis(:);
gPlotTraj_radm(g2Shots)

% gPlotTraj_radm(g20)
