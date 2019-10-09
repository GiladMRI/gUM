function g2Shots=OutSingleInnerShotPNSf(FOV,res_mm,VD,nLoops,UseMTK,Gmax_mTm,Smax_Tms,GradDwellTime_ms)
% addpath(genpath('/autofs/space/daisy_002/users/Gilad/mintgrad/'));
%%
% New out in
% Gmax_mTm=38;
% Smax_Tms=160;
% GradDwellTime_ms=10e-3;

% FullROTime_ms=50;
% FOV=200;
nForDesign=10000;
% res_mm=2;
% VD=1.3;
% nLoops=2.1;

Kmax=FOV/(2*res_mm);
kMax_radm=pi*1000/res_mm;
R=(linspace(0,1,nForDesign).^VD)*kMax_radm;
Phi=linspace(0,2*pi*nLoops,nForDesign);
EPhi=exp(1i*Phi);
EPhi=EPhi.*exp(-1i*Phi(end))*1i;
Traj=R.*EPhi;

radm2cm=1/(2*pi*100);

ktoradm=(kMax_radm/Kmax)*radm2cm;

[~,g,~]=minTimeGradient_radm(Traj,[], [], [], Gmax_mTm, Smax_Tms,GradDwellTime_ms);

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;

k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
PhiEnd=angle(k(end));
g=g.*exp(-1i*PhiEnd);
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

% % gPlotTraj_radm(g)
% %%now connect
% % % gEnd=g(end);
% % % gStartNext=-conj(gEnd);
% % % kEnd=k(end);
% % % Nest=100;
% % % Cg=mintimegrad_radm(Nest,gEnd,gStartNext,-2*kEnd,DwellTimeGrad_ms,Gmax_mTm, Smax_Tms);
% % %
% % gEnd=g(end);
% % gStartNext=[];
% % kEnd=k(end);
% % Nest=50;
% % Cg=mintimegrad_radm(Nest,gEnd,gStartNext,-kEnd,GradDwellTime_ms,Gmax_mTm, Smax_Tms);
% %
% AngleT=pi/6;
% Id=find(angle(k)<-AngleT,1,'last');
% dId=numel(k)-Id;
% Id2=Id-dId;
% 
% % figure;
% % plot(k);hold on;
% % plot(k(end),'*');hold on;
% % plot(k(end)*0,'*');
% % plot(k(Id),'*');hold on;
% % plot(conj(k(Id)),'*');hold on;
% % plot(k(Id2),'*');hold on;
% % plot(IxyC,'r')
% %
% CC=[k(Id2) k(Id) k(end) conj(k(Id)) k(end)*0];
% CC(3)=abs(CC(2)).*exp(1i*angle(CC(3)));
% CC(4)=real(CC(4))*0.8+1i*imag(CC(4))*0.8;
% xy=[real(CC(:)) imag(CC(:))].';
% npts=size(xy,2);
% 
% Curve=cscvn(xy);
% Ixy=fnval(Curve,linspace(Curve.breaks(2),Curve.breaks(end),100));
% IxyC=Ixy(1,:)+1i*Ixy(2,:);
% IxyC=IxyC.';
% L=sum(abs(diff(IxyC)));
% dPoint=abs(k(end)-k(end-1));
% Ixy=fnval(Curve,linspace(Curve.breaks(2),Curve.breaks(end),ceil(L/dPoint)));
% IxyC=Ixy(1,:)+1i*Ixy(2,:);
% IxyC=IxyC.';
% 
% % Out=Ixy(1,:)+1i*Ixy(2,:);
% % plot(IxyC,'r')
% %
% kAll1=[k(1:Id); IxyC(2:end); -flipud(IxyC(2:end-1)); -flipud(k(1:Id))];
% % kAll1=[k(1:Id); IxyC(2:end)];
% %%
% % figure;plot(k(1:Id),'b.');hold on
% % plot(IxyC(2:end),'ro');
% % plot(-flipud(IxyC(2:end-1)),'g.')
% % plot(-flipud(k(1:Id)),'k*');
% %%
% % kAll=interp1(1:numel(kAll1),kAll1,1:1:numel(kAll1)).';
% % kAll=interp1(1:numel(kAll1),kAll1,[1:0.5:5 6:numel(kAll1)]).';
% 
% kAll=kAll1;
% % DD=20;
% % kAll=interp1(1:numel(kAll1),kAll1,[1:0.6:DD DD+1:numel(kAll1)-DD-1 numel(kAll1)-DD:0.9:numel(kAll1)]).';
% 
% [~,g2,~]=minTimeGradient_radm(kAll,[], 0, 0, Gmax_mTm, Smax_Tms,GradDwellTime_ms);
% kAll=cumsum([0; g2])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
kAll=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
% kAll=[kAll;0];
% kAll(end)=0;
g2=diff(kAll)/GradDwellTime_ms/TwoPiGammaMHz;
% gPlotTraj_radm(g2) %% 6ms
g20=[0;g2];
% max(abs(diff(g20)))
%%
if(UseMTK)
    Nest=50;
    
    Cg=mintimegrad_radm(Nest,g20(end-Nest),complex(0),-kAll(end-Nest),GradDwellTime_ms,Gmax_mTm, Smax_Tms);
    
    % Cgs=mintimegrad_radm(Nest,complex(0),g20(Nest),kAll(Nest),GradDwellTime_ms,Gmax_mTm, Smax_Tms);
    %
    g20=[0;g2(1:end-Nest);Cg];
end

g2Shots=g20;