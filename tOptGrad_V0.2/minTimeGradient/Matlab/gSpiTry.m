% return values:
%   C       - reparametrized curve, sampled at T[ms]
%   time    - total time to get to the end
%   g       - gradiet waveform [G/cm]
%   s       - slew rate [G/cm/ms]
%   k       - exact k-space corresponding to gradient g (This function reparametrizes
%             C, then takes a derivative. Numerical errors in the derivative can lead to 
%             deviation.  
%   phi     - Geometry constraints on the amplitude vs. arclength
%   sta     - Solution for the forward ODE
%   stb     - Solution for the backward ODE

riv=1;
rv=0;
ds=[];
Smax=15;
Gmax=3.5;
T=2.5e-3;

VD=1.3;
Maxk=48;
N=200;
nLoops=5;
CSpi=(linspace(0,1,N).^VD)*5  .*exp(1i*linspace(0,2*pi*nLoops,N));
CIn=CSpi;
CIn = [real(CIn(:)), imag(CIn(:)), CIn(:)*0];

[Ck_rv,time_rv,g_rv,s_rv,k_rv, phi_rv, sta_rv, stb_rv]=minTimeGradient(CIn,rv, 0, 0, Gmax, Smax,T,ds,0);
[Ck_riv,time_riv,g_riv,s_riv,k_riv,phi_riv, sta_riv, stb_riv]=minTimeGradient(CIn,riv, 0, 0, Gmax, Smax,T,ds,0);

figure, subplot(2,2,1), plot(k_rv(:,1), k_rv(:,2)); title('k-space'); axis([-6 6 -6 6]);
subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Variant)')
hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
legend('gx', 'gy', 'Location', 'NorthEast');
subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');
%%
riv=0;
ds=[];
Smax=15;
Gmax=3.5;
T=2.5e-3;

% VD=1.3;
% Maxk=48;
% N=200;
% nLoops=5;

VD=2;
Maxk=1600; % rad/meters
res=192;
FOVm=0.192;
RadMeter=pi*res/FOVm % e.g res=128, FOVm=0.2
ICm=RadMeter/(2*pi)/100;
Maxk_Icm=ICm;
Maxk=Maxk_Icm;
nLoops=5;
N=500;
CSpi=(linspace(0,1,N).^VD)*Maxk_Icm  .*exp(1i*linspace(0,2*pi*nLoops,N));
% figure;plot(CSpi)
CIn=CSpi;
% CIn = [real(CIn(:)), imag(CIn(:)), CIn(:)*0];

[Ck_riv,time_riv,g_riv,s_riv,k_riv,phi_riv, sta_riv, stb_riv]=minTimeGradient(CIn,riv, 0, 0, Gmax, Smax,T,ds,0);

% L = max(length(s_riv), length(s_rv));
L = length(s_riv);

% figure, subplot(2,2,1), plot(k_riv(:,1), k_riv(:,2)); title('k-space'); axis([-6 6 -6 6]);
% subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms')
% plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
% plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 

figure, subplot(2,2,1), plot(k_riv); title('k-space'); axis([-1 1 -1 1]*Maxk);
subplot(2,2,2), plot(real(g_riv)); axis([0,L,-4.5,4.5]); title('gradient waveforms')
hold on, subplot(2,2,2), plot(imag(g_riv), 'r');
legend('gx', 'gy', 'Location', 'NorthEast');
subplot(2,2,3)
plot(abs(g_riv), 'r');  axis([0 L 0 6]);
subplot(2,2,4), 
plot(abs(s_riv), 'r'); 
title('slew-rate magnitude');  axis([0 L 0 20]);

% COut=Ck_riv(:,1)+1i*Ck_riv(:,2);
%%

KP=20;
KPNew=20;
C1=kTraj1(1:KP,1)+1i*kTraj1(1:KP,2);
U=unwrap(angle(C1));
CC=linspace(0,abs(C1(end)),KPNew).*exp(1i*linspace(0,U(end),KPNew));

figure;plot(kTraj1(1:KP,1),kTraj1(1:KP,2))
hold on;
plot(real(CC),imag(CC),'r')
%%
% kTraj=kTraj1(1:20,:);

% clear TT
% for i=1:2
%     TT(:,i)=interp1(1:KP,kTraj1(1:KP,i),1:0.5:KP);
% end
TT=[real(CC);imag(CC)].';
kTraj=[TT; kTraj1((KP+1):200,:)];

GradBuf=[0 0;diff(kTraj(:,:))*MaxGrad]/2.3833;

BARTTrajx=kTraj.'*FOVx/1000/2/pi;
nTrajG=size(BARTTrajx,2);
C=BARTTrajx(1,1:size(BARTTrajx,2)/2)+1i*BARTTrajx(2,1:size(BARTTrajx,2)/2);

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajxQ=kTrajQ.'*FOVx/1000/2/pi;
nTraj=size(BARTTrajxQ,2);
Traj=BARTTrajxQ;
%
figure;
subplot(2,2,1);
plot3(1:size(BARTTrajx,2),BARTTrajx(1,:),BARTTrajx(2,:))
%
HnTraj=floor(nTrajG/nInnerInterleaves);
CLRs='krbbcmkrgkrbbcmkrg';
subplot(2,2,2);plot(BARTTrajx(1,1:HnTraj),BARTTrajx(2,1:HnTraj),[CLRs(1) '-'])
hold on
for i=2:nInnerInterleaves
    plot(BARTTrajx(1,HnTraj*(i-1)+(1:HnTraj)),BARTTrajx(2,HnTraj*(i-1)+(1:HnTraj)),[CLRs(i) '-'])
end
axis equal
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
%
subplot(2,2,3);
Ts=(0:(size(GradBuf,1)-1))*10;
plot(Ts,GradBuf*MaxGrad*1000,'.-');title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
% setXaxis([1 nTraj]);CBARTTrajx=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel('time (ms)');
ylabel('mT/m');
%
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(Ts(1:end-1),SlewBuf*100);
MaxSlew=max(max(abs(SlewBuf((20+ExtraDelay):(HnTraj-1),:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
% setXaxis([1 nTraj]);dTraj(:,1)
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel(['time (ms), ' num2str(nTraj*TimePerAcqPoint_us/1000,'%.1f') 'ms']);
ylabel('mT/m/s');

CBARTTrajx=BARTTrajx(1,:)+1i*BARTTrajx(2,:);