function gPlotTraj_radm(g,FOV_mm,PNSsystems,GradPerm)
if(nargin<3)
    PNSsystems=[];
end
if(nargin<4)
    GradPerm=[];
end
if(size(g,1)==2 || size(g,1)==3 || size(g,1)==1)
    g=g.';
end
if(isreal(g))
    g=g(:,1)+1i*(g(:,2));
end

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_ms=10e-3;

k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

kK=k*FOV_mm/1000/2/pi;
figure;
subplot(2,3,1);
plot(kK);
axis square;
axis equal;
ylabel('k');
title(max(abs(kK)));

subplot(2,3,2);
plot(abs(kK));hold on;plot(angle(kK)*10);
setXaxis([1 numel(k)]);
ylabel('k');

subplot(2,3,4);
plot([real(g) imag(g) abs(g)]);
title(['Grad, max ' num2str(max(abs(g)),'%.2f')]);
setXaxis([1 numel(g)]);
ylabel('mT/m');
subplot(2,3,5);
plot([real(s) imag(s) abs(s)]);
title(['Slew, max ' num2str(max(abs(s)),'%.2f')]);
setXaxis([1 numel(s)]);
ylabel('T/m/s');

subplot(2,3,6);
gsafe_plot2(g,GradDwellTime_ms,PNSsystems,GradPerm);

subplot(2,3,3);
ForbiddenFreqsTestf(g,GradDwellTime_ms*1000)