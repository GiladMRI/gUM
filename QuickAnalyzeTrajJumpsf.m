function QuickAnalyzeTrajJumpsf(g,GradDwellTime_ms,TwoPiGammaMHz,FOV_mm)
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
kKA=k*FOV_mm/1000/2/pi;

% subplot(1,2,1);
% plot(kKA)
N=5;
for i=1:N
    kK=kKA*exp(1i*2*pi*(i-1)/N);
%     subplot(1,2,1);
    zIdxs=find((real(kK(1:end-1)).*real(kK(2:end)))<0);
%     hold on;
%     plot(kK(zIdxs)*exp(-1i*2*pi*(i-1)/N),'r.');
%     subplot(1,2,2);
    SI=sort(imag(kK(zIdxs)));
    dSI=diff(SI);
    cSI=(SI(1:end-1)+SI(2:end))/2;
    plot(abs(cSI),dSI,'.');hold on;
end
setXaxis([0 max(abs(kKA))*1.1]);
setYaxis([0 max(abs(dSI))*1.2]);