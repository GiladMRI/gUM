for m=1:3
figure;
for i=1:5
    GTrajaC=GAll(:,i,m);
    g=GTrajaC;
    k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
    kK=k*FOV_mm/1000/2/pi;
    subplot(5,2,i*2-1);
    plot(kK)
    zIdxs=find((real(kK(1:end-1)).*real(kK(2:end)))<0);
    hold on;
    plot(kK(zIdxs),'ro');
    subplot(5,2,i*2);
    SI=sort(imag(kK(zIdxs)));
    dSI=diff(SI);
    plot(dSI)
    
    Crossings=find(  (real(kK(1:end-3)).*real(kK(4:end))<0) & (imag(kK(3:end-1)).*imag(kK(4:end))<0) );
    title(numel(Crossings));
end
MaximizeFig;
end