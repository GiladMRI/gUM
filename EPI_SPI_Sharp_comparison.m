a=load('spi_perf_sli6.mat');
SPI=a.spi_perf_sli6;
a=load('epi_perf_sli6.mat');
EPI=a.epi_perf_sli6;
B=cat(4,SPI,EPI);
RGB=squeeze(B(:,:,1,:));
RGB(1,1,3)=0;
B=max(B,0);
for i=1:2
    for j=1:2
        K(:,:,i,j)=autocorr2d(B(:,:,i,j));
    end
end
close all
fgmontage(max(B,0))
figure;
CLRs='rg';
Mark='-:';
subplot(1,2,1);
plot(squeeze(K(49,:,1,1)))
k=1;
for i=1:2
    for j=1:2
        Cur=squeeze(K(49,:,i,j));
        Cur=Cur./max(Cur);
        h(k)=plot(Cur,[Mark(i) CLRs(j)]);
        hold on;
        k=k+1;
    end
end
legend(h,{'SPI band 1','EPI band 1','SPI band 2','EPI band 2'},'FontSize',20);
title('Axis 1','FontSize',20)
subplot(1,2,2);
plot(squeeze(K(49,:,1,1)))
k=1;
for i=1:2
    for j=1:2
        Cur=squeeze(K(:,49,i,j));
        Cur=Cur./max(Cur);
        h(k)=plot(Cur,[Mark(i) CLRs(j)]);
        hold on;
        k=k+1;
    end
end
legend(h,{'SPI band 1','EPI band 1','SPI band 2','EPI band 2'},'FontSize',20);
title('Axis 2','FontSize',20)