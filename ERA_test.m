a=load('ERA_epi_r1_VC.mat');
EPI=a.epi_ERA_VC;
a=load('ERA_spi_r1_e1_VC.mat');
SPI=a.spi_e1_ERA_VC;
%
Method='linear';
% Method='spline';
% Method='pchip';
Xs=1:32;
SPI_C=interp1(Xs(1:2:end),SPI(1:2:end),Xs,Method,'extrap');
SPI_T=interp1(Xs(2:2:end),SPI(2:2:end),Xs,Method,'extrap');
EPI_C=interp1(Xs(1:2:end),EPI(1:2:end),Xs,Method,'extrap');
EPI_T=interp1(Xs(2:2:end),EPI(2:2:end),Xs,Method,'extrap');

SPI_A=(SPI_C+SPI_T)/2;
SPI_AN=(SPI_A-SPI_A(3))/SPI_A(3);

EPI_A=(EPI_C+EPI_T)/2;
EPI_AN=(EPI_A-EPI_A(3))/EPI_A(3);

SPI_D=SPI_C-SPI_T;
SPI_DN=SPI_D/SPI_A(3);

EPI_D=EPI_C-EPI_T;
EPI_DN=EPI_D/EPI_A(3);

figure;
subplot(3,2,1);
plot(SPI);hold on;plot(SPI_C,'r');plot(SPI_T,'r');
ylabel(Method,'FontSize',20);
title('Spiral','FontSize',20);
subplot(3,2,2);
plot(EPI);hold on;plot(EPI_C,'r');plot(EPI_T,'r');
title('EPI','FontSize',20);
subplot(3,2,3);
plot(SPI_DN);
title(['Perfusion, max=' num2str(max(SPI_DN),'%.3f')],'FontSize',20);
subplot(3,2,4);
plot(EPI_DN);
title(['Perfusion, max=' num2str(max(EPI_DN),'%.3f')],'FontSize',20);
subplot(3,2,5);
plot(SPI_AN);
title(['BOLD, max=' num2str(max(SPI_AN),'%.3f')],'FontSize',20);
subplot(3,2,6);
plot(EPI_AN);
title(['BOLD, max=' num2str(max(EPI_AN),'%.3f')],'FontSize',20);