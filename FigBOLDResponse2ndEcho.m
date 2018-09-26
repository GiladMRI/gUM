% save('FigBOLDResponse2ndEcho_Data.mat','rmImResRX','TMsk','StimOn')
% load('FigBOLDResponse2ndEcho_Data.mat')
TMsk=rmImResRX>0.8;
TwoD=Reshape4d22d(ImResRX,TMsk);
figure;plot(mean(abs(TwoD),1)-min(mean(abs(TwoD),1)),'LineWidth',2)
hold on
plot(circshift(StimOn,-6,2)*5e-2)
xlabel('Time points (a.u.)');
ylabel('Signal');