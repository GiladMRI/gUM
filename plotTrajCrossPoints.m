TrajPeaks=quickFindPeaks(diff(abs(Traj),2,2),0.5);
figure;
DataToPlot=grmss(ADataIsL,2);
Mx=gmax(DataToPlot);
for i=1:numel(TrajPeaks)
    subplot(1,numel(TrajPeaks),i);
    plot(abs(Traj)/50000,'k','LineWidth',2);hold on
    plot([TrajPeaks(i) TrajPeaks(i)],gmax(DataToPlot(TrajPeaks(i)+(-100:100),:))*[0 1],'r','LineWidth',2);
    plot(DataToPlot);
    xlim( max(1,min(numel(Traj),TrajPeaks(i)+([-1 1]*100)))   );
    ylim([0,Mx]);
    title([ num2str(AcqTimePoints_us(TrajPeaks(i))/1000,'%.1f') 'ms'],'FontSize',20);
    set(gca,'YTick',[]);
    if(i==1)
        ylabel(BaseFN,'interpreter','None','FontSize',14);
    end
end
MaximizeFig;