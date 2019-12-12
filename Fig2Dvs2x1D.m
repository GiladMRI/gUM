BaseBGP='/media/a/DATA/DataForFig2DFTvs2x1DFT/';
D=dir(BaseBGP);
D=D([D.isdir]);
D=D(3:end);

PC={'M2DFT_HCP128x128ImagesWithPhase__2018-06-25_17-27-31_train', '2D, 64, 0.002';
'M2DFT_HCP128x128ImagesWithPhase__2018-06-25_17-15-06_train','2D, 64, 0.0002';
'M1DFTxy_HCP128x128ImagesWithPhase__2018-06-25_17-25-52_train','2x1D, 128, 0.0002';
'M1DFTxy_HCP128x128ImagesWithPhase__2018-06-25_17-24-32_train','2x1D, 128,0.002';
'M1DFTxy_HCP128x128ImagesWithPhase__2018-06-25_17-22-59_train','2x1D, 64,0.002';
'M1DFTxy_HCP128x128ImagesWithPhase__2018-06-25_17-21-07_train','2x1D, 64,0.0002'};
%%
clear ScrNC BatchNC MinNC
% for i=1:numel(D)
%     [ScrNC{i},BatchNC{i},MinNC{i}]=GraphOptFromFolderf([BaseBGP D(i).name filesep]);
% end
for i=1:size(PC,1)
    [ScrNC{i},BatchNC{i},MinNC{i}]=GraphOptFromFolderf([BaseBGP PC{i,1} filesep]);
end
%%
figure;
Clr={'r','r','b','b','b','b'};
Mrkr={'-','--','--','-','-','--'};
LW=[1,1,2,2,1,1]*2;
for i=1:6
%     plot(BatchNC{i},ScrNC{i},'.-');hold on;
    plot(MinNC{i},ScrNC{i},[Clr{i} '*' Mrkr{i}],'LineWidth',LW(i));hold on;
end
% legend({D.name})
legend(PC{:,2})
xlabel('Minutes');
ylabel('Score');
setYaxis([-0.05 0.25]);
setXaxis([-1 5]);
% %%
% figure;
% subplot(1,2,1);
% plot(BatchN,ScrN,'.-');xlabel('batch #');ylabel('Scr');
% % subplot(1,3,2);
% % plot(MinN,ScrN,'.-');xlabel('Minutes');ylabel('Scr');
% ax1 = gca; % current axes
% ax=axis();
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','origin',...
%     'Color','none');
% 
% line(MinN,ScrN*0,'Parent',ax2,'Color','k')
% setYaxis(ax(3:4));
% xlabel('Minutes')
% subplot(1,2,2);
% ParamsFNX=[BaseGP 'ParamsUsed.txt'];
% lines = getLines(ParamsFNX);
% plot([0 1],[0 1],'.')
% for i=1:numel(lines)
%     text(0.1,1-0.03*i,lines{i},'Interpreter','none');
% end