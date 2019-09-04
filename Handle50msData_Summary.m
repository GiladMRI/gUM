%%
BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';

BaseP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/';
FNs={'meas_MID855_gBP_ASL_SMS_Spi_TI_VD1_ST10_FID51469',...
    'meas_MID857_gBP_ASL_SMS_Spi_TI_VD1_ST11_FID51471',...
    'meas_MID863_gBP_ASL_SMS_Spi_TI_VD1_ST12_FID51477',...
    'meas_MID861_gBP_ASL_SMS_Spi_TI_VD1_ST13_FID51475',...
    'meas_MID865_gBP_ASL_SMS_Spi_TI_VD1_ST14_FID51479',...
    'meas_MID867_gBP_ASL_SMS_Spi_TI_VD1_ST15_FID51481',...
    'meas_MID869_gBP_ASL_SMS_Spi_TI_VD1_ST16_FID51483',...
    'meas_MID871_gBP_ASL_SMS_Spi_TI_VD1_ST17_FID51485',...
    'meas_MID873_gBP_ASL_SMS_Spi_TI_VD1_ST18_FID51487',...
    'meas_MID875_gBP_ASL_SMS_Spi_TI_VD1_ST19_FID51489',...
    'meas_MID883_gBP_ASL_SMS_Spi_TI_VD1_ST10_low_FID51497',...
    'meas_MID885_gBP_ASL_SMS_Spi_TI_VD1_ST13_low_FID51499',...
    'meas_MID887_gBP_ASL_SMS_Spi_TI_VD1_ST15_FID51501'};

RefFldMapPs=[repmat({'/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/meas_MID853_BP_fieldmap_9echos_2mm_Full_FID51467/'},[10,1]);...
    repmat({'/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/meas_MID881_BP_fieldmap_9echos_2mm_Full_low_FID51495/'},[3,1])];
%%
Fi=1;
close all;
for Fi=1:13
FN=FNs{Fi};
RefFldMapP=RefFldMapPs{Fi};

TrgP=[BaseP FN filesep];
%%
load([TrgP 'Traj.mat'],'g','k','kK');
gPlotTraj_radm(g,FOV_mm);MaximizeFig

load([TrgP 'Mg.mat'],'Mg');
load([TrgP 'B0Q2.mat'],'B0Q2');
subplot(2,3,3);
gmontage(Mg);

gprint(0,[BaseP num2str(Fi) '_aTraj.png']);close all
%%
load([TrgP 'Rec1.mat'],'Rec1');
load([TrgP 'Rec1c.mat'],'Rec1c');
load([TrgP 'RecTS_W.mat'],'RecTS_W');
load([TrgP 'RecTS_Wp.mat'],'RecTS_Wp');
%%
figure;
subplot(2,3,1);
gmontage(Rec1);removeTicks;title('Using 1 compressed channel');
subplot(2,3,2);
gmontage(Rec1c);removeTicks;title('Using multichannel');
subplot(2,3,3);
gmontage(B0Q2,[-200 200]);colorbar;removeTicks;title('B_0 Hz');
subplot(2,3,4);
gmontage(RecTS_W);removeTicks;title('Using GRE-based B_0');
subplot(2,3,[5 6]);
gmontage(RecTS_Wp);removeTicks;title('On first and 2nd halves of acquisition');
MaximizeFig;
gprint(0,[BaseP num2str(Fi) '_bBaseRecon.png']);close all
%%
load([TrgP 'RecTS_2LbY.mat'],'RecTS_2LbY');
load([TrgP 'RecTS_2LCbY.mat'],'RecTS_2LCbY');
%%
fgmontage(cat(4,RecTS_2LbY(:,:,1:9:end),RecTS_2LCbY(:,:,1:9:end)));removeTicks;
ylabel('3 components                 5 time-points','FontSize',16);MaximizeFig;
gprint(0,[BaseP num2str(Fi) '_cRecon.png']);close all
%%
end
