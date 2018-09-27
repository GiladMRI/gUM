M=rot90(loadniidata('MLN_mc.nii'));
E=loadniidata('resL1ESPIRiTCCS1A_mc.nii');

ME=M(:,:,:,2:2:end);
MO=M(:,:,:,1:2:end);
DM=ME-MO;
MDM=mean(DM,4);
SDM=std(DM,[],4);
tDM=MDM./SDM;

MM=grmss(M,4);
MMsk=MM>0.2;
MMsk=imfillholesBySlices(MMsk);
DMMsk=imerode(MMsk,ones(7));
fgmontage(MMsk);
fgmontage(DMMsk);

MtMsk=(tDM>0.4) & DMMsk ;
fgmontage(MtMsk);
TwoDM=Reshape4d22d(M,MtMsk);

figure;
rectangle('Position',[-30,-30,300,300],'FaceColor',[0 0 0],'EdgeColor','b')
hold on
plot(2:79,mean(TwoDM(:,2:79),1),'w*-','LineWidth',2);hold on;
plot(2:2:79,mean(TwoDM(:,2:2:79),1),'r-','LineWidth',2);
plot(3:2:79,mean(TwoDM(:,3:2:79),1),'g-','LineWidth',2);
setXaxis([-4 85]);
setYaxis([0.45 0.465]);
% setYaxis([0 0.5]);

print('FigMLN_ASL_Timeline.eps',['-f' num2str(get(gcf,'Number'))],'-deps');
print('FigMLN_ASL_Timeline.png',['-f' num2str(get(gcf,'Number'))],'-dpng');
print('FigMLN_ASL_Timeline.svg',['-f' num2str(get(gcf,'Number'))],'-dsvg');
