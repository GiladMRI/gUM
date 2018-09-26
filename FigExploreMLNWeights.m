VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/';
VP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RegridTry3C2_7TS_S04_Sli01__2018-07-19_01-47-17_train/';
VP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RegridTry3C2_7TS_S04_Sli05__2018-07-18_17-44-37_train/';
VP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RegridTry3C2_7TS_S04_Sli03__2018-07-18_21-46-02_train/'
VP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RegridTry3C2_7TS_S04_Sli10__2018-07-18_00-02-23_train/';
D=dir([VP 'Tra*.mat']);
D=D([D.bytes]>1000);
Q=load([VP D(end).name]);
disp(D(end).name);

G_LossV=Q.G_LossV;
var_list=Q.var_list;
Q=rmfield(Q,{'G_LossV','var_list'});

Q=CombineRIFlds(Q);

Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},65) num2str(size(Q.(SFlds{i})),'% 9d         ')]);
end
%%
X=reshape(Q.gene_GEN_L010_PixelwiseMultC_weightC,[128 128 size(Q.gene_GEN_L010_PixelwiseMultC_weightC,2)]);
% X=reshape(Q.gene_GEN_L006_PixelswiseMult_weightC,[128 128 7]);
Y=X./X(:,:,1);
% ShowAbsAngle(Y)
ShowAbsAngle(X)
subplot(1,2,1);
title('Channels/Time-segments? collapes layer')
xlabel(VP(end-32:end-28))
%%

%%
X=reshape(Q.gene_GEN_L007_PixelwiseMultC_weightC,[128 128 7 8 12]);
%%
[ScrNC,BatchNC,MinNC,LastFN]=GraphOptFromFolderf(VP);
Lines=getLines([VP 'ParamsUsed.txt']);
LinesC=regexp(Lines,' ','split');
St=struct();
for i=1:numel(LinesC)
    St.(LinesC{i}{1})=LinesC{i}{2};
end
% BaseP='/media/a/DATA/11Jul18/RL/meas_MID149_gBP_VD11_U19_G35S155_FID23846/Sli03/';
BaseP=St.RealDataFN(1:end-17);
S=load([BaseP 'SensCC1.mat']);
B=load([BaseP 'B0TS.mat']);
T=load([BaseP 'TSBF_TSC_Sens.mat']);
T2=load([BaseP 'TrajAndRealData.mat']);
%%
F=64/56;
fgmontage(grmss(X,3:5)); hold on;
plot(T2.BARTTrajAct(1,:)*F+64.5,T2.BARTTrajAct(2,:)*F+64.5,'.')
title('Regrid rms weight and trajectory')
axis equal
% plot(-T2.BARTTrajAct(1,:)*F+64.5,-T2.BARTTrajAct(2,:)*F+64.5,'r')
% plot(T2.BARTTrajAct(2,:)*F+64.5,T2.BARTTrajAct(1,:)*F+64.5,'m')
% plot(-T2.BARTTrajAct(2,:)*F+64.5,-T2.BARTTrajAct(1,:)*F+64.5,'y')
% plot([1,2],[1,2],'r')