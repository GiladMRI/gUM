MLNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_GL_S3__2018-07-10_19-13-55_train/';
ScanP='/media/a/DATA/13May18/Me/';
BaseFN='meas_MID409_gBP_VD11_U19_7ADCs_FID17798';
SliI=3;
YLbl='Sli03';

MLNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/';
ScanP='/media/a/DATA/11Jul18/RL/';
BaseFN='meas_MID149_gBP_VD11_U19_G35S155_FID23846';
SliI=3;
YLbl='Sli03';
%%
D=dir([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0*_CC.mat']);
ESPResCC=load([D(1).folder filesep D(1).name]);ESPResCC=ESPResCC.resL1ESPIRiTCC1;
ESPResCC=double(ESPResCC);

D=dir([ScanP BaseFN filesep YLbl '_L1ESPIRiT_B0*.mat']);
D=D(~strhas({D.name},'CC'));
ESPRes=load([D.folder filesep D.name]);ESPRes=ESPRes.resL1ESPIRiT1;
ESPRes=double(ESPRes);

D=dir([ScanP BaseFN filesep YLbl '_*].mat']);
SparseMRIRes=load([D.folder filesep D.name]);SparseMRIRes=SparseMRIRes.im_res;
SparseMRIRes=double(SparseMRIRes);

D=dir([ScanP BaseFN filesep 'BARTRecon_AllS_NoB0_W*.mat']);
BARTRes=load([D(1).folder filesep D(1).name]);BARTRes=BARTRes.RecS(:,:,SliI);
BARTRes=double(BARTRes);

D=dir([ScanP BaseFN filesep 'BARTRecon2Maps_AllS_NoB0_W*.mat']);
BART2MapsRes=load([D.folder filesep D.name]);BART2MapsRes=BART2MapsRes.RecSMM(:,:,SliI);
BART2MapsRes=double(BART2MapsRes);

[ScrNC,BatchNC,MinNC,LastFN]=GraphOptFromFolderf(MLNP);
X=imread([MLNP LastFN]);
MLNRes=X(1:128,128*5+(1:128),1);
MLNRes=double(MLNRes);

ESPRes=ESPRes*grmss(MLNRes)/grmss(ESPRes);
ESPResCC=ESPResCC*grmss(MLNRes)/grmss(ESPResCC);
SparseMRIRes=SparseMRIRes*grmss(MLNRes)/grmss(SparseMRIRes);
BARTRes=BARTRes*grmss(MLNRes)/grmss(BARTRes);
BART2MapsRes=BART2MapsRes*grmss(MLNRes)/grmss(BART2MapsRes);
disp('ok');
%%
clear AllRes
AllRes(:,:,:,1)=cat(3,BARTRes,BART2MapsRes);
AllRes(:,:,:,2)=cat(3,ESPRes, ESPResCC);
AllRes(:,:,3,1)=MLNRes;
AllRes(:,:,3,2)=SparseMRIRes;
AllRes=gflip(permute(AllRes,[2 1 3 4]),1);
fgmontage(AllRes)
title('');xlabel([ScanP ' ' BaseFN ' ' YLbl],'Interpreter','None');

XLim=get(gca,'XLim');
YLim=get(gca,'YLim');
NNNx=floor(XLim(2))/3;
NNNy=floor(YLim(2))/2;
gx=20;
gy=20;
FS=11;
Clr=[1 1 1];
text(gx,gy,'BART','Color',Clr,'FontSize',FS)
text(NNNx*1+0,gy,'BART 2 Maps','Color',Clr,'FontSize',FS)
text(NNNx*2+0,gy,'MLN','Color','red','FontSize',FS)

text(gx,NNNy+gy,'L_1 CG-ESPIRIT','Color',Clr,'FontSize',FS)
text(NNNx*1+0,NNNy*1+gy,'L_1 CG-ESPIRIT+CC, \lambda=10^{-4}','Color',Clr,'FontSize',FS)
text(NNNx*2+0,NNNy*1+gy,'SparseMRI','Color',Clr,'FontSize',FS)
%%
% gprint(g