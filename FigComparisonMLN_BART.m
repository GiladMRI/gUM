NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S1__2018-07-08_13-08-03_train/';
SliI=1;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S2__2018-07-07_19-49-15_train/';
SliI=2;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S3__2018-07-07_09-50-44_train/';
SliI=3;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S4__2018-07-07_23-50-30_train/';
SliI=4;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S5__2018-07-07_17-41-32_train/';
SliI=5;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S6__2018-07-08_15-14-42_train/';
SliI=6;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S9__2018-07-09_00-48-37_train/';
SliI=9;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S11__2018-07-08_17-35-59_train/';
SliI=11;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S12__2018-07-08_22-37-11_train/';
SliI=12;
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S10__2018-07-09_08-39-10_train/';
SliI=10; % Problematic run again
NNP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_AK_S7__2018-07-09_12-00-36_train/';
SliI=7;

[ScrNC,BatchNC,MinNC,LastFN]=GraphOptFromFolderf([NNP]);

X=imread([NNP LastFN]);

NNRes=double(rot90(X(1:128,128*5+(1:128),1)))/255;
%

BaseBP='/media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/';
X=load([BaseBP 'BARTRecon_AllS_NoB0_W1e-05.mat']);

BARTRes=double(rot90(X.RecS(:,:,SliI)));
BARTRes=BARTRes*grmss(NNRes)/grmss(BARTRes);
%
All=cat(3,BARTRes,NNRes);
fgmontage(All,[0 0.5])
title(['BART                                                                                    NN: ' LastFN],'Interpreter','None');
ylabel(['Sli #' num2str(SliI)]);
%%
BaseTP='/media/a/DATA/TrainedNetsOnAK/';
for i=1:12
    D=dir([BaseTP 'RegridTry3C2_7TS_AK_S' num2str(i) '_*train']);
    SliP{i}=D(1).name;
end
%%
ParamsSDefaults=gLines2Struct(getLines('/media/a/DATA/TrainedNetsOnAK/RegridTry3C2_7TS_AK_S1__2018-07-08_13-08-03_train/ParamsUsed.txt'));
%%
% TBaseP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';
%
% DataH=147; %size(NMap,1);
% DataW=147; %size(NMap,2);
% DataCh=192; %size(NMap,3)*nChToUseInNN*2;
% LabelsH=1;%Sz2(1);
% LabelsW=1;%Sz2(2);
% LabelsCh=1; %2;
% RealDataFN='x';
% NMapFN='/media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/NMapIndCB0X.mat'; %'x';
%
% ParamsSDefaults=struct('DataH',DataH,'DataW',DataW,'channelsIn',DataCh,'LabelsH',LabelsH,'LabelsW',LabelsW,'channelsOut',LabelsCh,...
%   'aDataH' 147
% aDataW 147
%   'dataset',TFDataP,'learning_rate_start',0.002,...
%   'learning_rate_half_life',30,... % in minutes if <1000
%   'summary_period',0.5,'checkpoint_period',20,...
%   'MapSize',3,'train_time',120,'batch_size',16,'NumFeatPerChannel',2,'NumTotalFeat',64,...
%   'WL1_Lambda',0,'WL2_Lambda',0,...
%   'QuickFailureTimeM',3,'QuickFailureThresh',0.3,'DiscStartMinute',500,...
%   'ShowRealData',1,'CmplxBias',0,...
%   'InputMode','RegridTry1',...
%   'NetMode','RegridTry1C2_TS',...
%   'SessionNameBase','RegridTry1C2_TS',...
%   'ImgMode','Cmplx',...
%   'nTimeSegments',7,...
%   'UseSharedWightesInRelaxedFT',1,...
%   'WPhaseOnly',0.001,...
%   'NMAP_FN',NMapFN,...
%   'RealDataFN',RealDataFN,...
%   'LoadAndRunOnData',1,...
%   'LoadAndRunOnData_checkpointP','x',...
%   'LoadAndRunOnData_OutP','x');
%
% ParamsS=ParamsSDefaults;
% Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
%%
for i=1:12
    ParamsS=ParamsSDefaults;
    ParamsS.dataset='/media/a/DATA/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/Sli01/TF/';
    ParamsS.LoadAndRunOnData=1;
    ParamsS.LoadAndRunOnData_checkpointP=[BaseTP SliP{i}(1:end-5) 'checkpoint'];
    ParamsS.LoadAndRunOnData_Prefix=['/media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/RealData/Sli' num2str(i) '_r'];
    ParamsS.LoadAndRunOnData_OutP=['/media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/RealDataOut/Sli' num2str(i,'%02d') filesep];
    Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
    disp(['ok ' num2str(i)]);
end
%%
% BaseP='/media/a/DATA/180628_AK/';
% ScanP='/media/a/DATA/180628_AK/';
% RefFldMapP='/media/a/DATA/180628_AK/meas_MID265_BP_fieldmap_5echosX_FID22460/';
% BaseFN='meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439';

% BaseP='/media/a/DATA/11Jul18/RL/';
% ScanP='/media/a/DATA/11Jul18/RL/';
% RefFldMapP='/media/a/DATA/11Jul18/RL/meas_MID149_gBP_VD11_U19_G35S155_FID23846/';
% BaseFN='meas_MID149_gBP_VD11_U19_G35S155_FID23846';
CurRealDataOutP='/media/a/DATA/ASLSubjData/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RealDataOut/';
nReps=80;
nSlices=12;
% CurRealDataOutP=[ScanP BaseFN filesep 'RealDataOut' filesep];
% CurRealDataOutP='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944/RealDataOut/';
% /media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439
for SliI=1:nSlices
    for r=1:nReps
        %     M=load(['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/Out/OnRealData' num2str(r,'%02d') '.mat']);
        M=load([CurRealDataOutP 'Sli' num2str(SliI,'%02d') filesep 'OnRealData' num2str(r,'%02d') '.mat']);
        
        MC=M.x(:,:,:,1)+1i*M.x(:,:,:,2);
        MC=permute(MC,[2 3 1]);
        
        AMCS(:,:,SliI,r)=double(MC(:,:,1));
    end
    disp(SliI);
end
disp('ok')
%%
BaseOutLoc='/media/a/DATA/';
% BaseOutLoc=ScanP;
CurOurP=[BaseOutLoc BaseFN filesep];

im_resS=load([CurOurP 'im_resS.mat'],'im_resS');
im_resS=im_resS.im_resS;

Lambda=1e-6;
T2SCompStr='';
RecBART=load([ScanP BaseFN filesep 'BARTRecon_NoB0_W' num2str(Lambda) T2SCompStr '.mat'],'RecBART');
RecBART=RecBART.RecBART;
%%
RecBARTN=RecBART*grmss(AMCS)/grmss(RecBART);
im_resSN=im_resS*grmss(AMCS(:,:,:,1))/grmss(im_resS);
%%
RecBARTNP=rot90(RecBARTN);
AMCSP=rot90(AMCS);
im_resSP=rot90(im_resSN);
%%

FirstEcho=load([RefFldMapP 'FirstEcho.mat']);
FirstEcho=FirstEcho.FirstEcho;

FirstEcho=gflip(FirstEcho(:,:,:,6+(1:nSlices)),1:2);
Mg=grmss(FirstEcho,3);
disp('ok');

MgN=imresizeBySlices(Mg,gsize(AMCS));
MgN=MgN*grmss(AMCS(:,:,:,1))/grmss(MgN);
MgP=rot90(MgN);
%%
SliI=11;
fgmontage(cat(3,im_resSP(:,:,SliI),RecBARTNP(:,:,SliI,1),AMCSP(:,:,SliI,1)),'Size',[1 3]);
title([PadStringWithBlanks('ESPIRIT with B0',80) PadStringWithBlanks('BART No B0',80) 'MLN']);
ylabel(['Slice #' num2str(SliI)]);
%% With Mg
SliI=6;
fgmontage(cat(3,im_resSP(:,:,SliI),RecBARTNP(:,:,SliI,1),AMCSP(:,:,SliI,1),MgP(:,:,SliI)),'Size',[1 4]);
title([PadStringWithBlanks('ESPIRIT with B0',60) PadStringWithBlanks('BART No B0',60) PadStringWithBlanks('MLN',60) 'ME']);
ylabel(['Slice #' num2str(SliI)]);
%% Noise
AMCSE=AMCSP(:,:,:,6:2:end-5);
RecBARTNPE=RecBARTNP(:,:,:,6:2:end-5);

SAMCSE=std(AMCSP,[],4);
SRecBARTNPE=std(RecBARTNPE,[],4);
fgmontage(cat(3,SRecBARTNPE(:,:,SliI),SAMCSE(:,:,SliI)),'Size',[1 2]);
title([PadStringWithBlanks('BART No B0',120) 'MLN']);
xlabel('Real Noise (From unlabeled part) amplification');
ylabel(['Slice #' num2str(SliI)]);

%%
MAMCSE=abs(mean(AMCSP(:,:,:,6:2:end-5),4));
MAMCSO=abs(mean(AMCSP(:,:,:,5:2:end-5),4));

MRecBARTNPE=abs(mean(RecBARTNP(:,:,:,6:2:end-5),4));
MRecBARTNPO=abs(mean(RecBARTNP(:,:,:,5:2:end-5),4));

PerfNN=squeeze(MAMCSE-MAMCSO);
PerfBART=squeeze(MRecBARTNPE-MRecBARTNPO);

fgmontage(cat(3,PerfBART(:,:,SliI),PerfNN(:,:,SliI,1)),[0 3e-3],'Size',[1 2]);
%% Time averaged
SliI=6;
fgmontage(cat(3,MRecBARTNPE(:,:,SliI),MAMCSE(:,:,SliI)),'Size',[1 2]);
title([PadStringWithBlanks('BART No B0',80) 'MLN']);
ylabel(['Slice #' num2str(SliI)]);
xlabel('Averaged over time');
%%
AMCSPS=SmoothBySlices(AMCSP,5,0.6);
% AMCSPS=AMCSP;
MAMCSES=abs(mean(AMCSPS(:,:,:,6:2:end-5),4));
MAMCSOS=abs(mean(AMCSPS(:,:,:,5:2:end-5),4));
PerfNNS=squeeze(MAMCSES-MAMCSOS)/1.2;

MAMCSS=abs(mean(AMCSPS,4));

MRecBARTNPE=abs(mean(RecBARTNPE,4));

Fore=PerfNNS;
Back=MAMCSS;
PStr='NN Smoothed';
ClimF=[10e-5 35e-4];
ClimB=[0 4e-1];

% Fore=PerfNN;
% Back=MAMCS;
% PStr='NN';

% Fore=PerfBART;
% Back=MRecBARTNPE;
% PStr='BART';
% ClimF=[0.0002 4e-3];
% ClimB=[0 4e-1];
%%
WhichSlicesShow=1:12;

Fore(Back<0.1)=0;
A=(min(max(Fore,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
for s=1:numel(WhichSlicesShow)
    RGB(:,:,:,s)=ind2rgb(round(squeeze(A(:,:,s))*255)+1,parula(256));
end
Msk=permute(Fore(:,:,WhichSlicesShow)>ClimF(1),[1 2 4 3]);

B=(min(max(Back,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+permute(B(:,:,WhichSlicesShow),[1 2 4 3]).*(1-Msk);
% figure;imshow(X)
% X2=CombineDims(X,[4 2]);
X1=PartitionDim(X,4,3);
X2=CombineDims(X1,[4 2]);
X2=CombineDims(X2,[4 1]);
figure;imshow(X2);title(PStr);
%%
ResNNS=X2;
X3=cat(1,ResNNS,ResBART);
figure;imshow(X3);title('top - smoothed NN, Bottom - BART');

%%

AllRes3=repmat(AllRes,[1 1 1 1 3]);
AllRes3(:,:,3,2,:)=X*ClimB(2);

AllResI=CombineDims(AllRes3,[3 2]);
AllResI=CombineDims(AllResI,[3 1]);
figure;imshow(abs(AllResI)/ClimB(2))

FS=11;
Clr=[1 1 1];
text(10,10,'BART 1 map','Color',Clr,'FontSize',FS)
text(128*1+0,10,'BART 1 map','Color',Clr,'FontSize',FS)
text(10,128*1+10,'ESPIRIT 2 maps, 15TS','Color',Clr,'FontSize',FS)
text(128*1+0,128*1+10,'ESPIRIT 2 maps, 15TS, CC->13','Color',Clr,'FontSize',FS)

text(128*2+0,10,'Linear net, 7TS','Color',Clr,'FontSize',FS)
text(128*2+0,128*1+10,'Perfusion, 7TS','Color',Clr,'FontSize',FS)
%%
BaseP='/media/a/DATA/11Jul18/RL/';
ScanP='/media/a/DATA/11Jul18/RL/';
RefFldMapP='/media/a/DATA/11Jul18/RL/meas_MID149_gBP_VD11_U19_G35S155_FID23846/';
BaseFN='meas_MID149_gBP_VD11_U19_G35S155_FID23846';
nReps=80;
nSlices=12;
CurRealDataOutP=[ScanP BaseFN filesep 'RealDataOut' filesep];
% CurRealDataOutP='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944/RealDataOut/';
% /media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439
WhichSlices=3;
for SliI=WhichSlices
    for r=1:nReps
        %     M=load(['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/Out/OnRealData' num2str(r,'%02d') '.mat']);
        M=load([CurRealDataOutP 'Sli' num2str(SliI,'%02d') filesep 'OnRealData' num2str(r,'%02d') '.mat']);
        
        MC=M.x(:,:,:,1)+1i*M.x(:,:,:,2);
        MC=permute(MC,[2 3 1]);
        
        AMCS(:,:,SliI,r)=double(MC(:,:,1));
    end
    disp(SliI);
end
disp('ok')
%%
AMCSP=rot90(AMCS);
AMCSE=AMCSP(:,:,:,6:2:end-5);
SAMCSE=std(AMCSP,[],4);
fgmontage(SAMCSE)
MAMCSE=abs(mean(AMCSP(:,:,:,6:2:end-5),4));
MAMCSO=abs(mean(AMCSP(:,:,:,5:2:end-5),4));
PerfNN=squeeze(MAMCSE-MAMCSO);
MAMCS=abs(mean(AMCSP,4));

AMCSPS=SmoothBySlices(AMCSP,5,0.6);
% AMCSPS=AMCSP;
MAMCSES=abs(mean(AMCSPS(:,:,:,6:2:end-5),4));
MAMCSOS=abs(mean(AMCSPS(:,:,:,5:2:end-5),4));
PerfNNS=squeeze(MAMCSES-MAMCSOS)/1.2;

MAMCSS=abs(mean(AMCSPS,4));
%%
% Fore=PerfNNS;
% Back=MAMCSS;
% PStr='NN Smoothed';

Fore=PerfNN;
Back=MAMCS;
PStr='MLN';

ClimF=[10e-4 8e-3];
ClimB=[0 8e-1];
BThreshForF=0.1;
WhichSlicesShow=1:12;

Fore(Back<BThreshForF)=0;
A=(min(max(Fore,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
for s=1:numel(WhichSlicesShow)
    RGB(:,:,:,s)=ind2rgb(round(squeeze(A(:,:,WhichSlicesShow(s)))*255)+1,parula(256));
end
Msk=permute(Fore(:,:,WhichSlicesShow)>ClimF(1),[1 2 4 3]);

B=(min(max(Back,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+permute(B(:,:,WhichSlicesShow),[1 2 4 3]).*(1-Msk);

% figure;imshow(X)
%
Z=PartitionDim(X,4,3);
Z=CombineDims(Z,[4 2]);
Z=CombineDims(Z,[4 1]);
figure;imshow(Z);title(PStr)
%% Noise comparison
ESPIRITRec=loadniidata('resL1ESPIRiTCCS1A.nii');
MLNRec=abs(AMCS);
%%
ESPIRITRec=ESPIRITRec*grmss(MLNRec)/grmss(ESPIRITRec);
MLNRecE=MLNRec(:,:,:,2:2:end);
ESPIRITRecE=ESPIRITRec(:,:,:,2:2:end);
SMLNRecE=std(MLNRecE(:,:,:,5:35),[],4);
SESPIRITRecE=std(ESPIRITRecE(:,:,:,5:35),[],4);
Rng=[0 10e-3];
figure;
subplot(1,2,1);
gmontage(rot90(SMLNRecE),Rng*1);colorbar;title('MLN')
xlabel('Std, a.u')
subplot(1,2,2);
gmontage(rot90(SESPIRITRecE),Rng);colorbar;title('ESPIRIT')
xlabel('Std, a.u')
%%
T=1e-4;
Edges=linspace(T,0.05,300);
NSESPIRITRecE=histc(SESPIRITRecE(SESPIRITRecE>T),Edges);
NSMLNRecE=histc(SMLNRecE(SMLNRecE>T),Edges);
figure;plot(Edges,[NSESPIRITRecE NSMLNRecE],'LineWidth',2);
legend({'ESPIRIT','MLN'},'FontSize',14)
xlabel('a.u.');
ylabel('Count')
%%


%%
close all
test = zeros(size(X));

test(5:8,5:8,5:7)=50;
%large cube at higher
test(7:11,7:11,8:11)=60;

test= (A~=0)*50;
[X,Y,Z] = meshgrid(1:size(A,1),1:size(A,2),1:size(A,3));
p= patch(isosurface(X,Y,Z,test,1));
% figure;
isonormals(X,Y,Z,test,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis([1,size(A,1),1,size(A,2),1,size(A,3)])
camlight 
lighting gouraud
grid on