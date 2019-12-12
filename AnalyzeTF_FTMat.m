% A=load('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_dataNeighborhoodR__2018-06-07_17-43-50_train/TrainSummary_028420.mat');
% Flds=fieldnames(A);
% 
% X=A.gene_GEN_L004_M2D_MC_weightR_0+1i*A.gene_GEN_L004_M2D_MC_weightI_0;
% Y=A.gene_GEN_L005_M2D_MC_weightR_0+1i*A.gene_GEN_L005_M2D_MC_weightI_0;
% %%
% fgmontage(fft1cg(X,2))
% %%
% fgmontage(fft1cg(Y,2))


%%
% VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS2_dataNeighborhoodRCB0__2018-06-08_16-17-56_train/';
% VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS__2018-06-27_17-45-22_train/'

% VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS__2018-06-27_20-12-01_train/';

VP='/media/a/DATA/First interesting NN regridding results/RegridTry1C2_TS2_dataNeighborhoodRCB0__2018-06-09_19-44-17_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS2__2018-06-28_17-35-13_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS2__2018-06-28_17-58-43_train/';

VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS2__2018-06-28_18-35-53_train/';
%%
BaseTFRes='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';
Prefix='RegridTry1C2_TS';
% Prefix='RegridTry1C2_TS2';
DB=dir([BaseTFRes Prefix '*']);
LastDir=DB(end).name;

VP=[BaseTFRes LastDir filesep];
%% Maps
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/1DFTxyCMaps_HCP128x128ImagesWithPhase__2018-06-25_19-35-05_train/';
%% SCC
% 1ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-33-46_train/';
% 2ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-37-15_train/';
% 3ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-41-18_train/';
% 5ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-28-30_train/';
%% GCC
% 2ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC2__2018-06-29_16-49-45_train/';
% 3ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC3__2018-06-29_17-00-15_train/';
%% YCC 2ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_YCC2__2018-06-29_17-16-12_train/';
% LR 0.0002
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_YCC2__2018-06-29_17-36-15_train/';
%%
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASHTry1_GCCF_dataNeighborhoodSMASH6F__2018-06-11_16-55-48_train/';
%% TSB
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_TSB__2018-07-02_15-00-13_train/';
% VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_TSB__2018-07-02_15-48-58_train/';
%%
D=dir([VP 'Tra*.mat']);
D=D([D.bytes]>1000);
Q=load([VP D(end).name]);
disp(D(end).name);

G_LossV=Q.G_LossV;
var_list=Q.var_list;
Q=rmfield(Q,{'G_LossV','var_list'});

%%
Q=CombineRIFlds(Q);

Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},65) num2str(size(Q.(SFlds{i})),'% 9d         ')]);
end
%%
aaTSB=Q.gene_GEN_L007_PixelswiseMultC_weightC;
aaRegrid=reshape(Q.gene_GEN_L015_PixelswiseMultC_weightC,[131,131,96]);
% aaRegrid=permute(aaRegrid,[2 1 3]);

aaFT1=Q.gene_GEN_L018_add_Mult2DMCyCSharedOverFeat_weightC;
aaFT2=Q.gene_GEN_L019_add_Mult2DMCxCSharedOverFeat_weightC;
aaTSC=Q.gene_GEN_L020_PixelswiseMultC_weightC;
aaTSC=reshape(aaTSC,[128 128 7]);
% aaTSC=permute(aaTSC,[2 1 3]);

Start=tmp;
AfterTSB=Start.*permute(aaTSB,[1 3 2]);

clear Ordered
for i=1:7
    CurData=Row(AfterTSB(:,:,i));
    Ordered(:,:,:,i)=CurData(NMapCX);
end
% Ordered=permute(Ordered,[2 1 3 4]);
AfterRegrid=squeeze(sum(Ordered.*aaRegrid,3));
AfterFT1=sum(AfterRegrid.*permute(aaFT1,[1 3 4 2]));
AfterFT2=sum(AfterFT1.*permute(aaFT2,[5 1 3 4 2]));
AfterFT2x=permute(AfterFT2,[4 5 3 1 2]);
AfterTSC=squeeze(sum(AfterFT2x.*aaTSC,3));
ShowAbsAngle(AfterFT2x)
%% CC
MapsForTF=load('/media/a/H1/maps128x128x8.mat');
X=squeeze(Q.gene_GEN_L004_C2D_weightC);
OutCh=permute(MultMatTensor(X.',permute(MapsForTF.maps,[3 1 2])),[2 3 1]);
% ShowAbsAngle(OutCh,[],'Size',[1 3])
%% GCC
X=squeeze(Q.gene_GEN_L004_einsum_weightC);
% Y=permute(X,[1 4 2 3]);
Y=permute(X,[4 1 2 3]);
Z=Y.*MapsForTF.maps;
F=squeeze(sum(Z,3));
ShowAbsAngle(F)
%% YCC
X=squeeze(Q.gene_GEN_L005_einsum_weightC);
% Y=permute(X,[1 4 2 3]);
% Y=permute(X,[4 1 2 3]);
Z=X.*MapsForTF.maps(1:2:end,:,:);
F=squeeze(sum(Z,3));
% ShowAbsAngle(F)
R=F(:,:,1)./F(:,:,2);
% ShowAbsAngle(R)
R(~MapsForTF.Msk(1:2:end,:))=NaN;
U=unwrap(angle(R(:,40:92)));
U=U-U(32,:);
figure;plot(U);hold on;plot([0 64],[pi -pi],'k','LineWidth',2)
%% Maps:
X=rot90(reshape(Q.gene_GEN_L005_MapsForMatPixelswiseMultC_weightC,[128 128 8]));
ShowAbsAngle(X,[],'Size',[2 4])
%%
X=reshape(Q.gene_GEN_L006_PixelswiseMultC_weightC,[128 128 7]);
% X=reshape(Q.gene_GEN_L006_PixelswiseMult_weightC,[128 128 7]);
Y=X./X(:,:,1);
ShowAbsAngle(Y)
%%
X=reshape(Q.gene_GEN_L003_PixelswiseMult_weightC,[131 131 12 8 7]);
% A=reshape(Q.gene_GEN_L003_PixelswiseMultC_weightC,[131 131 
%%
% size(NMapCR)
nTS=7;
PWK=double(Q.gene_GEN_L003_PixelswiseMult_weightR_0 + 1i*Q.gene_GEN_L003_PixelswiseMult_weightI_0);
PWK=reshape(PWK,[131 131 96 nTS]);
PWKB=squeeze(double(Q.gene_GEN_L003_PixelswiseMult_bias_0(:,:,:,:,1)+1i*Q.gene_GEN_L003_PixelswiseMult_bias_0(:,:,:,:,2)));
FTX=Q.gene_GEN_L004_M2D_MC_weightR_0+1i*Q.gene_GEN_L004_M2D_MC_weightI_0;
FTY=Q.gene_GEN_L005_M2D_MC_weightR_0+1i*Q.gene_GEN_L005_M2D_MC_weightI_0;

PWI=double(Q.gene_GEN_L006_PixelswiseMult_weightR_0 + 1i*Q.gene_GEN_L006_PixelswiseMult_weightI_0);
PWI=reshape(PWI,[100 100 nTS]);

PWIB=squeeze(double(Q.gene_GEN_L006_PixelswiseMult_bias_0(:,:,:,:,1)+1i*Q.gene_GEN_L006_PixelswiseMult_bias_0(:,:,:,:,2)));


MCDR=reshape(PWI,[100 100 nTS]);

Phi1=angle(MCDR);

Phi2=Phi1-Phi1(:,:,1);
Phi2=angle(exp(1i*Phi2));
% fgmontage(Phi2)
%% Run the net
nukData=ADataIsPy(:,:,SliI,13).';
nukData=nukData(:,3:end);
nukDataCC=MultMatTensor(sccmtx(:,1:ncc).',nukData);

DataV=Row(nukDataCC);
In=DataV(NMapCX);
In=permute(In,[2 1 3]);
AfterPWK=squeeze(sum(In.*PWK,3));
AfterPWKB=AfterPWK+PWKB;
% AfterFTX=MultMatTensor(FTX.',AfterPWKB);
AfterFTX=MultMatTensor(FTX.',permute(AfterPWKB,[2 1 3]));
AfterFTY=MultMatTensor(FTY.',permute(AfterFTX,[2 1 3]));
AfterPWI=sum(permute(PWI,[2 1 3]).*AfterFTY,3);
% AfterPWI=sum(PWI.*AfterFTY,3);
AfterPWIB=AfterPWI+PWIB;

ShowAbsAngle(AfterPWIB)

%%
M=load('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/OnRealData.mat');


MC=M.x(:,:,:,1)+1i*M.x(:,:,:,2);
MC=permute(MC,[2 3 1]);

%%
CurRealDataP=[ScanP BaseFN filesep 'RealData' filesep];
mkdir(CurRealDataP);
for r=1:nReps
    nukData=ADataIsPy(:,:,SliI,r).';
    nukData=nukData(:,3:end);
    nukDataCC=MultMatTensor(sccmtx(:,1:ncc).',nukData);

    CurIDataV=Row(nukDataCC.')*RealDataFac;
    CurIDataVR=[real(CurIDataV) imag(CurIDataV)];
        
    Data=repmat(single(CurIDataVR),[16 1]);
    RealDataFN=[CurRealDataP 'Sli' num2str(SliI) '_r' num2str(r,'%02d') '.mat'];
%     RealDataFN=['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RealData/b_Ben14May_Sli5_r' num2str(r,'%02d') '.mat'];
    save(RealDataFN,'Data');
end
%%
AMCWithBias=AMC;
%%
% CurRealDataOutP='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944/RealDataOut/';
CurRealDataOutP='/media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439/RealDataOut/';
% WhichSlices=8:12;
WhichSlices=8;
for SliI=WhichSlices
    for r=1:nReps
        %     M=load(['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/Out/OnRealData' num2str(r,'%02d') '.mat']);
        M=load([CurRealDataOutP 'Sli' num2str(SliI,'%02d') filesep 'OnRealData' num2str(r,'%02d') '.mat']);
        
        MC=M.x(:,:,:,1)+1i*M.x(:,:,:,2);
        MC=permute(MC,[2 3 1]);
        
        AMCS(:,:,SliI,r)=double(MC(:,:,1));
    end
end
%%
TmpP='/media/a/H1/tmp/';
for r=1:nReps
    Raw2Nii(abs(AMCS(:,:,[WhichSlices(1) WhichSlices(1) WhichSlices WhichSlices(end) WhichSlices(end)],r)),[TmpP 'AMCS_' num2str(r,'%02d') '.nii'],'float32');
end

clear rAMCS
for r=1:nReps
    rAMCS(:,:,:,r)=loadniidata([TmpP 'rAMCS_' num2str(r,'%02d') '.nii']);
end
Raw2Nii(abs(rAMCS),'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/rAMCS.nii','float32');

MCE=abs(rAMCS(:,:,:,2:2:end));
MCO=abs(rAMCS(:,:,:,1:2:end));
D=MCE-MCO;
D(isnan(D))=0;
%%
WhichSlicesShow=4:7;
WhichSlicesShow=8;
PerfMap=mean(D(:,:,WhichSlicesShow,3:end-3),4);

PerfMapR=rot90(PerfMap);
% NNResR=rot90(mean(abs(AMCS(:,:,WhichSlices(1:4),:)),4));
NNResR=rot90(mean(abs(rAMCS(:,:,WhichSlicesShow,:)),4));

PerfMapR(NNResR<0.1)=0;
ClimF=[0.0004 4e-3];
A=(min(max(PerfMapR,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
for s=1:numel(WhichSlicesShow)
    RGB(:,:,:,s)=ind2rgb(round(squeeze(A(:,:,s))*255)+1,parula(256));
end
ClimB=[0 7e-1];
Msk=permute(PerfMapR>ClimF(1),[1 2 4 3]);

B=(min(max(NNResR,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+permute(B,[1 2 4 3]).*(1-Msk);
% figure;imshow(X)
X2=CombineDims(X,[4 2]);
figure;imshow(X2)
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
Raw2Nii(abs(AMCS(:,:,WhichSlices,:)),'AMCS.nii','float32');
Y=abs(AMCS(:,:,WhichSlices,:));
X=loadniidata('AMCS.nii');
%%
MCE=abs(mean(AMCS(:,:,WhichSlices,6:2:36),4));
MCO=abs(mean(AMCS(:,:,WhichSlices,5:2:35),4));
D=squeeze(MCE-MCO);
% PerfMap=D;
fgmontage(rot90(D),[0 7e-3])

%%
CurRealDataOutP=[ScanP BaseFN filesep 'RealDataOut' filesep];
% CurRealDataOutP='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944/RealDataOut/';
% /media/a/DATA/180628_AK/meas_MID244_gBP_VD11_U19_G35S155_4min_FID22439
for r=1:nReps
%     M=load(['/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/Out/OnRealData' num2str(r,'%02d') '.mat']);    
    M=load([CurRealDataOutP 'Sli' num2str(SliI,'%02d') filesep 'OnRealData' num2str(r,'%02d') '.mat']);

    MC=M.x(:,:,:,1)+1i*M.x(:,:,:,2);
    MC=permute(MC,[2 3 1]);
    
    AMCS(:,:,SliI,r)=double(MC(:,:,1));
end
disp('ok')
%%
AMC=squeeze(AMCS(:,:,SliI,:));
MCE=abs(AMC(:,:,2:2:end));
MCO=abs(AMC(:,:,1:2:end));
D=MCE-MCO;
fgmontage(mean(D(:,:,10:20),3))
%%
fgmontage(rot90(mean(D(:,:,5:35),3)))
%%
MCE=abs(mean(AMC(:,:,6:2:36),3));
MCO=abs(mean(AMC(:,:,5:2:35),3));
D=MCE-MCO;
PerfMap=D;
fgmontage(rot90(D),[0 7e-3])

%% read For Out
OP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_TSB__2018-07-02_16-46-03_train/';
OFN=[OP 'MatOut_000410.mat'];
OFN=[OP 'MatOut_000470.mat'];
Q=load(OFN);

gene_EndForOut=double(Q.gene_EndForOut(:,:,:,1)+1i*Q.gene_EndForOut(:,:,:,2));
ShowAbsAngle(permute(gene_EndForOut,[2 3 1]))
%% Mat,Out
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_TSB__2018-07-03_09-49-03_train/';
MatFN=[VP 'ForMat_000160.mat'];
OutFN=[VP 'ForOut_000160.mat'];
QM=CombineRIFlds(load(MatFN));
QO=CombineRIFlds(load(OutFN));
Q=QO;
Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    tmp=Q.(Flds{i});
    Q.(Flds{i})=tmp(:,:,:,:,1)+1i*tmp(:,:,:,:,2);
end
for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},67) num2str(size(Q.(SFlds{i})),'%9d         ')]);
end

gene_AfterRegridP_ForOut=squeeze(Q.gene_AfterRegridP_ForOut(5,:,:,:,:));
gene_AfterFT_ForOut=squeeze(Q.gene_AfterFT_ForOut(5,:,:,:,:));
gene_AfterTSC_ForOut=squeeze(Q.gene_AfterTSC_ForOut(5,:,:,:,:));
%%
fgmontage(grmss(permute(Q.gene_AfterFT_ForOut,[2,3,4,1]),3));title('AfterFT');
fgmontage(permute(Q.gene_AfterTSC_ForOut,[2 3 1]));title('AfterTSC');
%%
gene_GEN_L023_TSC_ForMatPixelwiseMultC_weightCX=reshape(QM.gene_GEN_L023_TSC_ForMatPixelwiseMultC_weightC,[128 128 7]);
% AfterFTy=gene_AfterRegridP_ForOut.*permute(QM.gene_GEN_L020_FTy_ForMatadd_Mult2DMCyCSharedOverFeat_weightC,[1 3 4 2]);
X=gene_AfterFT_ForOut.*gene_GEN_L023_TSC_ForMatPixelwiseMultC_weightCX;
Y=sum(X,3);
