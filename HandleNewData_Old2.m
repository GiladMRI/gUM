FN='/home/a/SpiralData/27Apr18/meas_MID80_gSpiral_BOLD_1Sli_1Rep_1Shot_VD10_6ADCs_25_125_FID15656.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID78_gSpiral_BOLD_1Sli_1Rep_1Shot_VD13_6ADCs_25_125_FID15654.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID85_gSpiral_BOLD_1Sli_1Rep_2Shot_VD10_4ADCs_25_125_FID15660.dat';
FN='/home/a/SpiralData/27Apr18/meas_MID87_gSpiral_BOLD_1Sli_4Rep_2Shot_VD10_4ADCs_25_125_FID15662.dat';

FN='/media/a/DATA1/13May18/Phantom/meas_MID364_gBP_VD11_U19_FID17753.dat';
FN='/media/a/DATA1/13May18/Phantom/meas_MID366_gBP_VD11_U10_FID17755.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID407_gBP_VD11_U10_FID17796.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID405_gBP_VD11_U19_FID17794.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID417_gBP_VD15_U20_FID17806.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID426_gBP_VD11_U19_A35S155_FID17815.dat';
FN='/media/a/DATA1/13May18/Me/meas_MID409_gBP_VD11_U19_7ADCs_FID17798.dat';

% FN='/media/a/DATA/14May18/Ben/meas_MID107_gBP_VD11_U19_FID17942.dat';

FN='/media/a/DATA/14May18/Ben/meas_MID109_gBP_VD11_U19_4min_FID17944.dat';

ScanP='/media/a/DATA/14May18/Ben/';
BaseFN='meas_MID111_gBP_VD11_U19_G35S155_FID17946';

FN=[ScanP BaseFN '.dat'];
%%
mkdir([ScanP BaseFN]);
%% Read raw
AData = mapVBVD(FN);
% ADataI=AData.image();
ADataIx=AData.image(:,:,:,:,:,3,:,:,:,:,:,:,:);
% ADataIs=squeeze(ADataI);
% ADataIs=ADataIs(ADataIs(:,:,:,3,:)); % remove weird dimension with zeros
% ADataIs=squeeze(ADataIs);
ADataIsL=squeeze(ADataIx);

for i=1:numel(AData.hdr.Phoenix.sSliceArray.asSlice)
    SLoc(i)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dTra;
end

nSlices=numel(AData.hdr.Phoenix.sSliceArray.asSlice);
% Ord=[1:2:nSlices 2:2:nSlices];
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);


FOVx=AData.hdr.Meas.ReadFOV;
dFOV=FOVx/1000;

paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{20};
spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{6};
paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{8};
VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{5};
paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{10};

save([ScanP BaseFN filesep 'Data.mat'],'ADataIsL','AData');
%%
% save('StatusForTrajDesign1Band_2.mat');
%% given everything but VD and Acc
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);
EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;

BARTTrajMS=BARTTrajx;
BARTTrajxC=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
for i=2:paramLongInterleaves
    CurRotTrajC=BARTTrajxC*exp(1i*2*pi*(i-1)/paramLongInterleaves);
    BARTTrajMS(1,:,i)=real(CurRotTrajC);
    BARTTrajMS(2,:,i)=imag(CurRotTrajC);
end
BARTTrajMS(3,end)=0;
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
figure;plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
Acc=ceil(MaxK*2).^2/nTraj;
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);

gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'Traj'],[]) 
close(gcf);
save([ScanP BaseFN filesep 'Traj.mat'],'BARTTrajMS');
%%

%%
SensX=SensB(:,:,:,6+(1:12));
SensX=gflip(SensX,1:2);

B0Q2=B0Q(:,:,6+(1:12));
B0Q2=gflip(B0Q2,1:2);
%%
SliI=5;
%%
ADataIsP=double(permute(ADataIsL,[1 5 2 3 4]));
nReps=size( ADataIsP,5);
nChannels=size( ADataIsP,3);
%
try
    dx=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
catch
    disp('No x shift!');
    dx=0;
end
% dy=-15;
dy=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
% dy=-45;
read_raw=1;

kx=BARTTrajMS(1,:)*2*pi;
ky=BARTTrajMS(2,:)*2*pi;
modx=double(exp(1i*(dx*kx+dy*ky))');


ADataIsPx=reshape(ADataIsP,[paramLongROSamples,nChannels nSlices nReps]);
ADataIsPy=ADataIsPx(1:size(BARTTrajMS,2),:,:,:,:);

ADataIsPy=RepDotMult(ADataIsPy,modx);

clear ADataIsPx ADataIsP ADataIx
%%
%
% osf = 2; % oversampling: 1.5 1.25
% wg = 3; % kernel width: 5 7
% sw = 8; % parallel sectors' width: 12 16

% NTrg1=[96 96];
% nBands=1;
% nShots=paramLongInterleaves;
% 
% TSC_MB=ones([NTrg1 1 nBands]);

% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),NTrg1,osf,wg,sw,ones(1,size(BARTTrajMS,2),nBands),TSC_MB,SensX(:,:,:,SliI));
%%
nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(SensX(:,:,:,SliI),[1 2 4 3]));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(SensX(:,:,:,SliI),[1 2 4 3]));
% ScrfTV=@(x) gScoreImageSimilarity(RecIfTV(x),ARefImg,RefDx,RefDy);

Lambda=1e-7;
% Rec=RecIfTVs(Lambda);
Rec=RecIfWs(Lambda);

ShowAbsAngle(Rec)

gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'BARTRecon_NoB0_W' num2str(Lambda)],[]) 
close(gcf);
save([ScanP BaseFN filesep 'BARTRecon_NoB0_W' num2str(Lambda) '.mat'],'Rec');

%%
Lambda=1e-9;
for SliI=1:nSlices
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);

    RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(SensX(:,:,:,SliI),[1 2 4 3]));
    RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(SensX(:,:,:,SliI),[1 2 4 3]));

    RecS(:,:,SliI)=RecIfTVs(Lambda);
end
%%
fgmontage(RecS)
gprint(get(gcf,'Number'),[ScanP BaseFN filesep 'BARTRecon_AllS_NoB0_W' num2str(Lambda)],[]) 
close(gcf);
save([ScanP BaseFN filesep 'BARTRecon_AllS_NoB0_W' num2str(Lambda) '.mat'],'RecS');

%% Only one channel
NTrg1=[128 128];
nukData=ADataIsPy(:,20,SliI,1).';
% nukData=nukData(:,42:end);

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(ones(NTrg1),[1 2 4 3]));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(ones(NTrg1),[1 2 4 3]));
% ScrfTV=@(x) gScoreImageSimilarity(RecIfTV(x),ARefImg,RefDx,RefDy);

Lambda=1e-9;
Rec=RecIfTVs(Lambda);
ShowAbsAngle(Rec)

%% Only one channel
for c=1:nChannels
    nukData=ADataIsPy(:,c,SliI,1).';
    % nukData=nukData(:,42:end);
    
    RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(ones(NTrg1),[1 2 4 3]));
    RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajMS(:,1:size(nukData,2)), permute(nukData,[3 2 4 1]), permute(ones(NTrg1),[1 2 4 3]));
    % ScrfTV=@(x) gScoreImageSimilarity(RecIfTV(x),ARefImg,RefDx,RefDy);
    
    Lambda=1e-2;
    RecC(:,:,c)=RecIfTVs(Lambda);
end
%%
ShowAbsAngle(RecC)
ShowAbsAngle(SensX(:,:,:,SliI))
%%
dx=-0
kx=kxFull(B);
ky=kyFull(B);
modx=double(exp(1i*(dx*kx+dy*ky))');
N=470;

A = nuFTOperator([kx ky],[N N]);

%%
kx=BARTTrajMS(1,:).'/EffMaxRes*2*pi;
ky=BARTTrajMS(2,:).'/EffMaxRes*2*pi;

EffMaxRes
N=floor(EffMaxRes);
KK=[kx ky];
A = nuFTOperator(KK,[N N]);

close all

proj=ADataIsPy(:,5,SliI,1).';
% proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;
% proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]);
% figure(5)

QQ = regularizedReconstruction(A,proj',@L2Norm,lambda,'maxit', nit);

%%
figure(1);clf;QQ = regularizedReconstruction(A,proj',@L2Norm,1e+3,'maxit', nit);

%% 2D spiral recon
% close all

temp=squeeze(data(:,c,WhichShots,p,s,1,1,e,r,n,:,1,1,1,1,1));
% proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]).*modx;
proj=reshape(permute(temp,[1 3 2]),[1 NRead*nWhichShots]);
QQ = regularizedReconstruction(A,proj',@L2Norm,lambda,'maxit', nit);

%%
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1/4:(size(kTraj,1)-0.1));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1/4:(size(kTraj,1)-0.1));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;
ADataIsP=double(permute(ADataIs,[1 3 2]));
ADataIsPx=reshape(ADataIsP,[paramLongROSamples,3]);
ADataIsPy=ADataIsPx(1:size(BARTTrajx,2),:);
%%
clear recoAD
Sz3=[88 88];
recoAD(:,:,c) = bart('pics -r:0.1 -R T:7:0:0.0010 -t ',BARTTrajx, ADataIsPy(:,c).', ones(Sz3));
ShowAbsAngle(recoAD(:,:,c))
%%
clear recoAD

Sz3=[88 88];
recoAD(:,:,c) = bart('pics -r:0.1 -R T:7:0:0.0010 -t ',BARTTrajx(:,1:end-2), ADataIsPy(3:end,c).', ones(Sz3));
ShowAbsAngle(recoAD(:,:,c))




%%
nukData=nukData(:,3:end);
nTraj=size(nukData,2);
TimeInMs2=(0:nTraj-1)*2.5/1e3;

Sz2=gsize(SensX,1:2);
%% All B0 effects across time
Mgc=imresizeBySlices(gflip(Mg(:,:,SliI+6),1:2),Sz2);
Mskc=Mgc>7e-5;

B0M2=-B0Q2(:,:,SliI);
B0M2(~Mskc)=0;
%%
% B0M2=B0RealEx;

[U_TimeInMs2, IA_TimeInMs2, IB_TimeInMs2]=unique(TimeInMs2);
nU_TimeInMs2=numel(U_TimeInMs2);

AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
MgcN=Mgc./grmss(Mgc);
WE=Col(MgcN)*0+1;

WeightedE=WE.*E;
%% Fessler time segmentation
% nTS=7;
clear ErrTS
TS_Thresh=1e-5;
for nTS=13:15
    disp(nTS)
    FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    TSC2=reshape(TSC,prod(Sz2),nTS);
    WTSC2=WE.*TSC2;
    tic
%     TSB=(E.')/(TSC2.');% W in both sides
    TSB=(WeightedE.')/(WTSC2.');% W in both sides
    toc
%     ErrTS(nTS)=grmss(E-TSC2*(TSB.')); %include W
    ErrTS(nTS)=grmss(WeightedE-WTSC2*(TSB.')); %include W
    if(ErrTS(nTS)<TS_Thresh)
        disp(['Stopped at #TS=' num2str(nTS) ' err=' num2str(ErrTS(nTS))]);
        break;
    end
end
figure;plot(log10(ErrTS),'-*')
%% GPU TS
BARTTrajAct=BARTTrajMS(:,1:size(nukData,2));

% Sens=imresizeBySlices( squeeze(SensP2),Sz2);
osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

% TSC_MB=ones([NTrg1 1 nBands]);

% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),NTrg1,osf,wg,sw,ones(1,size(BARTTrajMS,2),nBands),TSC_MB,SensX(:,:,:,SliI));

NTrg1=Sz2;
nBands=1;
nShots=paramLongInterleaves;

TSBF=TSB(IB_TimeInMs2,:).';
Sens=squeeze(SensX(:,:,:,SliI));

GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,1:size(nukData,2)),NTrg1,osf,wg,sw,TSBF,TSC,Sens);
GOP_MC=GOP_MCSMBMS;
nTrajP2=nU_TimeInMs2;
%%
% GOP_MC = ggpuNUFT_TS_MCx(BARTTraj2,Sz2,osf,wg,sw,TSB(IB_TimeInMs2,:).',TSC,Sens);

x = randn(Sz2) + 1j*randn(Sz2);
y = randn([size(Sens,3) nTrajP2]) + 1j*randn([size(Sens,3) nTrajP2]);
Ax = GOP_MC*x;
Aty = GOP_MC'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))

%%
% if(nShots==1)
%     TVOP=TVOP_MSlice;
% else
%     TVOP=TVOP_MTC_W([1 1 0 1e1]);
% end
% TVOP=Wavelet;
XFM=1;
% XFM = Wavelet('Daubechies',4,4);	% Wavelet
DataP=nukData;

% AOdd = GOP_MCSMBMS;
AOdd = GOP_MC;

% TVW=0.1;
TVW=3e-6;


% filterType:   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies','Symmlet', 'Vaidyanathan','Battle'
% Use suffix _TI for translation invariance, for example, 'Daubechies_TI'
% filterSize: related to the support and vanishing moments of the particular wavelet (See MakeONFilter in wavelab)
% wavScale: 	scallest scale of wavelet decomposition

XFMStr='Daubechies';

filterSize=4;
wavScale=4;

if(strcmp(XFMStr,'None'))
    XFM=1;
    xfmWeight=0;
else
    XFM = Wavelet(XFMStr,filterSize,wavScale);
    xfmWeight = 3e-6;	% Weight for Transform L1 penalty
end


param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);

param.data =     DataP;

% param.pNorm=1;
% param.TVWeight=1e-5;
% param.TVWeight=0.00001;
% param.TVWeight=0.0001;
% param.TVWeight=0.01*2;

nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
% res=resC;
% res=resB;
if(nShots>1)
    res=repmat(resA,[1 1 1 nShots]);
    res=res+randn(size(res))*max(abs(resA(:)))/20;
end
% max(abs(resA(:)))


FigH=4000;
figure(FigH);close(FigH);

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
param.Verbose=false;
clear ObjConv Score
for n=1:nfnlCgIters
%     disp([Slis WhichRep n]);
    [res, CurObj] = fnlCg(res,param);
%     res=repmat(mean(res,4),[1 1 1 2]);
    ObjConv(n)=CurObj;
%     (ObjConv(end)-ObjConv(end-1))*2/(ObjConv(end)+ObjConv(end-1))
    im_res = param.XFM'*res;
%     figure(FigH), gmontage(abs(gflip(im_res,1))), drawnow;% title(qq)\
%     Score(n)=grmss(CurMBMSIs-im_res);
    if(param.ShowFig)
        figure(FigH); subplot(1,3,1);
        gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,3,2);
        gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
        subplot(1,3,3);
        plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
        
%         subplot(2,3,4);
%         gmontage(abs(CurMBMSIs-im_res),[-100 100]);
        
%         if(nShots>1)
%             subplot(2,3,5);
%             gmontage(RepDotMult(Msks(:,:,1:nBands),abs(cat(4, diff(CurMBMSIs,[],4),diff(im_res,[],4)))),[-100 100]);
%         end
        
%         subplot(2,3,6);
%         plot(Score);setYaxis([0 Score(n)*3]);if(n>1), setXaxis([1 n]);end
    end
%     t=toc;
    if(n>1)
        dObjP=(ObjConv(n-1)-CurObj)/ObjConv(n-1);
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
        if(dObjP<1e-3)
            disp('Not advancing. Stopping.');
            break;
        end
    else
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
    end
    
    if(nShots==1)
        resA=res;
    end
end

fgmontage(im_res,[0 7e-3]);
XLbl=['L' num2str(param.pNorm) ', TVW=' num2str(param.TVWeight) ', ' XFMStr ' Size=' num2str(filterSize) ' Scale=' num2str(wavScale) ' W=' num2str(xfmWeight)];
xlabel(XLbl)
gprint(get(gcf,'Number'),[FN(1:end-4) XLbl],[]) 
save([FN(1:end-4) XLbl '.mat'],'im_res');










%% Fit B0:
Sz1=Sz2;
Area=[5 5];
Area=[7 7];

I=Mgc;
B0Real=B0M2;

TmpMaskA=imfillholesBySlices(I>5e-5);
TmpMaskE=TmpMaskA;
% se = strel('disk',1);
% TmpMaskE=imerode(TmpMaskA,se);

Z=zeros(Sz1);
MidSz=floor(Sz1/2);
Reg=false(Sz1);
Reg((MidSz(1)-Area(1)):(MidSz(1)+Area(1)),(MidSz(2)-Area(2)):(MidSz(2)+Area(2)))=true;
F=find(Reg);
n=numel(F);
ParamToField=@(x) abs(fft2c(PutValsInIdxs(Z,x(1:n)+1i*x(n+(1:n)),F)));
GetValsInLoc=@(A,B) A(B);
FMsk=find(TmpMaskE);


ToAdd=-min(B0Real(:))+20;


% CostFunc=@(x) grmss(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+200,FMsk));

W=GetValsInLoc(I,FMsk).';
CostFunc=@(x) W*abs(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+ToAdd,FMsk));

options = optimset('Display','iter');
B0M=(B0Real+ToAdd).*TmpMaskE;
B0M(~isfinite(B0M))=0;
C=ifft2c(B0M);
X0=[real(C(Reg)).' imag(C(Reg)).'];
[BestX, BestCost]=fminsearch(CostFunc,X0,options);

% fgmontage(ParamToField(BestX)-ToAdd,[-200 300]);colorbar
% fgmontage((B0Real).*TmpMaskE,[-200 300]);colorbar
% 
% fgmontage(ParamToField(BestX).*TmpMaskA,[0 500]);colorbar
% 
% fgmontage(ParamToField(BestX),[0 500]);colorbar

B0RealEx=ParamToField(BestX)-ToAdd;
%%
fgmontage(B0Real,[-400 300])
fgmontage(B0RealEx,[-400 300])

fgmontage(B0Real-B0RealEx,[-100 100])
%%
fo = fitoptions('Weights',z*0+1);
    
load franke
sf = fit([x, y],z,'poly23',fo)
figure;plot(sf,[x,y],z)
%%
ToFit=B0Q2(:,:,SliI);
W=abs(Mgc);
[x,y]=meshgrid(1:128,1:128);
xF=x(:);
yF=y(:);
wF=W(:);
zF=ToFit(:);

DForFit=3;
x=xF(1:DForFit:end);
y=yF(1:DForFit:end);
z=zF(1:DForFit:end);
w=wF(1:DForFit:end);
fo = fitoptions('Weights',w,'Method','LowessFit','Span',0.01);
tic
sf = fit([x, y],z,'lowess',fo)
% sf = fit([x, y],z,fo)
toc
% figure;plot(sf,[x,y],z)
tic
X=sf([xF(1:1:end),yF(1:1:end)]);
toc
X2=reshape(X,Sz2);

fgmontage(ToFit,[-800 400])
fgmontage(X2,[-800 400])