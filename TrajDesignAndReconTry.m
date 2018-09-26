%% Later take consideration of noise
%% Get realistic sensitivity maps
Sens=load('/home/deni/GL/meas_MID399_BP_fieldmap_4echos_FID17788/Sens.mat');
Sens=Sens.SensB;
SensX=Sens(:,:,:,12,1);
Msk=grmss(SensX,3)>0.1;
%% get some example image
load brain_8ch
I=ifft2c(DATA);
SensTmp=RunESPIRiTForSensMaps(I);
Im=CalcSENSE1f(I,SensTmp);
Im=rot90(Im);
Im=Im/grmss(Im);
%%
Im=double(rgb2gray(imread('~/Downloads/maxresdefault.jpg')));
Im=rot90(Im);
Im=Im/grmss(Im);

%%
%OutInOut=false;
% Fixed parameters
NumberOfADCBlockss=5; % each block is 1024 points, 2.5ms
FOVToUse=192; % in mm
AcqBW=400000;
MaxGradToUse=35;
MaxSlewToUse=155;
InterleavesToUse=1;
% Parameters to play with
AccToTryNow=2.4;
VDToTryNow=1.3;


paramLongROSamples=NumberOfADCBlockss*1024;
dFOV=FOVToUse/1000;
spBW=AcqBW;
AccR=AccToTryNow;
paramLongInterleaves=InterleavesToUse;
VD=VDToTryNow;
FOVx=FOVToUse;

OutInOut=paramLongROSamples>10000;



if(OutInOut)
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/3,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);
else
    [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([FOVToUse/1000,NumberOfADCBlockss*1024,AcqBW,AccToTryNow,...
        1,VDToTryNow,MaxGradToUse,MaxSlewToUse,0]);
end
WhichE=1;
if(OutInOut)
    BaseLen=size(kTraj,1)/3;
    WhichE=2;
    TrajIdxs=(1:BaseLen)+BaseLen*(WhichE-1);
    kTraj=kTraj(TrajIdxs,:);
    GradBuf=GradBuf(TrajIdxs,:);
end


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
disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
ylabel(['Nominal acc: ' num2str(AccR) ' VD=' num2str(VDToTryNow)]);
ResInmm=FOVToUse/(MaxK*2);
xlabel([num2str(ResInmm) 'mm']);
subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
%% So far we got the trajectory
%%
ResWeAimAt=1.5;
ResBase=192/ResWeAimAt;

ResToUse=ceil(ResBase/8)*8; % Divisable by 8
%ResToUse=2^ceil(log2(ResBase)); % "dydic"
Sz2=[ResToUse ResToUse];

MskR=imresize(Msk,Sz2)>0.1;
SensR=imresizeBySlices(SensX,Sz2);
ImRsized=imresizeBySlices(Im,Sz2).*MskR;
%% Create the data
% Sens=imresizeBySlices( squeeze(SensP2),Sz2);
osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

% TSC_MB=ones([NTrg1 1 nBands]);

% GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),NTrg1,osf,wg,sw,ones(1,size(BARTTrajMS,2),nBands),TSC_MB,SensX(:,:,:,SliI));

NTrg1=Sz2;
nBands=1;
nShots=1;

TSBFX=ones(1,nTraj);

TSCX=ones(Sz2);

GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS,NTrg1,osf,wg,sw,TSBFX,TSCX,SensR);
GOP_MC=GOP_MCSMBMS;
% nTrajP2=nU_TimeInMs2;
disp('ok');
%% Get the simulated signal
SimData=GOP_MC*ImRsized;
%% Reconstruct the image from the signal
nukDataP=permute(SimData,[3 2 4 1]);

SensP=permute(SensR,[1 2 5 3 4]);

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajMS, nukDataP, SensP(:,:,:,:,1));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajMS, nukDataP, SensP(:,:,:,:,1));


Lambda=1e-5;

RegStr='TV';
Rec=RecIfTVs(Lambda);
% RegStr='Wavelet';
% Rec=RecIfWs(Lambda);
RecN=Rec*grmss(ImRsized)/grmss(Rec);

[ssimval, ssimmap] = ssim(double(abs(RecN)),double(abs(ImRsized)));

fgmontage(cat(3,ImRsized,RecN,(ImRsized-RecN)*10,ssimmap*0.7*max(abs(RecN(:)))),[0 max(abs(RecN(:)))*0.7]);
xlabel(['BART Lambda= ' num2str(Lambda)]);

% fgmontage(RecMM)
%
title(['BART ' RegStr ' Lambda=' num2str(Lambda) ', VD=' num2str(VDToTryNow) ' Acc=' num2str(Acc) ' Target res ' num2str(ResInmm) 'mm, SSIM=' num2str(ssimval,'%.2f')]);
%%

   
DataP=SimData;

% QQ=load('/media/a/H1/TFTSNUFTOut.mat');
% WithTSBR=double(QQ.Out)/100;
% DataP=WithTSBR.';

AOdd = GOP_MCSMBMS;
% AOdd = GOP_MC;

% TVW=0.1;
TVW=1e-5;


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
    xfmWeight = 1e-5;	% Weight for Transform L1 penalty
end


param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);

% XFM=SubSpaceOp(V,SdT,1);
% xfmWeight=1e-1;
% % param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);
% % param=ExtendStruct(struct('pNorm',2,'TVWeight',0,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',1,'xfmWeight',xfmWeight),init);
% param=ExtendStruct(struct('pNorm',2,'TVWeight',0.01,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);


param.data =     DataP;

nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
if(nShots>1)
    res=repmat(resA,[1 1 1 nShots]);
    res=res+randn(size(res))*max(abs(resA(:)))/20;
end

FigH=4000;
figure(FigH);close(FigH);

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
param.Verbose=false;
% res=double(MLN*grmss(im_resA)/grmss(MLN));
% res=im_resA;
% res=XFM*double(MLN*grmss(im_resA)/grmss(MLN));
% res=XFM*im_resA;
clear ObjConv Score
for n=1:nfnlCgIters
% for n=1:2
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
        if(dObjP<2e-3)
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
close(FigH);
disp('ok im_res');

% fgmontage(im_res);


RecNSparseMRI=im_res*grmss(ImRsized)/grmss(im_res);

[ssimvalSparseMRI, ssimmapSparseMRI] = ssim(double(abs(RecNSparseMRI)),double(abs(ImRsized)));

fgmontage(cat(3,ImRsized,RecNSparseMRI,(ImRsized-RecNSparseMRI)*10,ssimmapSparseMRI*0.7*max(abs(RecN(:)))),[0 max(abs(RecNSparseMRI(:)))*0.7]);
xlabel(['BART Lambda= ' num2str(Lambda)]);

% fgmontage(RecMM)
%
title(['SparseMRI ' XFMStr ' Lambda=' num2str(TVW) ', VD=' num2str(VDToTryNow) ' Acc=' num2str(Acc) ' Target res ' num2str(ResInmm) 'mm, SSIM=' num2str(ssimval,'%.2f')]);
