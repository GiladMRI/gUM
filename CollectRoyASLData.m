BaseRoyP='/media/a/DATA1/FromRoy/';
D=dir(BaseRoyP);
D=D([D.isdir]);
D=D(3:end);
%% load all
for i=1:numel(D)
    CurP=[BaseRoyP D(i).name '/niftis/asl/'];
    D2=dir([CurP '*.nii.gz']);
    D2N={D2.name};
    D2N=D2N{strhas(D2N,'_RS_')};
    RoyAll(:,:,:,:,i)=int32(loadniidata([CurP D2N]));
end
% save([BaseRoyP 'All.mat'],'RoyAll');
%% collect 2 timepoints, 2 images sets
nTP=size(RoyAll,4)/2;
nS=size(RoyAll,3)/2;
MaxTimeDiff=5;
nP=size(RoyAll,5);
CouplesAll=int32(zeros([nP*nS*2*(nTP-MaxTimeDiff),gsize(RoyAll,1:2)  ,2,2]));
c=0;
for i=1:nP
    for s=1:nS
        for k=1:2 % tagged and untagged
            for t=1:(nTP-MaxTimeDiff)
                tpd=randi(MaxTimeDiff);
                CurCouple=squeeze(RoyAll(:,:,[s s+nS],t*2-2+k+[0 tpd],i));
                c=c+1;
                CouplesAll(c,:,:,:,:)=CurCouple;
                disp(c);
            end
        end
    end
end
save([BaseRoyP 'AllT2S2.mat'],'CouplesAll');
%%
Ni=size(CouplesAll,1);
ToFlip=rand(Ni,1)>0.5;
CouplesAll(ToFlip,:,:,:,:)=flip(CouplesAll(ToFlip,:,:,:,:),1);
ToFlip=rand(Ni,1)>0.5;
CouplesAll(ToFlip,:,:,:,:)=flip(CouplesAll(ToFlip,:,:,:,:),2);
ToFlip=rand(Ni,1)>0.5;
CouplesAll(ToFlip,:,:,:,:)=permute(CouplesAll(ToFlip,:,:,:,:),[1 3 2 4 5]);
ToShift=randi(11,[Ni 1])-6;
for i=1:Ni
    CouplesAll(Ni,:,:,:,:)=circshift(CouplesAll(Ni,:,:,:,:),ToShift(i),1);
end
ToShift=randi(11,[Ni 1])-6;
for i=1:Ni
    CouplesAll(Ni,:,:,:,:)=circshift(CouplesAll(Ni,:,:,:,:),ToShift(i),2);
end
%%
% save([BaseRoyP 'AllT2S2.mat'],'CouplesAll');
load([BaseRoyP 'AllT2S2.mat'],'CouplesAll');
% load('CurStatusASL_2x2_Try2.mat');
% getSmallVars
% save('CurStatusASL_2x2_Try2.mat',SmallVars{:});
%% Masks
% load([BaseRoyP 'SensBoth.mat'],'SensBoth');

Msks=grmss(SensBoth,3)*sqrt(ncc)>0.3;
fgmontage(Msks);
Msksx=repmat(Msks,[1 1 1 2]);
%% Design trajectory
SzTK=[70 70];

% SzTK=[196 196];
% SzTK=[128 128];
FOVx=192;
dFOV=FOVx/1000;
nADCs=1;
paramLongROSamples=1024*nADCs;
spBW =4e5;
paramLongInterleaves=1;
AccR = (prod(SzTK)*pi/4)   /(paramLongInterleaves*paramLongROSamples);
VD =1.0;
paramLongSpGradAmp =25;
paramLongSpSlewRate =125;
[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate]);

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



% nShots=1;
nShots=paramLongInterleaves;
if(nShots==1)
%     BARTTrajMS=cat(2,BARTTrajMS(:,:,1),BARTTrajMS(:,:,2));
end
figure;
for i=1:nShots
    plot(BARTTrajMS(1,:,i),BARTTrajMS(2,:,i));
    hold on
end

nTraj=size(BARTTrajMS,2);
% title([BaseRes size(kTrajQ,1)])
title([BaseRes nTraj AccR])
%
%% CAIPI form
CAIPIBlockSize=10;
CAIPITraj=repmat([zeros(1, CAIPIBlockSize),ones(1, CAIPIBlockSize)],[1 ceil(nTraj/CAIPIBlockSize)]);
CAIPITraj=CAIPITraj(1,1:nTraj);
CAIPIFacSli2T1=1-CAIPITraj*2;
CAIPIFacSli2T2=-CAIPIFacSli2T1;
%
nBands=2;

BARTTrajMS(3,:,1)=CAIPITraj;
BARTTrajMS(3,:,2)=1-CAIPITraj;

osf = 2; % oversampling: 1.5 1.25
wg = 3; % kernel width: 5 7
sw = 8; % parallel sectors' width: 12 16

TSC_MB=ones([SzTK 1 nBands]);
TSB_MB=ones(1,nTraj,nBands);
GOP_MCSMBMS = ggpuNUFT_TS_MC_MB_MS(BARTTrajMS(:,:,1:nShots),SzTK,osf,wg,sw,TSB_MB,TSC_MB,SensBoth(:,:,:,1:nBands));
%%
CurMBMSIs=permute(double(CouplesAll(349,:,:,1:nBands,1:nShots)),[2 3 4 5 1]);
CurMBMSIs=RepDotMult(CurMBMSIs,Msks(:,:,1:nBands));
% CurMBMSIs=repmat(CurMBMSIs(:,:,:,1),[1 1 1 2]);

CurData=GOP_MCSMBMS*CurMBMSIs;
IT=GOP_MCSMBMS'*CurData;

SzIMBMS=[SzTK nBands nShots];
SzDataMBMS=[ncc nTraj nShots];
x = randn(SzIMBMS) + 1j*randn(SzIMBMS);
y = randn(SzDataMBMS) + 1j*randn(SzDataMBMS);
Ax = GOP_MCSMBMS*x;
Aty = GOP_MCSMBMS'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
%% Miki's nufft
tmpI=repmat(double(squeeze(CouplesAll(349,:,:,1,1))),[1 1 15]);
tmpI(1)=tmpI(1)+1i*0.01;
tmpK=BARTTrajMS(:,:,1);

writecfl('/home/a/gpu-nufft-master/bin/a', tmpI*10);
writecfl('/home/a/gpu-nufft-master/bin/k', tmpK);

~/gpu-nufft-master/bin$ ./nufft_main k a b

B=readcfl('/home/a/gpu-nufft-master/bin/b');

%% Recon
% TVOP=TVOP_MB_MT_W([1 1 0 1e7]);
% TVOP=TVOP_MT_W([1 1 0 1e7]);
if(nShots==1)
    TVOP=TVOP_MSlice;
else
    TVOP=TVOP_MTC_W([1 1 0 1e1]);
end
DataP=CurData;

AOdd = GOP_MCSMBMS;

TVW=0.1;

param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',1,'TV',TVOP,'xfmWeight',0),init);

param.data =     DataP;

param.pNorm=1;
% param.TVWeight=0.000001;
% param.TVWeight=0.0001;
param.TVWeight=0.01*2;

nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(SzIMBMS);
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
    Score(n)=grmss(CurMBMSIs-im_res);
    if(param.ShowFig)
        figure(FigH); subplot(2,3,1);
        gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(2,3,2);
        gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
        subplot(2,3,3);
        plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
        
        subplot(2,3,4);
        gmontage(abs(CurMBMSIs-im_res),[-100 100]);
        
        if(nShots>1)
            subplot(2,3,5);
            gmontage(RepDotMult(Msks(:,:,1:nBands),abs(cat(4, diff(CurMBMSIs,[],4),diff(im_res,[],4)))),[-100 100]);
        end
        
        subplot(2,3,6);
        plot(Score);setYaxis([0 Score(n)*3]);if(n>1), setXaxis([1 n]);end
    end
%     t=toc;
    disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
    
    if(nShots==1)
        resA=res;
    end
end
%%
% %% Mask (1,2) and Add random phase to data, different phase between slices, same phase for both timepoints
% % CouplesAll [Idx X Y sli time]
% save('CurStatusASL_2x2_Try1.mat');
% for all examples:
% CurIdx=7800;
% P1=GenerateRandomSinPhase(SzTK,5,.1);
% P2=GenerateRandomSinPhase(SzTK,5,.1);
% CurIs=double(squeeze(CouplesAll(CurIdx,:,:,:,:)));
% CurIs(:,:,1,:)=RepDotMult(CurIs(:,:,1,:),P1);
% CurIs(:,:,2,:)=RepDotMult(CurIs(:,:,2,:),P2);
% CurIs=CurIs.*Msksx;
% ShowAbsAngle(CurIs)
% 
% DataSli1T1=GOP_MCS{1,1}*CurIs(:,:,1,1);
% DataSli2T1=GOP_MCS{2,1}*CurIs(:,:,2,1);
% 
% DataSli1T2=GOP_MCS{1,2}*CurIs(:,:,1,2);
% DataSli2T2=GOP_MCS{2,2}*CurIs(:,:,2,2);
% 
% DataSli2T1C=DataSli2T1.*CAIPIFacSli2T1S;
% DataSli2T2C=DataSli2T1.*CAIPIFacSli2T2S;
% 
% DataT1=DataSli1T1+DataSli2T1C;
% DataT2=DataSli1T2+DataSli2T2C;
% %%
% ggpuNUFT_TS_MCx
% 
% %%
% Now send to TF and try to predict mag only. Cost on both slices at both time points, and on diff between time points
% %%
% 
% ToFlip=rand(Ni,1)>0.5;
% HCPData(:,:,ToFlip)=flip(HCPData(:,:,ToFlip),2);
% ToFlip=rand(Ni,1)>0.5;
% HCPData(:,:,ToFlip)=permute(HCPData(:,:,ToFlip),[2 1 3]);
% 
% 
% A=loadniidata('/media/a/DATA1/FromRoy/S021/niftis/asl/RoyHaa_270715_C001S021_20150727_001_033_BP_FAIR_2x2_z3_RS_BP_FAIR_2x2_z3_RS_E00_M_Echo_0.nii.gz');
% %%
% AE=mean(A(:,:,5:end,2:2:end),4);
% AO=mean(A(:,:,5:end,1:2:end),4);
% D2=AE-AO;
% fgmontage(D2(:,:,2:5:end));colorbar;title('2')
% 
% fgmontage(AE(:,:,2:5:end));colorbar;title('AE')
% fgmontage(AO(:,:,2:5:end));colorbar;title('AO')
% %%
% AE3=mean(A(:,:,5:end,1:80),4);
% AO3=mean(A(:,:,5:end,81:end),4);
% 
% fgmontage(AE3-AO3);colorbar;title('3')
% %%
% DSli=12;
% LowSli=5;
% TimeDiff=3;
% FirstIdx=35;
% A2=A(:,:,LowSli+[0 DSli],FirstIdx+[0 TimeDiff*2]);
% %%
% circshift
% Pretrain on everything
% then train on this data
% cost on difference
% flip both and rot
% 
% 5*12*2*70*5*8
% 
% 20000/5/12
% 
% 5*12*75
% 
% Choose sens maps and B0 map
% %%
% Then resahpe or otherwise enlarge to e.g. 100x100
%%
128*128*2*5*8*8/1e6
2*