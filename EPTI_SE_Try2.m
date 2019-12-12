%% Recon basic
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
BaseSEP='/autofs/cluster/kawin/Gilad/EPTI_fmri_SE/';
SliI=6;
RunI=1;
DynI=7;

OutP=[BaseSEP 'splitProx/'];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);

ES_ms=1.08;
FirstTE_ms=9;
%%
for SliI=6
for RunI=1:5
    for DynI=1:32
        disp([SliI RunI DynI]);
        clear RawData GRecon
RawData{SliI,RunI,DynI}=load([BaseSEP 'Undersampled/SE_Run' num2str(RunI) '_Dyn_' num2str(DynI) '_Slice' num2str(SliI) '.mat']);
GRecon{SliI,RunI,DynI}=load([BaseSEP 'Recon_B0GRAPPA/Recon_RawCalib_SE_Run' num2str(RunI) '_Dyn' num2str(DynI) '_Slice' num2str(SliI) '_SE_Noupdates.mat']);

MM=ifftshift(ifft2c(GRecon{SliI,RunI,DynI}.recon),2);
nSEEchos=size(MM,3);
SESz=gsize(MM,1:2);
if(~exist('CurSens','var'))
    CurSens=RunESPIRiTForSensMapsMultiMap(squeeze(MM(:,:,28,:)),0);
    CurSens=CurSens(:,:,:,1);
    SensMsk=grmss(CurSens,3)>0.01;
    CurSens=perm43(CurSens);
    CurSensS(:,:,SliI,:)=CurSens;
end
MMC=sum(MM.*conj(CurSens),4);
MMCS(:,:,:,SliI)=MMC;
% end
WhichEchosToUseA=11:28;
[PDBaseA, UpdatedB0Map_HzA, UpdatedT2SMap_msA, s_valsA, Fitted0A, PDBase0A]=...
    FitToModel_MPBD1CSf(MMC(:,:,1:28),WhichEchosToUseA,ES_ms,FirstTE_ms);
WhichEchosToUseB=1:18;
[PDBaseB, UpdatedB0Map_HzB, UpdatedT2SMap_msB, s_valsB, Fitted0B, PDBase0B]=...
    FitToModel_MPBD1CSf(MMC(:,:,29:end),WhichEchosToUseB,ES_ms,FirstTE_ms+28*ES_ms);

WhichEchosToUseC=11:28;
[PDBaseC, UpdatedB0Map_HzC, UpdatedT2SMap_msC, s_valsC, Fitted0C, PDBase0C]=...
    FitToModel_MPBD1CSf(flip(MMC(:,:,29:end),3),WhichEchosToUseC,ES_ms,FirstTE_ms);

WhichEchosToUseD=1:18;
[PDBaseD, UpdatedB0Map_HzD, UpdatedT2SMap_msD, s_valsD, Fitted0D, PDBase0D]=...
    FitToModel_MPBD1CSf(flip(MMC(:,:,1:28),3),WhichEchosToUseD,ES_ms,FirstTE_ms);

s_valsBoth=s_valsA+s_valsB;
Msk=s_valsBoth(:,:,1)>3.5;

UpdatedB0Map_HzM=(UpdatedB0Map_HzA+UpdatedB0Map_HzB)/2;

PDBase0Ax=min(abs(PDBase0A),getPercentile(abs(PDBase0A(:)),0.9)*2).*exp(1i*angle(PDBase0A));
PDBase0Bx=min(abs(PDBase0B),getPercentile(abs(PDBase0B(:)),0.9)*2).*exp(1i*angle(PDBase0B));
%%
RA=1./UpdatedT2SMap_msA; % R2-R2*
RB=1./UpdatedT2SMap_msB; % R2+R2'
RC=1./UpdatedT2SMap_msC; % 
RD=1./UpdatedT2SMap_msD; % R2'-R2
% RM=(RA+RB)/2;
% RD=(RB-RA)/2;
% T2p=1./RM;
% T2s=1./RD;

RM=(RB+RD)/2;
Rx=(RB-RD)/2;
T2p=1./RM;
T2=1./Rx;
T2s=UpdatedT2SMap_msB;
%%
% TimePoints=FirstTE_ms+(0:(nSEEchos-1))*ES_ms;
TimePoints=((1:nSEEchos)-29)*ES_ms;
TimePoints_ms3=perm32(TimePoints);

TimePointsSE=abs((1:nSEEchos)-29)*ES_ms;
TimePointsSE_ms3=perm32(TimePointsSE);
TSCM_T2=exp(-TimePoints_ms3./T2);
TSCM_T2s=exp(-TimePointsSE_ms3./T2p);
TSCP=exp(-1i*2*pi*UpdatedB0Map_HzM.*TimePoints_ms3/1000);
%%
ElemsL{1}=perm71(ones(nSEEchos,1));
ElemsL{2}=1i*perm71(ones(nSEEchos,1));
ElemsL{3}=perm73(-1i*2*pi*TimePoints_ms3/1000);
ElemsL{4}=perm73(-TimePoints_ms3);
ElemsL{5}=perm73(-TimePointsSE_ms3);
%%
% C0=(PDBase0Ax+PDBase0Bx)/2;
C0=(Fitted0A(:,:,end)+sum(Fitted0B(:,:,1:2),3))/3;
Est0=C0.*TSCP.*TSCM_T2.*TSCM_T2s;
%%
kdata_SE=permute(RawData{SliI,RunI,DynI}.kdata_SE,[2 3 1 4]);
% recon=GRecon{SliI,RunI,DynI}.recon;

mask_sample=perm73(abs(kdata_SE(:,:,:,1))>0);
sig=perm73(kdata_SE);
sig=sig.*( (-1).^(0:(SESz(2)-1)) );
%%
SensCSMap=CurSens;

MskRec=mask_sample;
% sig=single(MskCC);
% sig(MskCC)=SigMskd(:,SliI);
% sig=perm75(sig);

m0=abs(C0);
p0=angle(C0);
b0=UpdatedB0Map_HzM;
t0a=max(5,min(200,abs(T2)))*0+70;
t0b=max(5,min(2000,abs(T2p)))*0+800;

Elem0={m0,p0,b0,t0a,t0b};
%%
Mm0=ElemsL{1}.*m0;
Pp=exp(ElemsL{2}.*p0);
Bb=exp(ElemsL{3}.*b0);
Tta=exp(ElemsL{4}.*1./t0a);
Ttb=exp(ElemsL{5}.*1./t0b);
Est0x=Mm0.*Pp.*Bb.*Tta.*Ttb;
%%
SRange=[21 21];

SSigs=[0.0001 0.1:0.1:29];
m0b=m0; 
M0B0=m0b.*b0;

% SM0B0=SmoothBySlices(M0B0,SRange,SSig);
% SM0=SmoothBySlices(m0b,SRange,SSig);
% SB0=SM0B0./SM0;
clear SM0B0 SM0
for i=1:numel(SSigs)
    SM0B0(:,:,i)=SmoothBySlices(M0B0,SRange,SSigs(i));
    SM0(:,:,i)=SmoothBySlices(m0b,SRange,SSigs(i));
end
mThresh=4;
MskSM0=(SM0.*perm32(SSigs))>mThresh;
MskSM0(:,:,end)=1;
FirstTimeOverT=numel(SSigs)+1-sum(MskSM0,3);
SB0=SM0B0./SM0;
clear SB0x
for i=1:size(SB0,1)
    for j=1:size(SB0,2)
        SB0x(i,j)=SB0(i,j,FirstTimeOverT(i,j));
    end
end
SB0x(isnan(SB0x))=0;
% MaxB0=400;
% fgmontage(SB0x,[-400 400])
% fgmontage(b0a,[-400 400])
% fgmontage(SB0,[-400 400])
%% 




% CurSPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S4_'];
CurSPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S' num2str(SliI) '_'];

CurSDRPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S' num2str(SliI) '_' 'Run' num2str(RunI) '_Dyn' num2str(DynI) '_'];

CurPrefix='zz_';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

nccToUse=32;
disp('Prepare folders');
system('rm /tmp/*.cfl');
system('rm /tmp/*.hdr');
%%
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
% ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 1','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

ImSz16=FillOnesTo16([SESz 1 1 1 1 nSEEchos]);
BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap(:,:,:,1:nccToUse) mask_sample};

SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
% SigToUse=ifft1cg(SigToUse,1);

BARTS_Aop=BARTS_Aopx;

writecfl([OutP 'mask_sample'],mask_sample);
BARTS_Aop.Others{2}=[OutP 'mask_sample'];

writecfl([CurSPrefix 'SensCC'],SensCSMap(:,:,:,1:nccToUse));
BARTS_Aop.Others{1}=[CurSPrefix 'SensCC'];

SigFN=[CurSDRPrefix 'Sig'];
writecfl(SigFN,SigToUse);

ksp_adjFN=[CurSDRPrefix 'sig_adj'];

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
writecfl(ksp_adjFN,ksp_adj);

disp('got ksp_adj');
UseOS=false;
UseOSexplicit=false;
%
% ThroughPlaneDecayFN=[CurSPrefix 'ThroughPlaneDecay'];
% writecfl(ThroughPlaneDecayFN,perm73(ThroughPlaneDecayS(:,:,:,SliI)));
%%
% Elem0{1}=Elem0{1}*0;
Elem0{1}=m0;
Elem0{2}=p0;
%
% Elem0{2}([1:2 end-1:end],:)=0;
% Elem0{2}(:,[1:2 end-1:end])=0;
Elem0{3}=SB0x;
% Elem0{3}=b0a;
Elem0{4}=max(5,min(200,abs(Elem0{4})))*0+70;
Elem0{5}=max(5,min(2000,abs(Elem0{5})))*0+800;
disp('Improved initials');
%%
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+3,1e+4];
ElemsLambda=[1e-7,0.1,1e-3,1e-4,1e-4];
% Try 2
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+3,1e+4];
ElemsLambda=[1e-5,0.1,1e-3,1e-4,1e-4];
% Try 3
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+3,1e+3];
ElemsLambda=[1e-5,0.1,1e-3,1e-4,1e-4];
% Try 4
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
ElemsLambda=[1e-5,0.1,1e-3,1e-4,1e-4];
% Try 5
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
ElemsLambda=[1e-4,0.1,1e-3,1e-4,1e-4];
% Try 6
ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
ElemsLambda=[1e-2,0.1,1e-3,1e-4,1e-4];
% % Try 7 - too smooth
% ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
% ElemsLambda=[1e-0,0.1,1e-3,1e-4,1e-4];
% % Try 8
% ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
% ElemsLambda=[1e-1,0.1,1e-3,1e-4,1e-4];
% % Try 9
% ElemsAlphas=[1e-4, 1e-2, 1e+0, 1e+2,1e+4];
% ElemsLambda=[10^(-1.5),0.1,1e-3,1e-4,1e-4];

for i=1:numel(Elem0)
    writecfl([CurSDRPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
    writecfl([CurSDRPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
end

writecfl([CurSDRPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSDRPrefix 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4 4];
writecfl([CurSDRPrefix 'ElementTypes'],ElementTypes.');

nElements=5;
ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
%     ninneriterBART(i)=1;
end
writecfl([CurSDRPrefix 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([CurSDRPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSDRPrefix 'Elem' num2str(i-1) '.cfl']);
end
ContinueRun=false;
if(ContinueRun)
    for i=1:numel(Elem0)
        writecfl([CurSDRPrefix 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end
%%
QQQQ=bartCmd(['splitProx -i 5000 -s 600 -I 1000 -d 2 -g -f -F ' CurSDRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
system(QQQQ)
    end % Dyns
end % Runs
end % Slices 
%%

%%
ErrVec=readcfl([CurSPrefix 'ErrVec']);
ErrVec(isnan(ErrVec))=-1;
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
figure;plot(ErrVec)
%%
ErrVec=readcfl([CurSPrefix 'ErrVec']);
ErrVec(isnan(ErrVec))=-1;
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);

for i=1:5
    Maps{i}=readcfl([CurSPrefix 'Elem' num2str(i-1)]);
end

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtMa = exp(ElemsL{4} .* (1./Maps{4}));
expTtMb = exp(ElemsL{5} .* (1./Maps{5}));
RecM=MmM.*expPpM.*expBbM.*expTtMa.*expTtMb;
RecMX=squeeze(sum(RecM,CS_Dim));

disp('Loaded maps');
%
OSFac=1;
ElemRanges={[0 1e-0],[-pi pi],[-400 400],[0 100],[0 2000]};
mThresh=3e-8;
figure;
for i=1:5
    subplot(2,3,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); xlabel([num2str(numel(ErrVec)/sum(ninneriterBART)) ' : ' num2str(ErrVec(end),'%.7g')]); end
    if(i==3), xlabel(num2str(ElemsAlphas,' %.9g,')); end %ylabel(Pref,'Interpreter','None'); end
    if(i==4), xlabel(num2str(ElemsLambda,' %.9g,')); end
end
%%
for SliI=4:9

CurSPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S' num2str(SliI) '_'];

for i=1:5
    Maps{i}=readcfl([CurSPrefix 'Elem' num2str(i-1)]);
end
MapsS{:,SliI}=Maps;
MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtMa = exp(ElemsL{4} .* (1./Maps{4}));
expTtMb = exp(ElemsL{5} .* (1./Maps{5}));
RecM=MmM.*expPpM.*expBbM.*expTtMa.*expTtMb;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMXS(:,:,:,SliI)=RecMX;
end
%%
Both=cat(5,MMCS,RecMXS);
XX=perm43(squeeze(rot90(Both(:,:,29,4:9,:),3)));
YY=CombineDims(PartitionDim(XX,4,2),[5 3]);
fgmontagex(YY(20:end-20,:,:,:))
%%
Maps=zeros([168,96,5,32,5]);

RecMXSDR=zeros(168,96,56,32,5);
for SliI=6
    for RunI=1:5
        for DynI=1:32
            disp([SliI RunI DynI]);
            CurSPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S' num2str(SliI) '_'];
            CurSDRPrefix=['/autofs/cluster/kawin/FuyixueWang/Share/DataEPTI_20190719_SEfMRI_invivo/raw_and_recon/Try6S' num2str(SliI) '_' 'Run' num2str(RunI) '_Dyn' num2str(DynI) '_'];
            for i=1:5
                Maps(:,:,i,DynI,RunI)=readcfl([CurSDRPrefix 'Elem' num2str(i-1)]);
            end
            
            MmM=ElemsL{1}.*Maps(:,:,1,DynI,RunI);
            expPpM = exp(ElemsL{2} .* Maps(:,:,2,DynI,RunI));
            expBbM = exp(ElemsL{3} .* Maps(:,:,3,DynI,RunI));
            expTtMa = exp(ElemsL{4} .* (1./Maps(:,:,4,DynI,RunI)));
            expTtMb = exp(ElemsL{5} .* (1./Maps(:,:,5,DynI,RunI)));
            RecM=MmM.*expPpM.*expBbM.*expTtMa.*expTtMb;
            RecMX=squeeze(sum(RecM,CS_Dim));
            RecMXSDR(:,:,:,DynI,RunI)=RecMX;
        end
    end
end
%%
ZZ=loadniidata('/autofs/cluster/kawin/Gilad/EPTI_fmri_SE/zstat/zstat1_echo1_gradient.nii');
%%
Design_sec=[zeros(1,18) repmat([ones(1,35) zeros(1,53)],[1 3])];
TimePerVol=8.83;
Design_points=interp1(0:(numel(Design_sec)-1),Design_sec,0:TimePerVol:numel(Design_sec))>0.5;
PDX=rot90(Maps(:,:,1,:,:),3);
Squashed=grmss(PDX,3:5);
Msk=Squashed>0.2;
RGB=repmat(Squashed,[1 1 3]);
RGB=RGB./grms(RGB,1:2)/6;
ZZ1=rot90(ZZ,3);
ZZc=ZZ1(:,:,6)>2;
ZZc(:,[1:110 140:end])=0;
ZZc([1:50 70:end],:)=0;
ZZc=ZZc.*Msk;
RGBx=RGB;
ZZcf=flip(ZZc,2);
RGBx(:,:,1)=ZZc;
RGBx(:,:,3)=ZZcf;
figure;imagesc(RGBx);axis equal;removeTicks;
% TwoD=Reshape4d22d(PDX(:,:,:,:,1),ZZcf);
%%
RecMXSDRx=rot90(RecMXSDR,3);
%%
EchosToShow=[3 15 29 40];

figure;
subplot(1,2,1);
Clr=[0.3010 0.7450 0.9330];
rectangle('Position',[3,0,4,100],'FaceColor',Clr);hold on;
rectangle('Position',[13,0,4,100],'FaceColor',Clr);
rectangle('Position',[23,0,4,100],'FaceColor',Clr);
TClr='krgb';
WhichRunsToUse=1;
for i=1:numel(EchosToShow)
    TwoD=Reshape4d22d((mean(abs(RecMXSDRx(:,:,EchosToShow(i),:,WhichRunsToUse)),5)),ZZc);
    TwoDG=Reshape4d22d((mean(abs(GReconXx(:,:,EchosToShow(i),:,WhichRunsToUse)),5)),ZZc);
    h(i)=plot(mean(TwoD,1),TClr(i),'LineWidth',2);
    plot(mean(TwoDG,1),TClr(i),'LineWidth',1);
    Ttl{i}=['Echo #' num2str(EchosToShow(i))];
end
title('ROI');
axis([1 32 0.2 0.5] );
% figure;plot((TwoDx.'));hold on;plot(Design_points*0.05+0.3);
subplot(1,2,2);
rectangle('Position',[3,0,4,100],'FaceColor',Clr);hold on;
rectangle('Position',[13,0,4,100],'FaceColor',Clr);
rectangle('Position',[23,0,4,100],'FaceColor',Clr);
for i=1:numel(EchosToShow)
    TwoDf=Reshape4d22d((mean(abs(RecMXSDRx(:,:,EchosToShow(i),:,WhichRunsToUse)),5)),ZZcf);
    TwoDfG=Reshape4d22d((mean(abs(GReconXx(:,:,EchosToShow(i),:,WhichRunsToUse)),5)),ZZcf);        
%     plot(mean(TwoD,1),TClr(i),'LineWidth',2);
    plot(mean(TwoDf,1),[TClr(i) '--'],'LineWidth',2);
    plot(mean(TwoDfG,1),[TClr(i) '--'],'LineWidth',1);
end
title('Contralateral');
axis([1 32 0.2 0.5] );
subplot(1,2,1);
legend(h,Ttl);
%%
fgmontagex(RecMXSDRx(20:end-10,:,29,15,1)); title('Spin echo, single dynamic');
fgmontagex(mean(abs(RecMXSDRx(20:end-10,:,29,:,1)),4)); title('Spin echo, average over a run (32 dynamics)');
%%
% GReconXx
SpinE=mean(abs(GReconXx(:,:,29,:,1:5)),5);
% SpinE=mean(abs(RecMXSDRx(:,:,3,:,1:5)),5);
ValStart=mean(SpinE(:,:,1,1:3),4);
ValEnd=mean(SpinE(:,:,1,end-2:end),4);
W=perm42(linspace(0,1,32));
ValsFit=ValStart.*(1-W)+ValEnd.*W;
SpinEmFit=SpinE-ValsFit;

[~,tScrx]=ttest2(perm41(SpinE(:,:,1,Design_points)),perm41(SpinE(:,:,1,~Design_points)));
% [~,tScrx]=ttest2(perm41(SpinEmFit(:,:,1,Design_points)),perm41(SpinEmFit(:,:,1,~Design_points)));
% [~,tScrx]=ttest2(perm41(SpinEmFit(:,:,1,Design_points)),perm41(SpinEmFit(:,:,1,~Design_points)),'Vartype','unequal');
tScrx=perm41(tScrx);

OverlayMaskOnImage([], grmss(SpinE,3:5), (-log10(tScrx))>2,[1 0 0],0.7)

TwoD=abs(Reshape4d22d(SpinEmFit,ZZc));
%%
GReconX=zeros([168,96,56,32,5]);
for SliI=6
    for RunI=1:5
        for DynI=1:32
            disp([SliI RunI DynI]);
            tmp=load([BaseSEP 'Recon_B0GRAPPA/Recon_RawCalib_SE_Run' num2str(RunI) '_Dyn' num2str(DynI) '_Slice' num2str(SliI) '_SE_Noupdates.mat']);
            
            MM=ifftshift(ifft2c(tmp.recon),2);
            MMC=sum(MM.*conj(CurSens),4);
            GReconX(:,:,:,DynI,RunI)=MMC;
        end
    end
end

GReconXx=rot90(GReconX,3);
%%
mRecMXSDRx=mean(abs(RecMXSDRx),5);

Echo29=squeeze(mRecMXSDRx(:,:,29,:));

QQ=make_nii(mRecMXSDRx);
QQ.hdr.dime.pixdim(5)=8.83;
save_nii(QQ,'AllEchos.nii');