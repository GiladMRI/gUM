SliI=7;
%%
for SliI=1:nSlices
    for dyn=1:4
% nEchos=size(ksp_adj,7);
SensCSMap=SensCSMapS(:,:,SliI,:);
MskRec=mask_sample;
sig=single(MskCC);
sig(MskCC)=SigMskd(:,SliI,dyn);
sig=perm75(sig);

PDBase0=PDBase0SG(:,:,SliI);
UpdatedB0Map_Hz=UpdatedB0Map_HzSG(:,:,SliI);
UpdatedT2SMap_ms=UpdatedT2SMap_msSG(:,:,SliI);

m0=abs(PDBase0);
m0=min(m0,getPercentile(m0(:),0.9)*2);
p0=angle(PDBase0);
b0=UpdatedB0Map_Hz;
t0=max(5,min(200,abs(UpdatedT2SMap_ms)));

Elem0={m0,p0,b0,t0};
%
% load('ForMPBT.mat','Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
% load([OutP 'ForMPBT_Sli' num2str(SliI) '.mat'],'Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
m0a=Elem0{1};
b0a=Elem0{3};
disp('Loaded')
%%
SRange=[21 21];

SSigs=[0.0001 0.1:0.1:29];
m0b=m0a; 
M0B0=m0b.*b0a;

% SM0B0=SmoothBySlices(M0B0,SRange,SSig);
% SM0=SmoothBySlices(m0b,SRange,SSig);
% SB0=SM0B0./SM0;
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
CurSPrefix=[OutP 'Sli' num2str(SliI) '_'];
CurSDPrefix=[OutP 'Sli' num2str(SliI) '_dyn' num2str(dyn) '_'];

CurPrefix='zz';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

nccToUse=15;
disp('Prepare folders');
% system('rm /tmp/*.cfl');
% system('rm /tmp/*.hdr');
%%
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap(:,:,:,1:nccToUse) mask_sample};

SigToUse=sig(:,:,:,1:nccToUse,:,:,:);

BARTS_Aop=BARTS_Aopx;

writecfl([OutP 'mask_sample'],mask_sample);
BARTS_Aop.Others{2}=[OutP 'mask_sample'];

writecfl([CurSPrefix 'SensCC'],SensCSMap(:,:,:,1:nccToUse));
BARTS_Aop.Others{1}=[CurSPrefix 'SensCC'];

SigFN=[CurSDPrefix 'Sig'];
writecfl(SigFN,SigToUse);

ksp_adjFN=[CurSDPrefix 'sig_adj'];

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
writecfl(ksp_adjFN,ksp_adj);

disp('got ksp_adj');
UseOS=false;
UseOSexplicit=false;
%
% ThroughPlaneDecayFN=[CurSPrefix 'ThroughPlaneDecay'];
% writecfl(ThroughPlaneDecayFN,perm73(ThroughPlaneDecayS(:,:,:,SliI)));
%%
Elem0{1}=Elem0{1}*0;
%
Elem0{2}([1:2 end-1:end],:)=0;
Elem0{2}(:,[1:2 end-1:end])=0;
Elem0{3}=SB0x;
% Elem0{3}=b0a;
Elem0{4}=max(5,min(200,abs(Elem0{4})));
disp('Improved initials');
%%
ElemsAlphas=[0.0323, 1.5366e-07, 1.3036e-05, 0.0015];
ElemsAlphas=[0.323, 1e-03, 1e-1, 1e-0];
ElemsLambda=[0.015,0.05,0.05,0.0003]; % Very good for each one separately
ElemsLambda=[0.02,0.05,0.05,0.005];
ElemsLambda=[0.002,0.05,0.05,0.0005];
ElemsLambda=[0.002,0.05,0.05,0.0005];
ElemsLambda=[0.002,0.5,0.01,0.001];

% OS:
if(UseOS)
    ElemsAlphas=[1e-5, 1e-05, 1e-3, 1e-3];
    ElemsLambda=[1e-3,0.5,0.01,1e-4];
end
%
for i=1:numel(Elem0)
    writecfl([CurSPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
    if(UseOS)
        writecfl([CurSPrefix 'ElemsLOS_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
    else
        writecfl([CurSPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
    end
end

writecfl([CurSPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSPrefix 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([CurSPrefix 'ElementTypes'],ElementTypes.');

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
%     ninneriterBART(i)=1;
end
writecfl([CurSPrefix 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([CurSPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue
% for i=1:numel(Elem0)
%     writecfl([CurSPrefix 'ElemsWS_' num2str(i-1)],Maps{i});
% end
%%
if(UseOS)
    ImSz16Base=FillOnesTo16([Sz 1 1 1 1 nEchos]);
    BARTCmd{SliI,dyn}=bartCmd(['splitProx -i 10000 -s 200 -d 2 -o ' num2str(OS_Sz(1))  ' -g -F ' ToBARTP ' ' LS_ScriptFN],ImSz16Base,1,BARTS_Aop.Others{:});
else
    BARTCmd{SliI,dyn}=bartCmd(['splitProx -i 10000 -s 60 -d 2 -g -f -F ' CurSPrefix ' -O ' CurSDPrefix ' -S ' ksp_adjFN ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
end
    end % dyns
end % slices
%%
% BARTCmds=BARTCmd.';
% BARTCmds=BARTCmd(SliI,:).';
BARTCmds=BARTCmd(:);
BartCmdFN=[OutP 'splitProxRun'];
fid=fopen(BartCmdFN,'w+');
fprintf(fid,'# splitProx run\r\n');
fprintf(fid,'echo "Starting!"\r\n');
for i=1:numel(BARTCmds)
    fprintf(fid,'%s\n',BARTCmds{i});
end
fclose(fid);
[status,msg,msgID] = fileattrib(BartCmdFN,'+x');
disp(['Prepared script ' BartCmdFN]);
%%
QQ=system(BARTCmd{SliI});



%%
ErrVec=readcfl([CurSDPrefix 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
figure;plot(ErrVec)
%%
ErrVec=readcfl([CurSDPrefix 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);

for i=1:4
    Maps{i}=readcfl([CurSDPrefix 'Elem' num2str(i-1)]);
end

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtM = exp(ElemsL{4} .* (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));

disp('Loaded maps');
%
if(UseOS || UseOSexplicit)
    OSFac=OS_Sz(1)/Sz(1);
else
    OSFac=1;
end
ElemRanges={[0 10]/OSFac,[-pi pi],[-400 400],[0 100]};
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); xlabel([num2str(numel(ErrVec)/sum(ninneriterBART)) ' : ' num2str(ErrVec(end),'%.7g')]); end
    if(i==3), xlabel(num2str(ElemsAlphas,' %.9g,')); end %ylabel(Pref,'Interpreter','None'); end
    if(i==4), xlabel(num2str(ElemsLambda,' %.9g,')); end
end
%%
SliI=3;
for dyn=1:4
CurSPrefix=[OutP 'Sli' num2str(SliI) '_'];
CurSDPrefix=[OutP 'Sli' num2str(SliI) '_dyn' num2str(dyn) '_'];
disp(CurSDPrefix);
for i=1:4
    Maps{i}=readcfl([CurSDPrefix 'Elem' num2str(i-1)]);
end

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtM = exp(ElemsL{4} .* (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMXD(:,:,:,dyn)=RecMX;
end
%%
QQ=squeeze(RecMXD(:,:,50,:));
mQQ=mean(abs(QQ),3);
sQQ=std(abs(QQ),[],3);
tSNR=mQQ./sQQ;
fgmontagex(gflip(tSNR,2),[0 100]);colorbar;title('tSNR EPTI 3-shot');
figure;histogram(tSNR(mQQ>.5),0:5:200);title('tSNR EPTI 3-shot')
%%
OutP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns/';
nSlices=16;
SliI=1;
CurSPrefix=[OutP 'Sli' num2str(SliI) '_'];
for i=1:4
    ElemsL{i}= readcfl([CurSPrefix 'ElemsL_' num2str(i-1)]);
end
CS_Dim=5;
%%
for SliI=1:nSlices
for dyn=1:4
CurSPrefix=[OutP 'Sli' num2str(SliI) '_'];
CurSDPrefix=[OutP 'Sli' num2str(SliI) '_dyn' num2str(dyn) '_'];
disp(CurSDPrefix);
for i=1:4
    Maps{i}=readcfl([CurSDPrefix 'Elem' num2str(i-1)]);
end

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtM = exp(ElemsL{4} .* (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMXDS(:,:,:,SliI,dyn)=RecMX;
end
end
%%
aRecMXDS=single(abs(RecMXDS));
save([OutP 'aRecMXDS.mat'],'aRecMXDS');
%%
figure;
for SliI=1:nSlices
QQ=squeeze(RecMXDS(:,:,50,SliI,:));
mQQ=mean(abs(QQ),3);
sQQ=std(abs(QQ),[],3);
tSNR=mQQ./sQQ;
tSNRS(:,:,SliI)=tSNR;
% fgmontagex(gflip(tSNR,2),[0 100]);colorbar;title('tSNR EPTI 3-shot');
subplot(4,4,SliI);
HMsk(:,:,SliI)=mQQ>.5;
h=histogram(tSNR(HMsk(:,:,SliI)),0:5:200);title('tSNR EPTI 3-shot')
counts{SliI} = h.Values;
end
CountsM=cat(1,counts{1:14});
CountsV=sum(CountsM,1);
CountsV=CountsV./sum(CountsV);
%%
SKEPTIC_Counts=load('SKEPTIC_Counts.mat','counts');
SKEPTIC_Counts=SKEPTIC_Counts.counts;
SKEPTIC_CountsM=cat(1,SKEPTIC_Counts{1:14});
SKEPTIC_CountsV=sum(SKEPTIC_CountsM,1);
SKEPTIC_CountsV=SKEPTIC_CountsV./sum(SKEPTIC_CountsV);
%%
Xs=5:5:200;
figure;
plot(Xs,CountsV,'k');hold on;
plot(Xs,SKEPTIC_CountsV,'r');
legend({'EPTI 3-shot','SKEPTIC 3-shot'})
title(['EPTI ' num2str((5:5:200)*CountsV.','%.2f') ' SKEPTIC ' num2str((5:5:200)*SKEPTIC_CountsV.','%.2f')]);
%%
%%
figure;
for SliI=1:nSlices
QQ=squeeze(RecMXDS(:,:,50,SliI,:));
mQQ=mean(abs(QQ),3);
sQQ=std(abs(QQ),[],3);
tSNR=mQQ./sQQ;
tSNRS(:,:,SliI)=tSNR;
% fgmontagex(gflip(tSNR,2),[0 100]);colorbar;title('tSNR EPTI 3-shot');
subplot(4,4,SliI);
HMsk(:,:,SliI)=mQQ>.5;
h=histogram(tSNR(HMsk(:,:,SliI)),0:5:200);title('tSNR EPTI 3-shot')
counts{SliI} = h.Values;
end

%%
SmallVarT=50e6;
getSmallVars;
save([OutP 'SmallVars.mat'],SmallVars{:});
%%
EPTIRecP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Recon/Human/';

G6=load([EPTIRecP 'Recon_EPTI_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_6_GE.mat']);
G16=load([EPTIRecP 'Recon_EPTI_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_16_GE.mat']);
S6=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_6_GE.mat']);
S16=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_16_GE.mat']);
%%
FirstT_ms=9;
ES=0.93;
TE0_ms=0;
WhichTSToUs=11:60;
[PDBaseG6, UpdatedB0MapG6, UpdatedT2SMap_msG6, s_valsG6, Fitted0G6, PDBase0G6]=...
    FitToModel_MPBD1CSf(G6.im_EPTI_correct,WhichTSToUs,ES,TE0_ms+FirstT_ms);

[PDBaseG16, UpdatedB0MapG16, UpdatedT2SMap_msG16, s_valsG16, Fitted0G16, PDBase0G16]=...
    FitToModel_MPBD1CSf(G16.im_EPTI_correct,WhichTSToUs,ES,TE0_ms+FirstT_ms);

[PDBaseS6, UpdatedB0MapS6, UpdatedT2SMap_msS6, s_valsS6, Fitted0S6, PDBase0S6]=...
    FitToModel_MPBD1CSf(S6.im_recon,WhichTSToUs,ES,TE0_ms+FirstT_ms);

[PDBaseS16, UpdatedB0MapS16, UpdatedT2SMap_msS16, s_valsS16, Fitted0S16, PDBase0S16]=...
    FitToModel_MPBD1CSf(S16.im_recon,WhichTSToUs,ES,TE0_ms+FirstT_ms);
%%
Maps6=Maps;

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtM = exp(ElemsL{4} .* (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMX6=RecMX;
%%
Maps16=Maps;

MmM=ElemsL{1}.*Maps{1};
expPpM = exp(ElemsL{2} .* Maps{2});
expBbM = exp(ElemsL{3} .* Maps{3});
expTtM = exp(ElemsL{4} .* (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));

RecMX16=RecMX;
%%
Recs6=cat(4,squeeze(G6.im_EPTI_correct),Fitted0G6,S6.im_recon,Fitted0S6,RecMX6);
Recs6=Recs6./grms(Recs6,1:3);
fgmontagex(Recs6(:,:,11:31:end,:));caxis(caxis/1.5);
fgmontagex(rot90(Recs6(:,70:end,50,:),3));caxis(caxis/1);

T2S=cat(3,UpdatedT2SMap_msG6,UpdatedT2SMap_msS6,Maps{4});
fgmontagex(rot90(T2S,3),[0 200],'Size',[1 3]);colormap hot
%%
Recs16=cat(4,squeeze(G16.im_EPTI_correct),Fitted0G16,S16.im_recon,Fitted0S16,RecMX16);
Recs16=Recs16./grms(Recs16,1:3);
fgmontagex(Recs16(:,:,11:31:end,:));caxis(caxis/1.5);
fgmontagex(rot90(Recs16(:,70:end,50,:),3));caxis(caxis/1);

T2S=cat(3,UpdatedT2SMap_msG16,UpdatedT2SMap_msS16,Maps{4}*1.5);
fgmontagex(rot90(T2S,3),[0 300],'Size',[1 3]);colormap hot
%%
fgmontage(G6.im_EPTI_correct(:,:,1,11:11:end));title('GRAPPA');
fgmontage(S6.im_recon(:,:,11:11:end));title('Subspace');
fgmontage(RecMX(:,:,11:11:end));title('SplitProx');
%%



%% 
nShots=3;
MskCCShot=zeros([size(MskCC) nShots]);
ShotLocs={1:40,41:80,81:120};
for i=1:nShots
    MskCCShot(:,ShotLocs{i},:,:,:,i)=MskCC(:,ShotLocs{i},:,:,:);
end
MskCCShot=MskCCShot>0;
% IFMskCCShot=fft1cg(MskCCShot(:,:,:,1:nccToUse,:,:,:,:),2);
% IFMskCCShot=fft1cg(MskCCShot(1,:,:,1,:,:,:,:),2);
IFMskCCShot=fft1cg(MskCCShot(:,:,:,1,:,:,:,:),2);
IFMskCCShot=perm75(IFMskCCShot);

SensCSMap=SensCSMapS(:,:,SliI,:);
MskRec=mask_sample;
sig=single(MskCC);
sig(MskCC)=SigMskd(:,SliI);

sigShot=zeros([size(MskCC) nShots]);
for i=1:nShots
    sigShot(:,ShotLocs{i},:,:,:,i)=sig(:,ShotLocs{i},:,:,:);
end
sigShot=perm75(sigShot);
sigShot=sum(sigShot,2);

IFsigShot=ifft1cg(sigShot(:,:,:,1:nccToUse,:,:,:,:,:,:),1);

CurPrefix='zz';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
EPTI_ITS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC_NoPEFT.txt'];

% 0: Sens map
% 1: F map
ITS_Ops={'fmac 0 0','fmac 1 2'};
WriteLinopToFile(EPTI_ITS_ScriptFN,ITS_Ops);

% ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
ksp_adjx=bart(['linopScript -A -d 5 ' EPTI_ITS_ScriptFN],BARTS_Aop.ImSz16,IFsigShot,BARTS_Aop.Others{1},IFMskCCShot);






SensCSMap=SensCSMapS(:,:,SliI,:);
MskRec=mask_sample;
sig=single(MskCC);
sig(MskCC)=SigMskd(:,SliI);
sig=perm75(sig);

PDBase0=PDBase0SG(:,:,SliI);
UpdatedB0Map_Hz=UpdatedB0Map_HzSG(:,:,SliI);
UpdatedT2SMap_ms=UpdatedT2SMap_msSG(:,:,SliI);

m0=abs(PDBase0);
m0=min(m0,getPercentile(m0(:),0.9)*2);
p0=angle(PDBase0);
b0=UpdatedB0Map_Hz;
t0=max(5,min(200,abs(UpdatedT2SMap_ms)));

Elem0={m0,p0,b0,t0};
%
% load('ForMPBT.mat','Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
% load([OutP 'ForMPBT_Sli' num2str(SliI) '.mat'],'Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
m0a=Elem0{1};
b0a=Elem0{3};
disp('Loaded')
%%
SRange=[21 21];

SSigs=[0.0001 0.1:0.1:29];
m0b=m0a; 
M0B0=m0b.*b0a;

% SM0B0=SmoothBySlices(M0B0,SRange,SSig);
% SM0=SmoothBySlices(m0b,SRange,SSig);
% SB0=SM0B0./SM0;
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
CurSPrefix=[OutP 'Sli' num2str(SliI) '_'];

CurPrefix='zz';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

nccToUse=15;
disp('Prepare folders');
system('rm /tmp/*.cfl');
system('rm /tmp/*.hdr');
%%
% ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 1','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

NewMethod=false;
if(NewMethod)
    LS_ScriptFN=EPTI_ITS_ScriptFN;
end



BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap(:,:,:,1:nccToUse) mask_sample};





SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
% SigToUse=ifft1cg(SigToUse,1);

if(NewMethod)
    BARTS_Aopx.Others{2}=IFMskCCShot;
    SigToUse=IFsigShot;
end

BARTS_Aop=BARTS_Aopx;

writecfl([OutP 'mask_sample'],mask_sample);
BARTS_Aop.Others{2}=[OutP 'mask_sample'];



if(NewMethod)
    writecfl([OutP 'IFMskCCShot'],IFMskCCShot);
    BARTS_Aop.Others{2}=[OutP 'IFMskCCShot'];
end



writecfl([CurSPrefix 'SensCC'],SensCSMap(:,:,:,1:nccToUse));
BARTS_Aop.Others{1}=[CurSPrefix 'SensCC'];

SigFN=[CurSPrefix 'Sig'];
writecfl(SigFN,SigToUse);

ksp_adjFN=[CurSPrefix 'sig_adj'];

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
writecfl(ksp_adjFN,ksp_adj);

disp('got ksp_adj');
UseOS=false;
UseOSexplicit=false;
%
ThroughPlaneDecayFN=[CurSPrefix 'ThroughPlaneDecay'];
writecfl(ThroughPlaneDecayFN,perm73(ThroughPlaneDecayS(:,:,:,SliI)));
%%
Elem0{1}=Elem0{1}*0;
%
Elem0{2}([1:2 end-1:end],:)=0;
Elem0{2}(:,[1:2 end-1:end])=0;
Elem0{3}=SB0x;
% Elem0{3}=b0a;
Elem0{4}=max(5,min(200,abs(Elem0{4})));
disp('Improved initials');
%%
ElemsAlphas=[0.0323, 1.5366e-07, 1.3036e-05, 0.0015];
ElemsAlphas=[0.323, 1e-03, 1e-1, 1e-0];
ElemsLambda=[0.015,0.05,0.05,0.0003]; % Very good for each one separately
ElemsLambda=[0.02,0.05,0.05,0.005];
ElemsLambda=[0.002,0.05,0.05,0.0005];
ElemsLambda=[0.002,0.05,0.05,0.0005];
ElemsLambda=[0.002,0.5,0.01,0.001];
% ElemsLambda=[0.02,0.5,0.05,0.005];
% ElemsLambda=[0.02,0,0.05,0.005];
% OS in main linop (explicit OS 2)
% ElemsAlphas=[1e-1, 1e-03, 1e-1, 1e-0];
% ElemsLambda=[1e-4,1e-1,1e-2,1e-3];
% 
% % OS in main linop (explicit OS 3)
% ElemsAlphas=[1e-2, 1e-03, 1e-1, 1e-0];
% ElemsLambda=[1e-4,1e-1,1e-2,1e-3];
% ElemsAlphas=[1e-2, 1e-03, 1e-1, 1e-0];
% ElemsLambda=[1e-3,1e-1,1e-2,1e-4];

% OS:
if(UseOS)
    ElemsAlphas=[1e-5, 1e-05, 1e-3, 1e-3];
    ElemsLambda=[1e-3,0.5,0.01,1e-4];
end
%
for i=1:numel(Elem0)
    writecfl([CurSPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
    if(UseOS)
        writecfl([CurSPrefix 'ElemsLOS_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
    else
        writecfl([CurSPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
    end
end

writecfl([CurSPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSPrefix 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([CurSPrefix 'ElementTypes'],ElementTypes.');

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
%     ninneriterBART(i)=1;
end
writecfl([CurSPrefix 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([CurSPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue
% for i=1:numel(Elem0)
%     writecfl([CurSPrefix 'ElemsWS_' num2str(i-1)],Maps{i});
% end
%%
QQQQ=bartCmd(['splitProx -i 10000 -s 60 -d 2 -g -f -F ' CurSPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
system(QQQQ)

% if(UseOS)
%     ImSz16Base=FillOnesTo16([Sz 1 1 1 1 nEchos]);
%     BARTCmd{SliI}=bartCmd(['splitProx -i 10000 -s 200 -d 2 -o ' num2str(OS_Sz(1))  ' -g -F ' ToBARTP ' ' LS_ScriptFN],ImSz16Base,1,BARTS_Aop.Others{:});
% else
%     BARTCmd{SliI}=bartCmd(['splitProx -i 10000 -s 60 -d 2 -g -f -F ' CurSPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
% end
%%
FNB='/autofs/cluster/kawin/Gilad/EPTI_fmri_SEX/splitProx/splitProxRun_S11_All';
LL=getLines(FNB);
LLX=strrep(LL,'-g -f -F','-g -F');
putLines(FNB,LLX);