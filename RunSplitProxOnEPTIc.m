%%
% load('ForMPBT.mat','Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
load([OutP 'ForMPBT_Sli' num2str(SliI) '.mat'],'Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
m0a=Elem0{1};
b0a=Elem0{3};
disp('Loaded')
%%
SRange=[21 21];

SSigs=[0.0001 0.1:0.1:29];
m0b=min(m0a,median(m0a(:))*20);
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
fgmontage(SB0x,[-400 400])
% fgmontage(b0a,[-400 400])
% fgmontage(SB0,[-400 400])
%%
nEchos=size(ksp_adj,7);

SensCSMap=BARTS_Aopx.Others{1};

MskRec=BARTS_Aopx.Others{2};
% SensCSMap(1:2,:,:,:)=0;
% SensCSMap(end-1:end,:,:,:)=0;
% SensCSMap(:,1:2,:,:,:)=0;
% SensCSMap(:,end-1:end,:,:,:)=0;
disp('ok');
%% 
CurPrefix='zz';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
BaseSP=['/autofs/space/daisy_002/users/Gilad/' CurPrefix];
BaseFP=['/tmp/' CurPrefix];
LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

nccToUse=15;
disp('Prepare folders');
system('rm /tmp/*.cfl');
system('rm /tmp/*.hdr');
%%
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);
BARTS_Aopx.Others{1}=SensCSMap(:,:,:,1:nccToUse);
SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
disp('got ksp_adj');
UseOS=false;
UseOSexplicit=false;
%%
% SensCSMap is 0
% mask_sampleP is 1
% TSC is 2
% Copy with sens/CS mask, sum over CS, FFT and sample mask
ITS_TSC_Cmnds={'fmac 2 0',['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_TSC_Cmnds);
BARTS_Aopx.Others{3}=perm73(ThroughPlaneDecay);
BARTS_Aopx.Others{1}=SensCSMap(:,:,:,1:nccToUse);
BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});

disp('got ksp_adj Through-plane');
UseOS=false;
UseOSexplicit=false;
%% Simple OS
load('ForMPBT.mat','Elem0');
Elem0a=Elem0;
%%
nccToUse=15;
nccToUse=9;
OSFac=3;
OS_Sz=Sz.*OSFac;

NHalfDiff=(OS_Sz(1)-Sz(1))/2;

ITS_OS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3',['cropc 3 ' num2str(NHalfDiff)],'fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_OS_Cmnds);

% C0=Elem0a{1}.*exp(1i*Elem0a{2});
% C0r=imresizeBySlices(C0,OS_Sz);
% Elem0{1}=abs(C0r)/2;
% Elem0{2}=angle(C0r);
% Elem0{3}=imresizeBySlices(Elem0a{3},OS_Sz);
% Elem0{4}=max(5,min(200,imresizeBySlices(max(5,min(200,abs(Elem0a{4}))),OS_Sz)));
C=MapsT{1}.*exp(1i*MapsT{2});
CExt=imresize(C,OS_Sz);
Elem0{1}=abs(CExt);
Elem0{2}=angle(CExt);
Elem0{3}=imresize(MapsT{3},OS_Sz);
Elem0{4}=max(5,min(200,imresize(max(5,min(200,MapsT{4})),OS_Sz)));

% BARTS_OS_Aopx.Others{1}=imresizeBySlices(BARTS_Aopx.Others{1},OS_Sz); % Sens
% BARTS_OS_Aopx.Others{3}=imresizeBySlices(BARTS_Aopx.Others{3},OS_Sz); % TSC
BARTS_OS_Aopx.Others{1}=imresize(SensCSMap(:,:,:,1:nccToUse),OS_Sz); % Sens
BARTS_OS_Aopx.Others{3}=imresize(BARTS_Aopx.Others{3},OS_Sz); % TSC

% SigToUse=padarray(sig(:,:,:,1:nccToUse,:,:,:),[NHalfDiff NHalfDiff],'both');
SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
% MskRec=padarray(BARTS_Aopx.Others{2},[NHalfDiff NHalfDiff],'both');
MskRec=BARTS_Aopx.Others{2};

BARTS_OS_Aopx.Others{2}=MskRec;
BARTS_OS_Aopx.ImSz16=FillOnesTo16([OS_Sz 1 1 1 1 nEchos]);
BARTS_OS_Aopx.cmd=['linopScript -A ' LS_ScriptFN];

BARTS_Aop=WriteBARTStructToFiles(BARTS_OS_Aopx,BaseFP);
ksp_adj=bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
disp('got ksp_adj OS explicit');
UseOS=false;
UseOSexplicit=true;
%% OS using linop inside code
UseOSexplicit=false;
BaseFP='/autofs/cluster/kawin/Gilad/tmp/';
OSN=240;
OS_Sz=[OSN OSN];
UseOS=true;
NHalfDiff=(OS_Sz(1)-Sz(1))/2;

% OSx=squeeze(imresize(perm32(eye(Sz(1))),[OS_Sz(1) 1])).';
% OSy=squeeze(imresize(perm32(eye(Sz(2))),[OS_Sz(2) 1]));
% OSxP=perm52(OSx);
% OSyP=perm61(OSy);
% writecfl([ToBARTP 'OSx'],repmat(OSxP,[1 Sz(2) 1 1 1]));
% writecfl([ToBARTP 'OSy'],repmat(OSyP,[1 1 1 1 OS_Sz(1) 1]));
% SensOS=sum(OSxP.*SensCSMap(:,:,:,1:nccToUse),1);
% SensOS=sum(OSyP.*SensOS,2);

% ThroughPlaneDecayOS=sum(OSxP.*ThroughPlaneDecay,1);
% ThroughPlaneDecayOS=sum(OSyP.*ThroughPlaneDecayOS,2);
SensOS=imresize(SensCSMap(:,:,:,1:nccToUse),OS_Sz);
ThroughPlaneDecayOS=imresize(ThroughPlaneDecay,OS_Sz);
% SensCSMap is 0
% mask_sampleP is 1
% TSC is 2
% Copy with sens/CS mask, sum over CS, FFT and sample mask
% ITS_TSC_Cmnds_OS={'fmac 2 0',['fmac 0 ' num2str(CS_Flag)],'fftc 3',['cropc 3 ' num2str(NHalfDiff)],'fmac 1 0'};
% ITS_TSC_Cmnds_OS={'fmac 2 0','fmac 0 0','fftc 48',['cropc 48 ' num2str(NHalfDiff)],'fmac 1 0'};
ITS_TSC_Cmnds_OS={'fmac 2 0','fmac 0 0','fftc 3',['cropc 3 ' num2str(NHalfDiff)],'fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_TSC_Cmnds_OS);
BARTS_Aopx.Others{3}=perm73(ThroughPlaneDecayOS);
BARTS_Aopx.Others{1}=SensOS;
% BARTS_Aopx.Others{2}=permute(MskRec,[5 6 3 4 1 2 7]);
BARTS_Aopx.Others{2}=MskRec;

BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
% ImSz16OS=FillOnesTo16([1 1 1 1 OS_Sz nEchos]);
ImSz16OS=FillOnesTo16([OS_Sz 1 1 1 1 nEchos]);
BARTS_Aop.ImSz16=ImSz16OS;
% SigToUse=permute(sig(:,:,:,1:nccToUse,:,:,:),[5 6 3 4 1 2 7]);
SigToUse=sig(:,:,:,1:nccToUse,:,:,:);
ksp_adj=bart(['linopScript -d 4 -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});

disp('got ksp_adj Through-plane + OS');
%% with OS for in-plane
% nDbls=1;
% OS_Sz=Sz.*(2^nDbls);
% OS_Cmnds=[{'fftc 3'},repmat({'dblszc 3'},[1 nDbls]),{'ifftc 3'}];
% OSScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC_OS.txt'];
% OS_ITS_TSC_Cmnds=[OS_Cmnds,ITS_TSC_Cmnds];
% WriteLinopToFile(OSScriptFN,OS_ITS_TSC_Cmnds);
% disp('Wrote OS LinopScript');
% 
% BARTS_OS_Aopx.cmd=['linopScript -N ' OSScriptFN];
% BARTS_OS_Aopx.ImSz16=FillOnesTo16([Sz 1 1 1 1 nEchos]);
% 
% BARTS_OS_Aopx.Others{1}=imresizeBySlices(BARTS_Aopx.Others{1},OS_Sz); % Sens
% BARTS_OS_Aopx.Others{3}=imresizeBySlices(BARTS_Aopx.Others{3},OS_Sz); % TSC
% 
% NHalfDiff=(OS_Sz(1)-Sz(1))/2;
% 
% SigToUse=padarray(sig,[NHalfDiff NHalfDiff],'both');
% MskRec=padarray(BARTS_Aopx.Others{2},[NHalfDiff NHalfDiff],'both');
% 
% BARTS_OS_Aopx.Others{2}=MskRec;
% 
% BARTS_Aop=WriteBARTStructToFiles(BARTS_OS_Aopx,BaseFP);
% ksp_adj=bart(['linopScript -A ' OSScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
% disp('got ksp_adj');
%% Normalize
Elem0=MapsT;
Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTt0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
% Rec0X=squeeze(sum(Rec0,CS_Dim));

% psig=bart(['linopScript ' LS_ScriptFN],BARTS_Aop.ImSz16,Rec0,BARTS_Aop.Others{:});
AHA0=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Rec0,BARTS_Aop.Others{:});
% [grmss(ksp_adj) grmss(AHA0)]
NFac=grmss(ksp_adj)/grmss(AHA0)
Elem0{1}=Elem0{1}*NFac;
disp(['Normalized ' num2str(NFac) ' ' num2str(grmss(ksp_adj)) ' ' num2str(grmss(AHA0))]);
if(UseOS)
%     Elem0{1}=Elem0{1}/2.5720;
end
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

% ElemsLambda=[0.02,0.05,0.05,0.005]/2; % for OS2
% ElemsLambda=[0.01,0.1,0.05,0.0005]; % for OS2
% ElemsLambda=[0.01,0.1,0.0005,0.0005]; % for OS2
% ElemsLambda=[0.02,0.05,0.05,0.005]/3; % for OS3

for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
    if(UseOS)
%         writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[1 1 1 1 OS_Sz 1]));
        writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
    else
        writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
    end
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
%     ninneriterBART(i)=1;
end
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%% Continue
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
end
%%
if(UseOS)
    ImSz16Base=FillOnesTo16([Sz 1 1 1 1 nEchos]);
    QQ=bart(['splitProx -i 10000 -s 200 -d 2 -o ' num2str(OS_Sz(1))  ' -g -F ' ToBARTP ' ' LS_ScriptFN],ImSz16Base,1,BARTS_Aop.Others{:});
else
    QQ=bart(['splitProx -i 10000 -s 200 -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
end
% QQ=bart(['splitProx -i 1000 -s 60 -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
% QQ=bart(['splitProx -i 500 -s 60 -d 2 -g -F ' ToBARTP ' ' OSScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
% QQ=bart(['splitProx -i 500 -s 60 -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,PaddedData,BARTS_Aop.Others{:});
%%
tmp=readcfl([ToBARTP 'ElemsWS_1']);

tmp=readcfl([ToBARTP 'mDimsOut']);
%%
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
figure;plot(ErrVec)
%%
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);

for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
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
Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTt0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
Rec0X=squeeze(sum(Rec0,CS_Dim));
disp('WS maps');
%%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Elem0{i},ElemRanges{i});removeTicks;
    if(i==1), title('Start point'); end
end
%%
SpiRefPrefix='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1/Sli12_RA1_2_3A_';
SpiRefMaps=load([SpiRefPrefix 'Maps.mat']);
SpiRefMaps=SpiRefMaps.Maps;
SpiRefElemL3=readcfl([SpiRefPrefix 'ElemsL_3']);
SpiRefElem0=readcfl([SpiRefPrefix 'Elem0']);
SpiRefElem3=readcfl([SpiRefPrefix 'Elem3']);

SpiRef=SpiRefElem0.*exp(SpiRefElemL3./SpiRefElem3);
disp('Loaded');
%%
fgmontage(Maps{3}-Elem0{3},[-5 5])
%%
EchoIToShow=30;
% fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec_LLRX(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)),'Size',[1 3]);
fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)));
fgmontagex(angle(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow))));

fgmontagex(RecAllEchos(:,:,EchoIToShow));
fgmontagex(Rec0X(:,:,EchoIToShow),[0 10]);title('GRAPPA Fitted');
fgmontagex(RecMX(:,:,EchoIToShow),[0 10]);title('splitProx');
fgmontagex(RecMX(:,:,EchoIToShow),[0 10]);title('splitProx+Through plane');
fgmontagex(RecMX(:,:,EchoIToShow),[0 5]);title('splitProx+Through plane+InPlane');

fgmontagex(Rec0X(:,:,1:8:end),[0 10]);title('GRAPPA Fitted');
fgmontagex(RecMX(:,:,1:8:end),[0 10]);title('splitProx');
fgmontagex(RecMX(:,:,1:8:end),[0 10]);title('splitProx+Through plane');
fgmontagex(RecMX(:,:,1:8:end),[0 10/OSFac]);title('splitProx+Through plane+in-plane');

fgmontagex(RecAllEchosS(:,:,5:13:end,SliI),[0 4]);title('GRAPPA');
fgmontagex(FittedSG(:,:,5:13:end,SliI),[0 4]);title('Fitted GRAPPA');
fgmontagex(RecMX(:,:,5:13:end),[0 4]);title('splitProx');

fgmontagex(RecMX(:,:,1:13:end),[0 10]);title('splitProx+Through plane');
fgmontagex(RecMX(:,:,1:13:end),[0 10/OSFac]);title('splitProx+Through plane+in-plane');
%%
fgmontagex(G6.im_EPTI_correct(:,:,1,1:8:end));title('GRAPPA');
fgmontagex(S6.im_recon(:,:,1:8:end));title('Subspace');

fgmontagex(G6.im_EPTI_correct(:,:,1,EchoIToShow));title('GRAPPA');
fgmontagex(S6.im_recon(:,:,EchoIToShow));title('Subspace');
fgmontagex(G6.im_EPTI_correct(:,:,1,1:13:end),[0 10]);title('GRAPPA');
fgmontagex(S6.im_recon(:,:,1:13:end),[0 10]);title('Subspace');
%%
SeveralRecs=cat(3,G6.im_EPTI_correct(:,:,1,EchoIToShow),S6.im_recon(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow),RecMX(:,:,EchoIToShow));
fgmontagex(SeveralRecs,[0 7]);
%%
SeveralRecsFull=cat(4,squeeze(G6.im_EPTI_correct),S6.im_recon,Rec0X,RecMX);
for i=1:4
        [PDBaseSR(:,:,i), UpdatedB0Map_HzSR(:,:,i), UpdatedT2SMap_msSR(:,:,i), s_valsSR(:,:,:,i),...
        ~, PDBase0SR(:,:,i)]=...
        FitToModel_MPBD1CSf(SeveralRecsFull(:,:,:,i),WhichEchosToUse,ES_ms,FirstTE_ms);
end
fgmontagex(UpdatedT2SMap_msSR,[0 200]);colormap hot
%%
fgmontage(G16.im_EPTI_correct(:,:,1,11:11:end))
fgmontage(S16.im_recon(:,:,11:11:end))
%%
S6Rec={S6Grappa S6Subspace S6SplitProx S6SplitProxT S6SplitProxTI};
Ttls={'GRAPPA' 'Subspace' 'SplitProx' 'SplitProx+Throuh plane' '+oversampling'};
RR=[1 1 1 1 2];
% tight_subplot(Nh, Nw, gap, marg_h, marg_w)
figure;
ha = tight_subplot(2,3,[.0 .0],[.01 .0],[.0 .0]);
for ii = 1:5
    axes(ha(ii));
    gmontage(S6Rec{ii},[0 10]/RR(ii));title(Ttls{ii});
    removeTicks;daspect([1 1 1]);
end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

%%
save('EPTI_splitProx_recons.mat','Rec0X','RecMX','Maps','Maps1','Maps50');