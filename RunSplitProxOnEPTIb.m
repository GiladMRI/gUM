%%
load('ForMPBT.mat','Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
nEchos=size(ksp_adj,7);
nElements=numel(Elem0);
%% 
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/aa';

BaseSP='/autofs/space/daisy_002/users/Gilad/aa';

BaseFP='/tmp/aa';
%%
LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% TSC is 2
% Copy with sens/CS mask, sum over CS, FFT and sample mask
ITS_TSC_Cmnds={'fmac 2 0',['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_TSC_Cmnds);

BARTS_Aopx.Others{3}=perm73(ThroughPlaneDecay);

BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
disp('got ksp_adj resular');
%% with OS for in-plane
nDbls=2;
OS_Sz=Sz.*(2^nDbls);
OS_Cmnds=[{'fftc 3'},repmat({'dblszc 3'},[1 nDbls]),{'ifftc 3'}];
OSScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC_OS.txt'];
OS_ITS_TSC_Cmnds=[OS_Cmnds,ITS_TSC_Cmnds];
WriteLinopToFile(OSScriptFN,OS_ITS_TSC_Cmnds);
disp('Wrote OS LinopScript');

BARTS_OS_Aopx.cmd=['linopScript -N ' OSScriptFN];
BARTS_OS_Aopx.ImSz16=FillOnesTo16([Sz 1 1 1 1 nEchos]);

BARTS_OS_Aopx.Others{1}=imresizeBySlices(BARTS_Aopx.Others{1},OS_Sz); % Sens
BARTS_OS_Aopx.Others{3}=imresizeBySlices(BARTS_Aopx.Others{3},OS_Sz); % TSC

NHalfDiff=(OS_Sz(1)-Sz(1))/2;

PaddedData=padarray(sig,[NHalfDiff NHalfDiff],'both');
MskRec=padarray(BARTS_Aopx.Others{2},[NHalfDiff NHalfDiff],'both');

BARTS_OS_Aopx.Others{2}=MskRec;

BARTS_Aop=WriteBARTStructToFiles(BARTS_OS_Aopx,BaseFP);
ksp_adj=bart(['linopScript -A ' OSScriptFN],BARTS_Aop.ImSz16,PaddedData,BARTS_Aop.Others{:});
disp('got ksp_adj');
%%



%%
ElemsAlphas=[0.0323, 1.5366e-07, 1.3036e-05, 0.0015];
ElemsAlphas=[0.323, 1e-03, 1e-1, 1e-0];
ElemsLambda=[0.015,0.05,0.05,0.0003]; % Very good for each one separately
ElemsLambda=[0.02,0.05,0.05,0.005];

for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[Sz 1 1 1 1 1]));
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
end
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%% Continue?
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
end
%%
% QQ=bart(['splitProx -i 500 -s 60 -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
QQ=bart(['splitProx -i 500 -s 60 -d 2 -g -F ' ToBARTP ' ' OSScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
%%
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
figure;plot(ErrVec)
%%
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
%%
Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTt0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
Rec0X=squeeze(sum(Rec0,CS_Dim));
disp('WS maps');
%%
ElemRanges={[0 10],[-pi pi],[-400 400],[0 100]};
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); end
end
%%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Elem0{i},ElemRanges{i});removeTicks;
    if(i==1), title('Start point'); end
end
%%
fgmontage(Maps{3}-Elem0{3},[-5 5])
%%
EchoIToShow=30;
% fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec_LLRX(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)),'Size',[1 3]);
fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)));
fgmontagex(angle(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow))));

fgmontagex(RecAllEchos(:,:,EchoIToShow));
fgmontagex(Rec0X(:,:,EchoIToShow),[0 10]);title('GRAPPA Fitted');
fgmontagex(RecMX(:,:,EchoIToShow),[0 10]);title('splitProx+Through plane');

fgmontage(Rec0X(:,:,1:8:end),[0 10]);title('GRAPPA Fitted');
fgmontage(RecMX(:,:,1:8:end),[0 10]);title('splitProx+Through plane');
%%
save('EPTI_splitProx_recons.mat','Rec0X','RecMX','Maps','Maps1','Maps50');