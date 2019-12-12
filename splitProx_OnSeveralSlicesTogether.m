nccToUseSeveral=23;
%%
for SliI=1:nSlices
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    
    clear SelfSens1 sccmtx SensCC
    CurSensFN=[SliPrefix 'SelfSens.mat'];
    load(CurSensFN,'SelfSens1');
    sccmtxFN=[SliPrefix 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccToUseSeveral),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');
    SensCCS(:,:,:,:,SliI)=SensCC;
end
disp('Collected sens maps');
% SensCCS=SensCCS(:,:,:,1:nccToUseSeveral,:);
%% 
for SliI=1:nSlices
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    
    sccmtxFN=[SliPrefix 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,1:MaxRepsToUse,:,:,:,:,:,:,:,:);
    DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
    DataCCP=CombineDims(DataCCP,[4 1]);

    DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccToUseSeveral),[3 4 1 2]),3));
    DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,1:MaxRepsToUse);
    DataCCP=permute(DataCCP,[1:3 8 4:7]);
    DataCCPS(:,:,:,:,SliI)=DataCCP;
end
disp('read and CC sig all slices');
% DataCCPS=DataCCPS(:,:,:,1:nccToUseSeveral,:);
%%
for SliI=1:nSlices
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');

for rs=2 % 1:numel(RepSets)
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

        CurSRPrefix=[mainP filesep 'gB0Sli' num2str(SliI) '_R' RepsStr '_'];

% c0=PDBase0THLR_RS(:,:,rs);
% b0=UpdatedB0MapTHLR_RS(:,:,rs);
% t0=UpdatedT2SMap_msTHLR_RS(:,:,rs);

% c0=PDBase0THLR_RS(:,:,rs);
c0=PDBase0THLR_RS(:,:,end);
b0=UpdatedB0MapTHLR_RS(:,:,end);
t0=UpdatedT2SMap_msTHLR_RS(:,:,rs);

c0=min(abs(c0),getPercentile(abs(c0(:)),0.9)*2).*exp(1i*angle(c0));
%
SRange=[21 21];

SSigs=[0.0001 0.1:0.1:29];
m0a=abs(c0);
b0a=b0;
b0a=min(max(b0a,-MaxB0),MaxB0);

m0b=min(m0a,median(m0a(:))*20);
M0B0=m0b.*b0a;

for i=1:numel(SSigs)
    SM0B0(:,:,i)=SmoothBySlices(M0B0,SRange,SSigs(i));
    SM0(:,:,i)=SmoothBySlices(m0b,SRange,SSigs(i));
end
mThresh=median(m0b(m0b>0));
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
%


m0=abs(c0);
% m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

% c0=SmoothBySlices(c0,[20 20],5);
% WW=SmoothBySlices(abs(c0),[20 20],5);

t0=t0*0+50;
t0=min(max(abs(t0),T2SRange(1)),T2SRange(2));
b0=min(max(b0,-MaxB0),MaxB0);
m0=abs(c0);
m0=min(m0,median(m0(:))*M0MedianFac);
p0=angle(c0);
c0=m0.*exp(1i*p0);

b0=SB0x;

m0=m0*0;

Elem0={m0,p0,b0,t0};
Elems0S{SliI,rs}=Elem0;
end
end
disp('Got starting point');
%%
% Mm0=M*m0;
% expPp0 = exp(P * p0);
% expBb0 = exp(B * b0);
% expTt0 = exp(TT * (1./t0));
% 
% Est0=Mm0.*expPp0.*expBb0.*expTt0;

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
ElemL={M,P,B,TT};
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
disp('ok initials');
%% Now specific to a slice set
% SliIs=ROrd(4:7);
% SliIs=ROrd(10:13);
SliIs=ROrd(3:5);
% SliIs=ROrd(4:11);
% SliIs=ROrd(1:16);
nSlicesTogether=numel(SliIs);
% WhichRepsToUseS=repmat([1 13 25],[nSlicesTogether 1]);
WhichRepsToUseS=repmat([1 13],[nSlicesTogether 1]);
% WhichRepsToUseS=WhichRepsToUseS+(((1:nSlicesTogether)-1)*29).';
WhichRepsToUseS=WhichRepsToUseS+(((1:nSlicesTogether)-1)*17).';
WhichRepsToUseS=mod(WhichRepsToUseS-1,36)+1;

WhichRepsToUseS=[25 7; 1 13; 19 31];
RepPhis=(WhichRepsToUseS-1)*110*pi/180;
figure;
CLRS={'r','g','b','k','r--','g--','b--','k--' 'r.-','g.-','b.-','k.-','r-o','g-o','b-o','k-o'};
for i=1:nSlicesTogether
    for j=1:size(WhichRepsToUseS,2)
        plot([0 exp(1i*RepPhis(i,j))+1i*eps],CLRS{i});hold on;
    end
end
title(numel(unique(WhichRepsToUseS)));
%%
% WhichRepsToUseS=[1 13 25; 7 19 31];
% SliIs=ROrd([4 5]);
clear SigToUseS STraj3MS SumKernsPMMedS ksp_adjS
for i=1:numel(SliIs)
    SliI=SliIs(i);
    disp(SliI);
    WhichRepsToUse=WhichRepsToUseS(i,:);
    SigToUse=DataCCPS(:,:,WhichRepsToUse,:,SliI);
    SigToUseS(:,:,:,:,i)=SigToUse;

    BARTS_Aopx.ImSz16=ImSz16FN;
    BARTS_Aopx.Others={SensCCS(:,:,:,:,SliI) STraj3M(:,:,WhichRepsToUse) TSBPMed 1 sum(KernsPMMed(:,:,WhichRepsToUse,:,:,:,:),3)};
    BARTS_Aop=BARTS_Aopx;
    
    STraj3MS(:,:,:,:,i)=STraj3M(:,:,WhichRepsToUse);
    SumKernsPMMedS(:,:,:,:,i,:,:)=sum(KernsPMMed(:,:,WhichRepsToUse,:,:,:,:),3);
    
    ksp_adj=bart(['linopScript -A ' LS_ScriptFN],ImSz16,SigToUse,BARTS_Aop.Others{:});
    ksp_adjS(:,:,:,:,SliI)=ksp_adj;
end
disp('got all ksp_adj');
%%
clear Elem0
for i=1:numel(SliIs)
    SliI=SliIs(i);
    for j=1:4
        tmp=Elems0S{SliI,2};
        Elem0{j}(:,:,:,:,:,:,:,i)=tmp{j};
    end
end
%%
% CurSRPrefix=[mainP filesep 'Several' '_'];
SliIsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((SliIs).'),1),' ','')),'[','S'),']','S');
CurSRPrefix=[mainP filesep 'Several' SliIsStr '_'];

BARTS_AopS=BARTS_Aop;
BARTS_AopS.ImSz16=FillOnesTo16([Sz 1 1 1 1 nTSMed nSlicesTogether]);
BARTS_AopS.Others{1}=perm85(SensCCS(:,:,:,:,SliIs));
BARTS_AopS.Others{2}=perm85(STraj3MS);
BARTS_AopS.Others{3}=repmat(BARTS_Aop.Others{3},[ones(1,7) nSlicesTogether]);
BARTS_AopS.Others{5}=perm85(perm73(SumKernsPMMedS));

ksp_adjFN=[CurSRPrefix 'sig_adj'];
% if(~exist([ksp_adjFN '.cfl'],'file'))
    writecfl(ksp_adjFN,perm85(perm73(ksp_adjS(:,:,:,1,SliIs))));
% end

for i=1:numel(Elem0)
    writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
%     writecfl([CurSRPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
    writecfl([CurSRPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1 nSlicesTogether]));
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
end
%
ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
ElemsLambda=[1e-5,1e-3,1e-9,1e-10];
ElemsLambda(1)=10^(-4.5);

CurSRPrefixOut=[CurSRPrefix 'Lm45_'];



writecfl([CurSRPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSRPrefix 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([CurSRPrefix 'ElementTypes'],ElementTypes.');

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
end
% ninneriterBART(2:3)=0;
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([CurSRPrefix 'ninneriter'],ninneriterBART);
disp(['saved all ' CurSRPrefix]);
%
LS_ScriptFNOnlyN='/autofs/space/daisy_002/users/Gilad/gUM/nuftAllTC_OnlyN.txt';
for i=1:4
    delete([CurSRPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSRPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue?
ContinueRun=false;
if(ContinueRun)
    for i=1:numel(Elem0)
        writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end

BARTS_AopSx=BARTS_AopS;

SensSFN=[CurSRPrefix 'Sens'];
writecfl(SensSFN,BARTS_AopS.Others{1});
BARTS_AopSx.Others{1}=SensSFN;

TrajSFN=[CurSRPrefix 'TrajS'];
writecfl(TrajSFN,BARTS_AopS.Others{2});
BARTS_AopSx.Others{2}=TrajSFN;

TSBSFN=[CurSRPrefix 'TSBS'];
writecfl(TSBSFN,BARTS_AopS.Others{3});
BARTS_AopSx.Others{3}=TSBSFN;

BARTS_AopSx.Others{4}=OneFN;

KernsSFN=[CurSRPrefix 'KernsS'];
writecfl(KernsSFN,BARTS_AopS.Others{5});
BARTS_AopSx.Others{5}=KernsSFN;


ImSz16FN=[CurSRPrefix 'ImSz16'];
writecfl(ImSz16FN,BARTS_AopS.ImSz16);
BARTS_AopSx.ImSz16=ImSz16FN;

% disp('Wrote maps');
% BARTCmdSeveral=bartCmd(['splitProx -i 10000 -I 1000 -s 600 -d 2 -g -f -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_AopS.ImSz16,OneFN,BARTS_AopS.Others{:});
BARTCmdSeveral=bartCmd(['splitProx -i 10000 -I 1000 -s 60 -d 2 -g -f -F ' CurSRPrefix ' -O ' CurSRPrefixOut ' ' LS_ScriptFNOnlyN],BARTS_AopSx.ImSz16,OneFN,BARTS_AopSx.Others{:});
% system(BARTCmd{SliI,rs});
%%
SliIs=ROrd(4:7);
% SliIs=ROrd(10:13);
SliIsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((SliIs).'),1),' ','')),'[','S'),']','S');
CurSRPrefix=[mainP filesep 'Several' SliIsStr '_'];
CurSRPrefixOut=[CurSRPrefix 'Lm45_'];
%%
CurSRPrefixOutA=CurSRPrefixOut;
%%
ErrVec=readcfl([CurSRPrefixOut 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
% figure;plot(ErrVec);

for i=1:4
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
%     Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1) '_iter3000']);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));

disp('Loaded maps');
%% Multi-shot ref
for SliI=1:nSlices
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');

THLRMultiShot_RSS(:,:,:,SliI)=THLRMultiShot_RS(:,:,:,end);
end
%%
fgmontagex(perm43(SpiRefSRgB0(:,:,2:5:end,10:13,2)));title('Separate per slice');