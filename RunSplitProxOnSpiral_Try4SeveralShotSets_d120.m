mainPRef='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1/';
%%
% CurPrefix=[HostName '_SKEPTIC_'];
PrefixBase='SKEPTIC_';
% CurPrefix=[HostName '_' PrefixBase];
% ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
% ToBARTP=[mainP filesep];
% BaseSP=['/autofs/space/daisy_002/users/Gilad/' CurPrefix];
% BaseFP=['/tmp/' CurPrefix];
LS_ScriptFN=ScriptFN_AllTS;

disp('Prepared folders');
%%
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
MaxRepsToUse=36;
RepsForKernsX=1:MaxRepsToUse;
AcqDwellTime_us=WipMemBlock.adFree{13}/1000; % 1.1
% Sz=gsize(UpdatedB0Map_RSS,1:2);
nCS=1;

dTS_planned_ms=2.5;

nTSMed=ceil((nTraj+1)*AcqDwellTime_us/1000/dTS_planned_ms);
nEchos=nTSMed;

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nTraj)/1000,nTSMed);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
%
TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

% nPointsNoNav=floor(50000/AcqDwellTime_us);
% NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
% NoNavB=zeros(1,nTraj);
% NoNavB(1:nPointsNoNav)=1;

[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nTraj,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*TotalAcqTime_ms;
TimePoints_Med_ms=TimePointsR_Med*TotalAcqTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
TSBMed(nTraj,1)=0;
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
KernsPMMed=getKernsFromTrajM(TrajM(RepsForKernsX,:),Sz,TSBMed);
disp('ok KernsPMMed');
%%
RepsStrA='A1_2_3A';
%% Now do several THLRs
load([mainP 'modx.mat'],'modx');

% RepSets={1, 1:2, 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:MaxRepsToUse};
RepSets={1, [1 13 25], 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:MaxRepsToUse};
nccForTHLR=21;
for SlicesToRead=1:nSlices
% for SlicesToRead=15
    SliI=SlicesToRead;
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    if(exist(THLR_RS_FN,'file')), disp('exist, Skipping'); continue; end
    
    SliPrefixRef=[mainPRef 'Sli' num2str(SliI) '_'];
    CurSRPrefixRef=[mainPRef 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    CurSRPrefixARef=[mainPRef 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    
%     clear SelfSens1 sccmtx SensCC
%     CurSensFN=[SliPrefixRef 'SelfSens.mat'];
%     load(CurSensFN,'SelfSens1');
%     sccmtxFN=[SliPrefixRef 'sccmtx.mat'];
%     load(sccmtxFN,'sccmtx');
%     
%     SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccForTHLR),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
%     SensCC=permute43(SensCC);
%     disp('ok SensCC');
%   
%     DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,1:MaxRepsToUse,:,:,:,:,:,:,:,:);
%     DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
%     DataCCP=CombineDims(DataCCP,[4 1]);
% 
%     DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccForTHLR),[3 4 1 2]),3));
%     DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,1:MaxRepsToUse);
%     DataCCP=permute(DataCCP,[1:3 8 4:7]);
%     disp('read and CC sig');
    
    clear UpdatedB0MapTHLR_RS UpdatedT2SMap_msTHLR_RS s_valsTHLR_RS PDBase0THLR_RS
    clear THLRMultiShot_RS;

    THLR_RS_FNRef=[CurSRPrefixARef 'THLRres_RS.mat'];
    load(THLR_RS_FNRef);

    save(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');
    disp(['saved ' THLR_RS_FN]);
end
%%
for SliI=1:nSlices
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    load(THLR_RS_FN);
    disp(SliI);
    THLRMultiShot_RSS(:,:,:,SliI)=THLRMultiShot_RS(:,:,:,end);
end
%%



%% ScriptFN_AllTS
% FMAC 0 16      Multiply with sens, sum over maps
% NUFFT 1 -1 -1 7 0 0 0 0
% FMAC 2 64       Multiply with TSB and sum over TS
% NORMAL 123
% f 0
% dblsz 3
% fft 3
% fmac 4 0
% ifft 3
% halfsz 3
% a 0
%%
TEs_ms=TE0_ms+TimePointsMed_ms.';

% TEs_ms=(FirstTE_ms+(0:(nEchos-1))*ES_ms).';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms=TEs_ms.';
EchoTimes_ms3=permute32(EchoTimes_ms);
%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%%
% RepsStrA=RepsStr;
%%
% nccToUse=15;
WhichRepsToUse=1:3;
nccToUse=23;
% WhichRepsToUse=1;

MaxB0=400;
M0MedianFac=20;
T2SRange=[4 200];

RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

OneFN='/autofs/cluster/kawin/Gilad/One';
writecfl(OneFN,1);

ImSz16=FillOnesTo16(Sz);
ImSz16(TS_Dim)=nEchos;
ImSz16FN=[mainP filesep 'ImSz16'];
writecfl(ImSz16FN,ImSz16);
%%
clear BARTCmd
%%
for SliI=1:nSlices
% for SliI=[2 14]
    %
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    
    SliPrefixRef=[mainPRef 'Sli' num2str(SliI) '_'];
    CurSRPrefixRef=[mainPRef 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    CurSRPrefixARef=[mainPRef 'Sli' num2str(SliI) '_R' RepsStrA '_'];

    clear SelfSens1 sccmtx SensCC
    CurSensFN=[SliPrefixRef 'SelfSens.mat'];
    load(CurSensFN,'SelfSens1');
    sccmtxFN=[SliPrefixRef 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccToUse),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');
    SensCCS(:,:,SliI,:)=SensCC;
    
    sccmtxS(:,:,SliI)=sccmtx;
end
%%
RepSets={1:3 4:6 7:9 10:12 13:15 16:18 19:21 22:24 25:27 28:30 31:33 34:36};
nRepsSets=numel(RepSets);
%%
SliI=7;
for SliI=1:nSlices
% for SliI=[2 14]
    %
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    
    SensCC=SensCCS(:,:,SliI,:);
    sccmtx=sccmtxS(:,:,SliI);
  
    DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,1:MaxRepsToUse,:,:,:,:,:,:,:,:);
    DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
    DataCCP=CombineDims(DataCCP,[4 1]);

    DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccToUse),[3 4 1 2]),3));
    DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,1:MaxRepsToUse);
    DataCCP=permute(DataCCP,[1:3 8 4:7]);
    disp('read and CC sig');
    %
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');
%
for rs=1:nRepsSets
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

%         CurSRPrefix=[mainP filesep 'gB0Sli' num2str(SliI) '_R' RepsStr '_'];
        CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];

% c0=PDBase0THLR_RS(:,:,rs);
% b0=UpdatedB0MapTHLR_RS(:,:,rs);
% t0=UpdatedT2SMap_msTHLR_RS(:,:,rs);

% c0=PDBase0THLR_RS(:,:,rs);
c0=PDBase0THLR_RS(:,:,end);close
b0=UpdatedB0MapTHLR_RS(:,:,end);
t0=c0*0+50;

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

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
Elem0={m0,p0,b0,t0};
ElemL={M,P,B,TT};
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
disp('ok initials');
%
SigToUse=DataCCP(:,:,WhichRepsToUse,:);

% SigToUse=NoNavB.*DataCCP(:,TrajPartMed,WhichRepsToUse,1:nccToUse);
SigFN=[CurSRPrefix 'Sig'];
writecfl(SigFN,SigToUse);

BARTS_Aopx.ImSz16=ImSz16FN;
BARTS_Aopx.Others={SensCC STraj3M(:,:,WhichRepsToUse) TSBPMed 1 sum(KernsPMMed(:,:,WhichRepsToUse,:,:,:,:),3)};
% BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
BARTS_Aop=BARTS_Aopx;
% Suffixes={'Sens' ['STraj' RepsStr] 'TSBPMed' 'One' ['sumKerns' RepsStr]};
Suffixes={'Sens' 'STraj' 'TSBPMed' 'One' 'sumKerns'};
OtherFNs=cell(1,5); 
OtherFNs{4}=[ToBARTP Suffixes{4}];
for i=[1 3]
    OtherFNs{i}=[CurSPrefix Suffixes{i}];
end
for i=[2 5]
    OtherFNs{i}=[CurSRPrefix Suffixes{i}];
end
for i=1:5
    writecfl(OtherFNs{i},BARTS_Aopx.Others{i});
    BARTS_Aop.Others{i}=OtherFNs{i};
end

ksp_adjFN=[CurSRPrefix 'sig_adj'];
if(~exist([ksp_adjFN '.cfl'],'file'))
    ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
    disp('got ksp_adj');
    writecfl(ksp_adjFN,ksp_adj);
end

ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
% ElemsLambda=[1e-5,1e-4,1e-9,1e-10];
% ElemsLambda(1)=10^(-4.5);
ElemsLambda=[1e-4,1e-4,1e-9,1e-10];
        
for i=1:numel(Elem0)
    writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end
for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
    writecfl([CurSRPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[OS_Sz 1 1 1 1 1]));
end

writecfl([CurSRPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSRPrefix 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([CurSRPrefix 'ElementTypes'],ElementTypes.');

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
end
writecfl([CurSRPrefix 'ninneriter'],ninneriterBART);
disp(['saved all ' CurSRPrefix]);

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

% CurSRPrefixOut=[CurSRPrefix 'Lm5_'];
CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];

% disp('Wrote maps');
BARTCmd{SliI,rs}=bartCmd(['splitProx -i 10000 -I 1000 -s 600 -d 2 -g -f -F ' CurSRPrefix ' -O ' CurSRPrefixOut ' ' LS_ScriptFN],BARTS_Aop.ImSz16,OneFN,BARTS_Aop.Others{:});
% system(BARTCmd{SliI,rs});
end
end
disp('Prepared BARTCmd');
%%
% for SliI=[2 14] % 1:nSlices
for SliI=1:nSlices
    for rs=2 % 1:numel(RepSets)
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

%         CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
        CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];
        
%         CurSRPrefixOut=[CurSRPrefix 'Lm45b_'];
        CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];

        ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
        % ElemsLambda=[1e-5,1e+3,1e-9,1e-10];
%         ElemsLambda=[1e-5,1e-2,1e-9,1e-10];
%         ElemsAlphas=[1e-5, 1e-0, 1e+4, 1e+6];
%         ElemsLambda=[1e-5,1e-4,1e-9,1e-10];
        ElemsLambda=[1e-4,1e-4,1e-9,1e-10];
%         ElemsLambda=[1e-5,1e-4,1e-9,1e-10];
%         ElemsLambda(1)=10^(-4.5);

        writecfl([CurSRPrefix 'ElemsAlpha'],ElemsAlphas.');
        writecfl([CurSRPrefix 'ElemsLambda'],ElemsLambda.');
        disp(['ok ' CurSRPrefixOut]);
        
        F1=strfind(BARTCmd{SliI,rs},'-O ');
        F2=strfind(BARTCmd{SliI,rs},'/autofs/space/daisy_002/users/Gilad/gUM/nuftAllTSC_N.txt');
        BARTCmdx{SliI,rs}=[BARTCmd{SliI,rs}(1:F1+2) CurSRPrefixOut BARTCmd{SliI,rs}(F2-1:end)];
        BARTCmdx{SliI,rs}=strrep(BARTCmdx{SliI,rs},'10000','20000');
    end
end
clipboard('copy',BARTCmdx{7,rs});
%%
clear RanRec
for SliI=1:nSlices
    for rs=1:numel(RepSets)
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
        CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
        RanRec(SliI,rs)=exist([CurSRPrefix 'Elem0.cfl'],'file')>0;
    end
end
disp('got RanRec');
%%
% BARTCmds=BARTCmd(:);
for SliI=1:nSlices
% BARTCmds=BARTCmd(SliI,:).';
BARTCmds=BARTCmd(SliI,~RanRec(SliI,:)).';
% BARTCmds=BARTCmd(2,~RanRec(2,:)).';
BARTCmds=BARTCmds(cellNumel(BARTCmds)>0);
% BARTCmd=strrep(BARTCmd.','-s 60 ','-s 600 ');
%
% BartCmdFN=[mainP filesep 'splitProxRun_RS'];
BartCmdFN=[mainP filesep 'splitProxRun_RS' num2str(SliI)];
fid=fopen(BartCmdFN,'w+');
fprintf(fid,'# splitProx run\r\n');
fprintf(fid,'echo "Starting!"\r\n');
for i=1:numel(BARTCmds)
    fprintf(fid,'%s\n',BARTCmds{i});
end
fclose(fid);
[status,msg,msgID] = fileattrib(BartCmdFN,'+x');
disp(['Prepared script ' BartCmdFN]);
end
% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
%%
BARTCmds=BARTCmd.';
BARTCmds=BARTCmds(:);

BARTCmds=BARTCmd(:,1);
BARTCmds=BARTCmd(14,2:end);
BARTCmds=BARTCmd([1:13 15:16],:);
BARTCmds=BARTCmds(:);

BARTCmds=strrep(BARTCmds,'-g -f -F','-g -F');
BartCmdFN=[mainP filesep 'splitProxRun_RS_All'];
% BartCmdFN=[mainP filesep 'splitProxRun_RS1_All'];
delete(BartCmdFN);
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
SliI=2;
for SliI=1:nSlices
    
    rs=1;
    CurReps=RepSets{rs};
        
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    BARTS_Aop.Others=strrep(BARTS_Aop.Others,'Sli16_',['Sli' num2str(SliI) '_']);
%     QQ=bart(['splitProx -i 10000 -s 60 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
    
%     ErrVec=readcfl([CurSRPrefix 'ErrVec']);
%     ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
    
    for i=1:4
        Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1)]);
    end

    MapsOutFN=[CurSRPrefix 'Maps.mat'];
    save(MapsOutFN,'Maps');
    disp(['Saved ' MapsOutFN]);
end
%%
% ErrVec=readcfl([CurSRPrefix 'ErrVec']);
% ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
% figure;plot(ErrVec);setYaxis([0 max(ErrVec)*1.1]);
%%
SliI=2;

rs=2;
CurReps=RepSets{rs};
WhichRepsToUse=CurReps;
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
        
CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];

CurSRPrefixOut=[CurSRPrefix 'Lm5_'];
%%
SliI=2;
SliI=3;
% SliI=14;

% CurSRPrefix=[mainP filesep 'gB0Sli' num2str(SliI) '_R' RepsStr '_'];
CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];

% CurSRPrefixOut=[CurSRPrefix 'Lm5_'];
CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];
% CurSRPrefixOut=[CurSRPrefix 'Lm5b_'];
% CurSRPrefixOut=[CurSRPrefix 'Lm45b_'];

ErrVec=readcfl([CurSRPrefixOut 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);

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
%%
for SliI=1:nSlices
CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];
% CurSRPrefixOut=[CurSRPrefix 'Lm45b_'];
CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];
disp(CurSRPrefixOut);
for i=1:4
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter3000']);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMXS(:,:,:,SliI)=RecMX;

MapsS(:,:,:,SliI)=cat(3,Maps{:});
end
disp('Loaded maps');
%%
%%
for SliI=1:nSlices
for rs=1:nRepsSets
    CurReps=RepSets{rs};
WhichRepsToUse=CurReps;
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];
% CurSRPrefixOut=[CurSRPrefix 'Lm45b_'];
CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];
disp(CurSRPrefixOut);
try
for i=1:4
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter3000']);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));
RecMXSR(:,:,:,SliI,rs)=RecMX;

MapsS(:,:,:,SliI,rs)=cat(3,Maps{:});
catch
end
end
end
disp('Loaded maps');
%%
QQ=squeeze(RecMXSR(:,:,13,14,:));
mQQ=mean(abs(QQ),3);
sQQ=std(abs(QQ),[],3);
tSNR=mQQ./sQQ;
fgmontagex(tSNR,[0 100]);colorbar;title('tSNR SKEPTIC 3-shot');

figure;histogram(tSNR(mQQ>0.00013),0:5:200);title('tSNR SKEPTIC 3-shot')
%%
aRecMXSR=single(abs(RecMXSR));
save([mainP filesep 'aRecMXSR.mat'],'aRecMXSR');
%%
figure;
for SliI=1:nSlices
    QQ=squeeze(RecMXSR(:,:,13,SliI,:));
    RSran(SliI)=sum(grmss(QQ,1:2)>0);
    QQ=QQ(:,:,grmss(QQ,1:2)>0);
    mQQ=mean(abs(QQ),3);
    sQQ=std(abs(QQ),[],3);
    mRecMXSR(:,:,SliI)=mQQ;
    tSNRS(:,:,SliI)=mQQ./sQQ;
    tmp=tSNRS(:,:,SliI);
    subplot(4,4,ROrd(SliI));
    HMsk(:,:,SliI)=mQQ>0.00013;
    h=histogram(tmp(HMsk(:,:,SliI)),0:5:200);title('tSNR SKEPTIC 3-shot');hold on;
    counts{ROrd(SliI)} = h.Values;
end
fgmontagex(tSNRS(:,:,ROrd),[0 100]);colorbar;title('tSNR SKEPTIC 3-shot');

figure;histogram(tSNR(mQQ>0.00026),0:5:200);title('tSNR SKEPTIC 3-shot')

figure;histogram(tSNR(mQQ>0.00013),0:5:200);title('tSNR SKEPTIC 3-shot')
save('SKEPTIC_Counts.mat','counts');
%%
[~,MskB,MskBN]=CalcSlicesSNR(mRecMXSR,false,3);
MskB=imfillholesBySlices(MskB);
MskBN=imfillholesBySlices(~MskBN)>0.5;
for i=1:nSlices
    MskBN(:,:,i)=getLargestComponent(MskBN(:,:,i));
end

fgmontagex(mRecMXSR);
fgmontagex(MskB);
fgmontagex(MskBN);

hS=histogram(tSNRS(MskBN),0:5:200);title('tSNR SKEPTIC 3-shot');hold on;

figure;histogram(tSNRS(MskBN),0:5:200);title('tSNR SKEPTIC 3-shot')

%%
SmallVarT=50e6;
getSmallVars;
save([mainP filesep 'SmallVars.mat'],SmallVars{:});
%%
clear MapsI
for j=1:9
for i=1:4
    
%     Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1)]);
    MapsI{i}(:,:,j)=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter' num2str(j*1000)]);
    end
end

%%
Rec214=cat(4,RecMX2,RecMX14);
%%
AllRecs=cat(4,EPTIRefS(4:end-3,7:end,3:25:end,[4 11]),Rec214(:,:,2:5:end,:));
AllRecs=AllRecs./grms(AllRecs,1:3);
%%
AllRecs=cat(4,EPTIRefS(3:end-2,5:end,3:25:end,14),RecMX(:,:,2:5:end,:));
AllRecs=gflip(AllRecs,1)./grms(AllRecs,1:3);
%%
CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');

Recs=cat(4,THLRMultiShot_RS(:,:,2:4:end,rs),RecMX(:,:,2:6:end));
% fgmontagex(perm43(Recs));caxis(caxis/1.5);
fgmontagex(perm43(padarray(Recs,[3 3],'Both')));caxis(caxis/1.5);title(['splitProx ' num2str(rs) '-shot']);
%%
% SliI=15;
SliI=7;

CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');
CollectedRS=[];
clear RecMX_RS
for rs=1:numel(RepSets)
    % rs=4;
    CurReps=RepSets{rs};
    WhichRepsToUse=CurReps;
    RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    
    if(~exist([CurSRPrefix 'Elem0.cfl'],'file')), continue, end
    
    ErrVec=readcfl([CurSRPrefix 'ErrVec']);
    ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
    
    for i=1:4
        Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1)]);
    end
    
    MmM=ElemL{1}*Maps{1};
    expPpM = exp(ElemL{2} * Maps{2});
    expBbM = exp(ElemL{3} * Maps{3});
    expTtM = exp(ElemL{4} * (1./Maps{4}));
    RecM=MmM.*expPpM.*expBbM.*expTtM;
    RecMX=squeeze(sum(RecM,CS_Dim));
    
    disp('Loaded maps');
    
    RecMX_RS(:,:,:,rs)=RecMX;
    CollectedRS=[CollectedRS rs];
end
disp('Collected RecMX_RS');
fgmontagex(RecMX_RS(:,:,2:4:end,CollectedRS),[0 1e-3]);title(num2str(cellNumel(RepSets(CollectedRS))));
%%
THLR_FN=[CurSRPrefix 'THLRres.mat'];
load(THLR_FN,'THLRMultiShot');

Recs=cat(4,THLRMultiShot(:,:,2:4:end),RecMX(:,:,2:6:end));
fgmontagex(perm43(Recs));caxis(caxis/1.5);

%%
% MapsOutFN=[mainP filesep PrefixBase 'Maps_S' num2str(SliI) '.mat'];
MapsOutFN=[CurSRPrefix 'Maps.mat'];
save(MapsOutFN,'Maps');
disp(['Saved ' MapsOutFN]);

% end
%%
[~,Pref,~]=fileparts(CurSRPrefix);
ElemRanges={[0 11e-4],[-pi pi],[-400 400],[0 100]};
mThresh=60e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); xlabel([num2str(numel(ErrVec)/sum(ninneriterBART)) ' : ' num2str(ErrVec(end),'%.7g')]); end
    if(i==3), xlabel(num2str(ElemsAlphas,' %.9g,')); ylabel(Pref,'Interpreter','None'); end
    if(i==4), xlabel(num2str(ElemsLambda,' %.9g,')); end
end
%%
THLR_FN=[CurSRPrefix 'THLRres.mat'];
load(THLR_FN,'THLRMultiShot');

% fgmontagex(THLRMultiShot(:,:,2:4:end),[0 1e-3]);title('THLR3Shot');
% 
% fgmontagex(RecMX(:,:,2:6:end),[0 1e-3]);title('splitProx');

Recs=cat(4,THLRMultiShot(:,:,2:4:end),RecMX(:,:,2:6:end));
fgmontagex(perm43(Recs));caxis(caxis/1.5);
%%
Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTt0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
Rec0X=squeeze(sum(Rec0,CS_Dim));
disp('WS maps');
%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Elem0{i},ElemRanges{i});removeTicks;
    if(i==1), title('Start point'); end
end
%%
c1=PDBase0_RS(:,:,1);
b1=UpdatedB0Map_RS(:,:,1);
t1=UpdatedT2SMap_ms_RS(:,:,1);

t1=min(max(abs(t1),5),200);
b1=min(max(b1,-400),400);
m1=abs(c1);
m1=min(m1,median(m1(:))*20);
p1=angle(c1);
c1=m1.*exp(1i*p1);

Elem1={m1 p1 b1 t1};
Mm1=ElemL{1}*Elem1{1};
expPp1 = exp(ElemL{2} * Elem1{2});
expBb1 = exp(ElemL{3} * Elem1{3});
expTt1 = exp(ElemL{4} * (1./Elem1{4}));
Rec1=Mm1.*expPp1.*expBb1.*expTt1;
%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Elem1{i},ElemRanges{i});removeTicks;
    if(i==1), title('Multi-shot'); end
end
%% Collect
WhichRepsToUse=1:3;
WhichRepsToUse=1;
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
clear MapsS
for SliI=1:nSlices
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    for i=1:4
        Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1)]);
    end
    MapsS(:,:,:,SliI)=cat(3,Maps{:});
end
disp('ok');
%%
SplitProxRes=squeeze(MapsS(:,:,1,ROrd).*exp(-1*42./MapsS(:,:,4,ROrd)));
%%
ForBorder=abs(circshift(EPTIRefS(3:end-2,1:end-4,40,11),-6,2));
ForBorder=abs(circshift(EPTIRefS(4:end-3,2:end-5,40,11),-6,2));
ForBorder=abs(circshift(EPTIRefS(4:end-3,2:end-5,40,4),-6,2));
ForBorder=abs(RecMX(:,:,11,:));

[~,B,BN]=CalcSlicesSNR(ForBorder,true,9);
B1=getLargestComponent(B);
B2=imfillholesBySlices(B1);
BW3 = bwperim(B2);

% BWnobord = imclearborder(BWdfill,4);

% BW1 = edge(B2,'sobel');
% BW2 = edge(B2,'canny');

% figure;
% imshowpair(BW1,BW3,'montage')
%%
XX=cat(3,abs(circshift(EPTIRefS(4:end-3,2:end-5,40,4),-6,2)),abs(RecMX2(:,:,11,:)));
% XX=cat(3,abs(circshift(EPTIRefS(4:end-3,2:end-5,40,11),-6,2)),abs(RecMX(:,:,11,:)));
% XX=cat(3,abs(circshift(EPTIRefS(4:end-3,2:end-5,40,6),-6,2)),abs(RecMX(:,:,11,:)));
XX=XX./grms(XX,1:2);

BB=XX*0;
for i=1:prod(gsize(XX,3:max(3,ndims(XX))))
    ForBorder=XX(:,:,i);
    [~,B,BN]=CalcSlicesSNR(ForBorder,false,9);
    B1=getLargestComponent(B);
%     B1=getLargestComponent(~BN);
    B2=imfillholesBySlices(B1);
    BB(:,:,i) = bwperim(B2);
end
ZZ=BB;
ZZ(:,:,3)=0;

%%
XXC=abs(CombineDims(XX,[3 2]));
BBC=abs(CombineDims(repmat(BB(:,:,1),[1 1 2]),[3 2]));
BBFC=abs(CombineDims(repmat(BB(:,:,2),[1 1 2]),[3 2]));
OverlayMaskOnImage([],XXC,BBC,[1 0 0],0.5,BBFC,[0 1 0],0.7);removeTicks
%%
XX=cat(3,abs(circshift(EPTIRefS(4:end-3,2:end-5,40,:),-6,2)),abs(RecMXS(:,:,11,ROrd)));
XX=XX./grms(XX,1:2);

BB=XX*0;
for i=1:prod(gsize(XX,3:max(3,ndims(XX))))
    ForBorder=XX(:,:,i);
%     [~,B,BN]=CalcSlicesSNR(ForBorder,false,9);
%     B1=getLargestComponent(B);
%     B1=getLargestComponent(~BN);
    B1=getLargestComponent(ForBorder>0.8);
    B2=imfillholesBySlices(B1);
    BB(:,:,i) = bwperim(B2);
end
%%
XXC=CombineDims(CombineDims(PartitionDim(squeeze(abs(RecMXS(:,:,11,ROrd))),3,4),[4 1]),[3 2]);XXC=XXC./grmss(XXC);
EEC=CombineDims(CombineDims(PartitionDim(squeeze(abs(circshift(EPTIRefS(4:end-3,2:end-5,40,:),-6,2))),3,4),[4 1]),[3 2]);EEC=EEC./grmss(EEC);
BBC1=CombineDims(CombineDims(PartitionDim(squeeze(abs(BB(:,:,1,:))),3,4),[4 1]),[3 2]);
BBC2=CombineDims(CombineDims(PartitionDim(squeeze(abs(BB(:,:,2,:))),3,4),[4 1]),[3 2]);
OverlayMaskOnImage([],XXC,BBC1,[1 0 0],0.5,BBC2,[0 1 0],0.7);removeTicks;daspect([1 1 1]);
OverlayMaskOnImage([],EEC,BBC1,[1 0 0],0.5,BBC2,[0 1 0],0.7);removeTicks;daspect([1 1 1]);

OverlayMaskOnImage([],[EEC XXC],[BBC1 BBC1],[1 0 0],0.5,[BBC2 BBC2],[0 1 0],0.7);removeTicks;daspect([1 1 1]);
%%
SliI=11;
XXC=abs(CombineDims(XX(:,:,:,SliI),[3 2]));
BBC=abs(CombineDims(repmat(BB(:,:,1,SliI),[1 1 2]),[3 2]));
BBFC=abs(CombineDims(repmat(BB(:,:,2,SliI),[1 1 2]),[3 2]));
OverlayMaskOnImage([],XXC,BBC,[1 0 0],0.5,BBFC,[0 1 0],0.7);removeTicks
