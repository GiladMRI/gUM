%%
PrefixBase='SKEPTIC_';
ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
LS_ScriptFN=ScriptFN_AllTS;
disp('Prepared folders');
%%
nCS=1;

CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;
M_Dim=8;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
M_Flag=2^(M_Dim-1);
%%
TimePoints_ms_GESE=[TimePoints_Med_ms TimePoints_Medb_ms];
TEs_ms=TimePoints_ms_GESE;
TEs_ms3=perm32(TEs_ms);
nEchos_GESE=numel(TEs_ms);
NTEs=Col(1:nEchos_GESE);

EchoTimes_ms_GESE_T2sActive=[TimePoints_Med_ms max(0,TimePoints_Medb_ms-TE_SE)];
EchoTimes_ms_GESE_T2qActive=[TimePoints_Med_ms*0 max(0,TE_SE-TimePoints_Medb_ms)];

EchoTimes_ms=TEs_ms.';
EchoTimes_ms3=permute32(EchoTimes_ms);
%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos_GESE,TS_Dim);
T=cat(8,perm71([ones(nTS_THLR,1);zeros(nTS_THLRb,1)]),perm71([zeros(nTS_THLR,1);ones(nTS_THLRb,1)]));
T2Sx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(Col(EchoTimes_ms_GESE_T2sActive),[TS_Dim 1]);
T2Qx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(Col(EchoTimes_ms_GESE_T2qActive),[TS_Dim 1]);
Tx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(Col(EchoTimes_ms),[TS_Dim 1]);
M = fmac(ones(nEchos_GESE,1), T,[Ch_Dim M_Dim],[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos_GESE), 1i*T,[Ch_Dim M_Dim],[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos_GESE)/1000, -1i*2*pi*Tx/1000,Ch_Dim,[CS_Dim TS_Dim]);
T2S = fmac(NTEs, -T2Sx,Ch_Dim,[CS_Dim TS_Dim]);
T2Q = fmac(NTEs, -T2Qx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%%
MaxB0=400;
M0MedianFac=20;
T2SRange=[4 200];
T2QRange=[-200 200];

WhichRepsToUse=RepsToUse;

RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

OneFN='/autofs/cluster/kawin/Gilad/One';
writecfl(OneFN,1);

ImSz16=FillOnesTo16(Sz);
ImSz16(TS_Dim)=nEchos_GESE;
ImSz16FN=[mainP filesep 'ImSz16'];
writecfl(ImSz16FN,ImSz16);
%%
clear BARTCmd
%%
m0=abs(c0);
% m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

t0q=-UpdatedT2SMap_ms_SEQ;
t0s=UpdatedT2SMap_ms_SES;

b0=B0VarAvg;

% m0=m0*0;

% Mm0=M*m0;
% expPp0 = exp(P * p0);
% expBb0 = exp(B * b0);
% expTt0 = exp(TT * (1./t0));
% 
% Est0=Mm0.*expPp0.*expBb0.*expTt0;

ElemNames={'m' 'p' 'b' 't' 't'};
nElements=numel(ElemNames);
Elem0={m0,p0,b0,t0s,t0q};
ElemL={M,P,B,T2S,T2Q};
ElemsL={T,1i*T,-1i*2*pi*Tx/1000,-T2Sx,-T2Qx};
disp('ok initials');

Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTs0 = exp(ElemL{4} * (1./Elem0{4}));
expTq0 = exp(ElemL{5} * (1./Elem0{5}));
Rec0=Mm0.*expPp0.*expBb0.*expTs0.*expTq0;
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

%

% %% Normalize
% if(grmss(Elem0{1})<eps)
%     disp('Skipping normalize, starting from 0s');
% else
%     Mm0=ElemL{1}*Elem0{1};
%     expPp0 = exp(ElemL{2} * Elem0{2});
%     expBb0 = exp(ElemL{3} * Elem0{3});
%     expTt0 = exp(ElemL{4} * (1./Elem0{4}));
%     Rec0=Mm0.*expPp0.*expBb0.*expTt0;
%     
%     AHA0=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Rec0,BARTS_Aop.Others{:});
%     
%     NFac=grmss(ksp_adj)/grmss(AHA0)
%     Elem0{1}=Elem0{1}*NFac;
%     disp(['Normalized ' num2str(NFac) ' ' num2str(grmss(ksp_adj)) ' ' num2str(grmss(AHA0))]);
% end
%
% ElemsAlphas=[1e-4, 1e+3, 1e+4, 1e+6]; reasonable for 2000
% ElemsLambda=[1e-7,1e-1,1e-9,1e-11]; reasonable for 2000
% ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
% ElemsLambda=[1e-5,1e+3,1e-9,1e-10];
% ElemsLambda=[1e-5,1e-2,1e-9,1e-10];
% ElemsLambda=[1e-5,1e-3,1e-9,1e-10];
% ElemsLambda=[1e-4,1e-4,1e-9,1e-10];
% ElemsLambda=[1e-5,1e-4,1e-9,1e-10];

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
% ninneriterBART(2:3)=0;
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
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
BartCmdFN=[mainP filesep 'splitProxRun_RS1_All'];
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
%%
WhichRepsToUse=[1 13 25];
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
clear MapsS
for SliI=1:nSlices
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    CurSRPrefix=[mainP filesep 'gB0x' 'Sli' num2str(SliI) '_R' RepsStr '_Lm45b_'];
    CurSRPrefix=[mainP filesep 'gB0x' 'Sli' num2str(SliI) '_R' RepsStr '_Lm4c_'];
%     [SpiP 'gB0x' 'Sli' num2str(SliI) '_R' RepsStr '_Lm45b_'];
    for i=1:4
        Maps{i}=readcfl([CurSRPrefix 'Elem' num2str(i-1)]);
    end
    MapsS(:,:,:,SliI)=cat(3,Maps{:});
end
disp('ok');
%%
EPTIRef=squeeze(abs(circshift(EPTIRefS(4:end-3,2:end-5,40,:),-6,2)));

SplitProxRes=squeeze(MapsS(:,:,1,ROrd).*exp(-1*42./MapsS(:,:,4,ROrd)));
Both=rot90(cat(4,EPTIRef./grms(EPTIRef,1:2),SplitProxRes./grms(SplitProxRes,1:2)));
%%
SmallVarT=200e6;
getSmallVars;
save('EPTI_SKEPTI_Perim_Comparison.mat',SmallVars{:});
%%
fgmontagex(perm43(Both(:,:,[2 6 8 11],:)));
%%
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
%%
Both=rot90(cat(4,EPTIRef./grms(EPTIRef,1:2),SplitProxRes./grms(SplitProxRes,1:2)));
XX=perm43(Both);
XX=XX./grms(XX,1:2);
%%
BB=XX*0;
for i=1:prod(gsize(XX,3:max(3,ndims(XX))))
    ForBorder=XX(:,:,i);
%     [~,B,BN]=CalcSlicesSNR(ForBorder,false,9);
    B1=getLargestComponent(ForBorder>0.5);
%     B1=getLargestComponent(B);
%     B1=getLargestComponent(~BN);
    B2=imfillholesBySlices(B1);
%     B2=~BN;
    BB(:,:,i) = bwperim(B2);
end
ZZ=BB;
ZZ(:,:,3)=0;
%%
save('EPTI_SEPTI_Perim.mat','EPTIRef','SplitProxRes');
%%
XXC=CombineDims(CombineDims(PartitionDim(squeeze(abs(XX(:,:,2,:))),3,4),[4 1]),[3 2]);XXC=XXC./grmss(XXC);
EEC=CombineDims(CombineDims(PartitionDim(squeeze(abs(XX(:,:,1,:))),3,4),[4 1]),[3 2]);EEC=EEC./grmss(EEC);

BBC1=CombineDims(CombineDims(PartitionDim(squeeze(abs(BB(:,:,1,:))),3,4),[4 1]),[3 2]);
BBC2=CombineDims(CombineDims(PartitionDim(squeeze(abs(BB(:,:,2,:))),3,4),[4 1]),[3 2]);
OverlayMaskOnImage([],XXC,BBC1,[1 0 0],0.5,BBC2,[0 1 0],0.7);removeTicks;daspect([1 1 1]);
OverlayMaskOnImage([],EEC,BBC1,[1 0 0],0.5,BBC2,[0 1 0],0.7);removeTicks;daspect([1 1 1]);

OverlayMaskOnImage([],[EEC XXC],[BBC1 BBC1],[1 0 0],0.5,[BBC2 BBC2],[0 1 0],0.7);removeTicks;daspect([1 1 1]);
%%
SliI=1;
XXC=squeeze(abs(XX(:,:,2,SliI)));XXC=XXC./grmss(XXC);
EEC=squeeze(abs(XX(:,:,1,SliI)));EEC=EEC./grmss(EEC);

BBC1=squeeze(abs(BB(:,:,1,SliI)));
BBC2=squeeze(abs(BB(:,:,2,SliI)));
OverlayMaskOnImage([],[EEC XXC],[BBC1 BBC1],[1 0 0],0.5,[BBC2 BBC2],[0 1 0],0.7);removeTicks;daspect([1 1 1]);
%%
SliI=11;
XXC=abs(CombineDims(XX(:,:,:,SliI),[3 2]));
BBC=abs(CombineDims(repmat(BB(:,:,1,SliI),[1 1 2]),[3 2]));
BBFC=abs(CombineDims(repmat(BB(:,:,2,SliI),[1 1 2]),[3 2]));
OverlayMaskOnImage([],XXC,BBC,[1 0 0],0.5,BBFC,[0 1 0],0.7);removeTicks
