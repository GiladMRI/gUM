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
%% Now do several THLRs
load([mainP 'modx.mat'],'modx');

RepSets={1, 1:2, 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:MaxRepsToUse};
nccForTHLR=21;
for SlicesToRead=1:nSlices
% for SlicesToRead=15
    SliI=SlicesToRead;
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    if(exist(THLR_RS_FN,'file')), continue; end
    
    clear SelfSens1 sccmtx SensCC
    CurSensFN=[SliPrefix 'SelfSens.mat'];
    load(CurSensFN,'SelfSens1');
    sccmtxFN=[SliPrefix 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccForTHLR),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');
  
    DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,1:MaxRepsToUse,:,:,:,:,:,:,:,:);
    DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
    DataCCP=CombineDims(DataCCP,[4 1]);

    DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccForTHLR),[3 4 1 2]),3));
    DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,1:MaxRepsToUse);
    DataCCP=permute(DataCCP,[1:3 8 4:7]);
    disp('read and CC sig');

    clear UpdatedB0MapTHLR_RS UpdatedT2SMap_msTHLR_RS s_valsTHLR_RS PDBase0THLR_RS
    clear THLRMultiShot_RS;
    for rs=1:numel(RepSets)
        disp(['Running THLR, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        
        tmp=bart(['picsS -w 1 -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16_THLR,NoNavB.*DataCCP(:,:,CurReps,1:nccForTHLR),...
            SensCC(:,:,:,1:nccForTHLR),STraj3M(:,:,CurReps),NoNavB.*TSB_THLRP,1,...
            sum(KernsP_TSTHLR(:,:,CurReps,:,:,:,:),3));
        THLRMultiShot_RS(:,:,:,rs)=squeeze(tmp);
        
        WhichTSToUs=2:12;
        
        [~, UpdatedB0MapTHLR_RS(:,:,rs), UpdatedT2SMap_msTHLR_RS(:,:,rs), s_valsTHLR_RS(:,:,:,rs), ~, PDBase0THLR_RS(:,:,rs)]=...
            FitToModel_MPBD1CSf(THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
        
        %     PDBase0THLRX=min(median(abs(PDBase0THLR))*10,abs(PDBase0THLR)).*exp(1i*angle(PDBase0THLR));
        %     SPDBase0THLRX=SmoothBySlices(PDBase0THLRX,[20 20],5)./SmoothBySlices(abs(PDBase0THLRX),[20 20],5);
    end
    save(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');
    disp(['saved ' THLR_RS_FN]);
end

%% LLR
% T2sCenters=1:100;
% T2Sbins=[-Inf T2sCenters-0.5 Inf];
% 
% nComponentsToUse=2;
% 
% Decays=exp(-TimePointsMed_ms./T2sCenters.');
% [Ud,Sd,Vd]=svd(Decays,'econ');
% CompsP(1,1,1,1,1,:,:)=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
% 
% clear TSCxPMed_RS CompsP_RS
% % for rs=1:nRepsSets
% rsRefI=3;
%     TSCMed_RS(:,:,:,rs)=exp(-TimePoints_Med_ms3./UpdatedT2SMap_msTHLR_RS(:,:,rsRefI)).*exp(-1i*2*pi*UpdatedB0MapTHLR_RS(:,:,rsRefI).*TimePoints_Med_ms3/1e3);
%     TSCxPMed_RS(:,:,rs,1,1,1,:)=permute(TSCMed_RS(:,:,:,rs),[1 2 7 6 5 4 3]);
%     
% % end
% disp('Prepared for LLR RepSets');
% 
% LLR_lambda=0.1;
% RhoStr=[' -u ' num2str(1e-3) ' '];
% % RhoStr=[' -u ' num2str(1e-1) ' '];
% BlkSz=4;
% ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
% Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);
% 
% TSCxPMedOnlyB0_RS=exp(1i.*angle(TSCxPMed_RS));
% %%
% clear Rec_CompgB0_C_RS
% for rs=1:nRepsSets
%     CurRep=RepSets{rs};
%     Rec_CompgB0_C_RS{rs}=bart(['picsS -m -w 1 ' RhoStr ' -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
%         Sz16CompgB0,NoNavB.*DataCCP(:,1:nTraj,CurRep,1:nccForTHLR),...
%         SensCC(:,:,:,1:nccForTHLR),STraj3M(:,:,CurRep),TSBPMed,TSCxPMedOnlyB0_RS(:,:,rs,1,1,1,:),...
%         sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP(1,1,1,1,1,:,:));
% end
% %%
% Rec_CompgB0_RS_M=cat(4,Rec_CompgB0_C_RS{:});
% Rec_CompgB0_RS_MX=perm43(squeeze(sum(Rec_CompgB0_RS_M.*gpermute(CompsP,[4 1]),6)));
% Rec_CompgB0_RS_MX=Rec_CompgB0_RS_MX.*perm43(squeeze(TSCxPMedOnlyB0_RS(:,:,rs,:,:,:,:,:,:,:)));
% 
% WhichTSToUs_LLR=2:12;
% for i=1:size(Rec_CompgB0_RS_MX,3)
% [~, UpdatedB0MapLLR_RS(:,:,i), UpdatedT2SMap_msLLR_RS(:,:,i), s_valsLLR_RS(:,:,:,i), ~, PDBase0LLR_RS(:,:,i)]=...
%             FitToModel_MPBD1CSf(squeeze(Rec_CompgB0_RS_MX(:,:,i,:)),WhichTSToUs_LLR,dT_Med_ms,TE0_ms+FirstT_Med_ms);
% end
%% end LLR
%% MLN
BaseBaseOutP=[mainP filesep 'MLN/'];
mkdir(BaseBaseOutP);
system(['chmod -R 777 ' BaseBaseOutP]);
disp(['Created ' BaseBaseOutP]);
%%
clear MSLines
SliI=16;
% for SliI=1:nSlices
    disp(SliI);
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);

CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
load(THLR_RS_FN,'UpdatedB0MapTHLR_RS');

Sensr=QQ.CurSens;
sig=QQ.CurSig;

% Sensr=imresize(Sensr,[116 116]);
Sz128=gsize(Sensr,1:2);

% B0_Hz=UpdatedB0Map_RS(:,:,WhichRSToUse);
% B0_Hz=RefB0MrS(:,:,SliI);
% B0_Hz=UpdatedB0Map_RS_LLRS(:,:,2,SliI);
B0_Hz=UpdatedB0MapTHLR_RS(:,:,3);
B0_Hz=imresize(B0_Hz,Sz128);
%
BaseOutDir=[BaseBaseOutP 'Sli' num2str(SliI) filesep];
mkdir(BaseOutDir);
system(['chmod -R 777 ' BaseOutDir]);
nCh=13;
% batchSize=8;
batchSize=16;
nTS=15;
SensCC=squeeze(Sensr(:,:,1:nCh)); % [X Y Ch]
SensMsk=grmss(SensCC,3)>0.01; % [X Y]
save([BaseOutDir 'SensCC1.mat'],'SensCC','SensMsk');

WhichRep=1;
% AcqDwellTime_us=2*1.1;
AcqDwellTime_us=3*1.1;
% TrajPartToUse=1:24000;
TrajPartToUse=1:45501;
% CurSig=     sig(1,TrajPartToUse(1:2:end),WhichRep,1:nCh)+...
%             sig(1,TrajPartToUse(2:2:end),WhichRep,1:nCh);
CurSig=     sig(1,TrajPartToUse(1:3:end),WhichRep,1:nCh)+...
            sig(1,TrajPartToUse(2:3:end),WhichRep,1:nCh)+...
            sig(1,TrajPartToUse(3:3:end),WhichRep,1:nCh);
CurSig=squeeze(CurSig);
tmp=Row(CurSig);
tmp2=[real(tmp) imag(tmp)]*600;
Data=tmp2;
Data(batchSize,end)=0;
save([BaseOutDir 'RealDataForNN.mat'],'Data');
nTrajA=size(CurSig,1);
TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TS_TimePoints3=permute(TS_TimePoints,[1 3 2]);
TSBF=GetTSCoeffsByLinear(nTrajA,nTS).';
WhichRSToUse=1;
% TSC=exp(-TS_TimePoints3./UpdatedT2SMap_ms_RS(:,:,WhichRSToUse)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,WhichRSToUse).*TS_TimePoints3/1e3);
% TSBF: [15×7162 double]
% TSC: [128×128×15 double]
save([BaseOutDir 'B0TS.mat'],'TSBF','B0_Hz');
%
% Traj=(TrajM(WhichRep,TrajPartToUse(1:2:end))+TrajM(WhichRep,TrajPartToUse(2:2:end)))/2;
Traj=(TrajM(WhichRep,TrajPartToUse(1:3:end))+TrajM(WhichRep,TrajPartToUse(2:3:end))+TrajM(WhichRep,TrajPartToUse(3:3:end)))/3;
clear Trajm2
Trajm2(1,:)=real(Traj);
Trajm2(2,:)=imag(Traj);
[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
save([BaseOutDir 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');

TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TS_TimePoints=TE0_ms+linspace(0,TimePoints_ms(end),nTS);
TSstr=strrep(num2str(TS_TimePoints,'%3.5f,'),' ','');
TSstr=TSstr(1:end-1);
% TSstr=['TimePoints_ms ' TSstr(1:end-1)];

TS_TimePointsForRec=linspace(0,TimePoints_ms(end),8);
TS_TimePointsForRec_Str=strrep(num2str(TS_TimePointsForRec,'%3.5f,'),' ','');
TS_TimePointsForRec_Str=TS_TimePointsForRec_Str(1:end-1);

St=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForS10.txt');
St.LabelsW=Sz128(1);
St.LabelsH=Sz128(2);
St.aDataW=Sz128(1);
St.aDataH=Sz128(2);
St.nTraj=size(CurSig,1);
St.nccInData=nCh;
St.RealDatancc=nCh;
St.DataH=size(CurSig,1)*nCh*2;
St.TimePoints_ms=TSstr;
St.batch_size=batchSize;
St.nToLoad=10000;
St.TimePointsForRec_ms=TS_TimePointsForRec_Str;
St.SessionNameBase=[mainP filesep 'MLN/S' num2str(SliI)];
St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
St.BaseTSDataP=BaseOutDir;
St.BaseNUFTDataP=BaseOutDir;
% if(strcmp(HostName,'tiger'))
% end
% disp(TSstr)
% ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForMES' num2str(SliI) '.txt'];
ParamsOutFn=[BaseBaseOutP 'ParamsForMES' num2str(SliI) '.txt'];
Txt=gStruct2txt(St,ParamsOutFn);
disp('ok')
MSLines{SliI}=['/usr/pubsw/bin/python3 /autofs/cluster/kawin/Gilad/TF/srezN/srez_main1.py ' ParamsOutFn];
% end
disp('Prepared for MLN 1-shot');
%%
MSFN=[BaseBaseOutP 'MLNRun'];
fid=fopen(MSFN,'w+');
% fprintf(fid,'#!/bin/bash\r\n');
fprintf(fid,'# MLN run\r\n');
fprintf(fid,'echo "Starting!"\r\n');
for i=1:numel(MSLines)
    fprintf(fid,'%s\n',MSLines{i});
end
fclose(fid);
[status,msg,msgID] = fileattrib(MSFN,'+x');
disp(['Prepared script for MLN 1-shot runs: ' MSFN]);
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
RepsStrA=RepsStr;
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
SliI=7;
for SliI=1:nSlices
    %%
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    
    clear SelfSens1 sccmtx SensCC
    CurSensFN=[SliPrefix 'SelfSens.mat'];
    load(CurSensFN,'SelfSens1');
    sccmtxFN=[SliPrefix 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccToUse),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');
  
    DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,1:MaxRepsToUse,:,:,:,:,:,:,:,:);
    DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
    DataCCP=CombineDims(DataCCP,[4 1]);

    DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccToUse),[3 4 1 2]),3));
    DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,1:MaxRepsToUse);
    DataCCP=permute(DataCCP,[1:3 8 4:7]);
    disp('read and CC sig');
    %%    
    CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
    THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
    
    load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');
%%
for rs=1:numel(RepSets)
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

        CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];

c0=PDBase0THLR_RS(:,:,rs);
b0=UpdatedB0MapTHLR_RS(:,:,rs);
t0=UpdatedT2SMap_msTHLR_RS(:,:,rs);

%%
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
%%


m0=abs(c0);
m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

c0=SmoothBySlices(c0,[20 20],5);
WW=SmoothBySlices(abs(c0),[20 20],5);

t0=t0*0+50;
t0=min(max(abs(t0),T2SRange(1)),T2SRange(2));
b0=min(max(b0,-MaxB0),MaxB0);
m0=abs(c0);
m0=min(m0,median(m0(:))*M0MedianFac);
p0=angle(c0);
c0=m0.*exp(1i*p0);

b0=SB0x;

m0=m0*0;

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Est0=Mm0.*expPp0.*expBb0.*expTt0;

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
Elem0={m0,p0,b0,t0};
ElemL={M,P,B,TT};
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
disp('ok initials');
%%
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
ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
% ElemsLambda=[1e-5,1e+3,1e-9,1e-10];
ElemsLambda=[1e-5,1e-2,1e-9,1e-10];
ElemsLambda=[1e-5,1e-3,1e-9,1e-10];

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
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([CurSRPrefix 'ninneriter'],ninneriterBART);
disp(['saved all ' CurSRPrefix]);

for i=1:4
    delete([CurSRPrefix 'Elem' num2str(i-1) '.hdr']);
    delete([CurSRPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue?
% for i=1:numel(Elem0)
%     writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
% end
% disp('Wrote maps');
BARTCmd{SliI,rs}=bartCmd(['splitProx -i 10000 -s 600 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,OneFN,BARTS_Aop.Others{:});
%%
end
end
%%
for SliI=1:nSlices
    for rs=1:numel(RepSets)
        disp(['Preparing for splitProx, Slice ' num2str(SliI) ' Reps set ' num2str(rs) ': ' datestr(now)]);
        CurReps=RepSets{rs};
        WhichRepsToUse=CurReps;
        RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

        CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];

        ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
        % ElemsLambda=[1e-5,1e+3,1e-9,1e-10];
        ElemsLambda=[1e-5,1e-2,1e-9,1e-10];
        ElemsAlphas=[1e-5, 1e-0, 1e+4, 1e+6];
        ElemsLambda=[1e-5,1e-4,1e-9,1e-10];

        writecfl([CurSRPrefix 'ElemsAlpha'],ElemsAlphas.');
        writecfl([CurSRPrefix 'ElemsLambda'],ElemsLambda.');
    end
end
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

BartCmdFN=[mainP filesep 'splitProxRun_RS_All'];
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
SliI=13;
for SliI=1:nSlices
    
    rs=1;
    CurReps=RepSets{rs};
        
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    BARTS_Aop.Others=strrep(BARTS_Aop.Others,'Sli16_',['Sli' num2str(SliI) '_']);
    QQ=bart(['splitProx -i 10000 -s 60 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,1,BARTS_Aop.Others{:});
    
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

rs=1;
CurReps=RepSets{rs};
WhichRepsToUse=CurReps;
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
        
CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];


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

CurSRPrefixA=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStrA '_'];
THLR_RS_FN=[CurSRPrefixA 'THLRres_RS.mat'];
load(THLR_RS_FN,'THLRMultiShot_RS','UpdatedB0MapTHLR_RS','UpdatedT2SMap_msTHLR_RS','s_valsTHLR_RS','PDBase0THLR_RS');

Recs=cat(4,THLRMultiShot_RS(:,:,2:4:end,rs),RecMX(:,:,2:6:end));
fgmontagex(perm43(Recs));caxis(caxis/1.5);
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
