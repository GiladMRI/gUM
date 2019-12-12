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
SliI=7;

disp(SliI);
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
WhichTSToUse_LLR=2:20;
Fitted0_RS_LLRS(:,:,:,SliI)=QQ.Fitted0_RS_LLR(:,:,:,2);
Fitted0_RSS(:,:,:,:,SliI)=QQ.Fitted0_RS;
UpdatedB0Map_RS_LLRS(:,:,:,SliI)=QQ.UpdatedB0Map_RS_LLR;
UpdatedT2SMap_ms_RS_LLRS(:,:,:,SliI)=QQ.UpdatedT2SMap_ms_RS_LLR;
rs=2;
[PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs),...
    Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
    FitToModel_MPBD1CSf(QQ.Rec_CompgB0_RS_MX(:,:,:,rs),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
PDBase0_RS_LLRS(:,:,SliI)=PDBase0_RS_LLR(:,:,rs);

WhichTSToUs=2:12;
clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
for rs=1:nRepsSets
    [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=...
        FitToModel_MPBD1CSf(QQ.THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end
UpdatedB0Map_RSS(:,:,:,SliI)=UpdatedB0Map_RS;
UpdatedT2SMap_ms_RSS(:,:,:,SliI)=UpdatedT2SMap_ms_RS;
disp('ok');
%% test
% WhichRepsToUse=1:39;
% nccToUse=11;
% tmp=bart(['picsS -w 1 -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16_THLR,...
%     NoNavB.*DataCCP(:,TrajPartMed,WhichRepsToUse,1:nccToUse),...
%     Sensr(:,:,:,1:nccToUse),STraj3MMed(:,:,WhichRepsToUse),TSB_THLRP,1,...
%     sum(KernsP_TSTHLR(:,:,WhichRepsToUse,:,:,:,:),3));
%%
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
RepsForKernsX=1:3;

AcqDwellTime_us=WipMemBlock.adFree{13}/1000; % 1.1

Sz=gsize(UpdatedB0Map_RSS,1:2);
nCS=1;

dTS_planned_ms=2.5;

nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);
nEchos=nTSMed;

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
%
nTraj=numel(TrajPartMed);
TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

nPointsNoNav=floor(50000/AcqDwellTime_us);
NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
NoNavB=zeros(1,nTraj);
NoNavB(1:nPointsNoNav)=1;

[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
TSBMed(nPointsMed,1)=0;
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
KernsPMMed=getKernsFromTrajM(TrajM(RepsForKernsX,TrajPartMed),Sz,TSBMed);
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
% nccToUse=15;
% WhichRepsToUse=1:3;
nccToUse=23;
WhichRepsToUse=1;
%%
SliI=11;
for SliI=1:nSlices
% c0=PDBase0_RS(:,:,1);
% b0=UpdatedB0Map_RS(:,:,1);
% t0=UpdatedT2SMap_ms_RS(:,:,1);

c0=PDBase0_RS_LLRS(:,:,SliI);
b0=UpdatedB0Map_RSS(:,:,2,SliI);
t0=UpdatedT2SMap_ms_RSS(:,:,2,SliI);

m0=abs(c0);
m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

c0=SmoothBySlices(c0,[20 20],5);
WW=SmoothBySlices(abs(c0),[20 20],5);

t0=t0*0+50;

t0=min(max(abs(t0),5),200);
b0=min(max(b0,-400),400);
m0=abs(c0);
m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

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


%%
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
Sensr=QQ.CurSens;
%
ADataIsL=ADatax.image(:,:,:,:,SliI,3,:,:,WhichRepsToUse,:,:,:,:,:,:,:,:);
ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);

Ch2D=CombineDims(ADataIsL,[3 1]);
[~,S,sccmtx] = svd(Ch2D(1:end,:),'econ');
clear Ch2D 

ncc=31;
DataC=perm43(sum(perm32(ADataIsL).*permute(sccmtx(:,1:ncc),[3 4 1 2]),3));
DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,WhichRepsToUse);
DataCCP=permute(DataPC,[1:3 8 4:7]);
clear DataC DataPC ADataIsL
disp('got sig');

%%
RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];

SigToUse=NoNavB.*DataCCP(:,TrajPartMed,WhichRepsToUse,1:nccToUse);
SigFN=[CurSRPrefix 'Sig'];
writecfl(SigFN,SigToUse);

ImSz16=FillOnesTo16(Sz);
ImSz16(TS_Dim)=nEchos;
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={Sensr(:,:,:,1:nccToUse) STraj3MMed(:,:,WhichRepsToUse) TSBPMed 1 sum(KernsPMMed(:,:,WhichRepsToUse,:,:,:,:),3)};
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

ksp_adj=bart(['linopScript -d 5 -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
disp('got ksp_adj');
%%

%% Normalize
if(grmss(Elem0{1})<eps)
    disp('Skipping normalize, starting from 0s');
else
    Mm0=ElemL{1}*Elem0{1};
    expPp0 = exp(ElemL{2} * Elem0{2});
    expBb0 = exp(ElemL{3} * Elem0{3});
    expTt0 = exp(ElemL{4} * (1./Elem0{4}));
    Rec0=Mm0.*expPp0.*expBb0.*expTt0;
    
    AHA0=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Rec0,BARTS_Aop.Others{:});
    
    NFac=grmss(ksp_adj)/grmss(AHA0)
    Elem0{1}=Elem0{1}*NFac;
    disp(['Normalized ' num2str(NFac) ' ' num2str(grmss(ksp_adj)) ' ' num2str(grmss(AHA0))]);
end
%%
% ElemsAlphas=[1e-4, 1e+3, 1e+4, 1e+6]; reasonable for 2000
% ElemsLambda=[1e-7,1e-1,1e-9,1e-11]; reasonable for 2000
ElemsAlphas=[1e-4, 1e+2, 1e+4, 1e+6];
ElemsLambda=[1e-5,1e+3,1e-9,1e-10];

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

writecfl([CurSRPrefix 'sig_adj'],ksp_adj);

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
%%
end
% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
%%
SliI=13;
for SliI=1:nSlices
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    BARTS_Aop.Others=strrep(BARTS_Aop.Others,'Sli16_',['Sli' num2str(SliI) '_']);
    SigFN=[CurSRPrefix 'Sig' RepsStr];
    QQ=bart(['splitProx -i 10000 -s 600 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
    
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
