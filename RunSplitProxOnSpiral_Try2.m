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
RepsForKernsX=1:3;
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
SliI=7;
for SliI=1:nSlices
    %%
    SliPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    disp(SliPrefix);
    CurSPrefix=[mainP filesep 'Sli' num2str(SliI) '_'];
    CurSRPrefix=[mainP filesep 'Sli' num2str(SliI) '_R' RepsStr '_'];
    
    
    THLR_FN=[CurSRPrefix 'THLRres.mat'];
    
    clear SelfSens1 sccmtx SensCC
    CurSensFN=[SliPrefix 'SelfSens.mat'];
    load(CurSensFN,'SelfSens1');
    sccmtxFN=[SliPrefix 'sccmtx.mat'];
    load(sccmtxFN,'sccmtx');
    
    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:nccToUse),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');
  
    DataCCP=ADatax.image(:,:,:,:,SliI,3,:,:,WhichRepsToUse,:,:,:,:,:,:,:,:);
    DataCCP=permute(DataCCP,[1 2 9 11 5 3:4 6:8 10]);
    DataCCP=CombineDims(DataCCP,[4 1]);

    DataCCP=perm43(sum(perm32(DataCCP).*permute(sccmtx(:,1:nccToUse),[3 4 1 2]),3));
    DataCCP=permute(DataCCP(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx(:,:,WhichRepsToUse);
    DataCCP=permute(DataCCP,[1:3 8 4:7]);
    disp('read and CC sig');
    %%
    load(THLR_FN,'UpdatedB0MapTHLR','UpdatedT2SMap_msTHLR','PDBase0THLR');

c0=PDBase0THLR;
b0=UpdatedB0MapTHLR;
t0=UpdatedT2SMap_msTHLR;

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
SigToUse=DataCCP;

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

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
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
ElemsAlphas=[1e-4, 1e-2, 1e+4, 1e+6];
% ElemsLambda=[1e-5,1e+3,1e-9,1e-10];
ElemsLambda=[1e-5,1e-2,1e-9,1e-10];

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
BARTCmd{SliI}=bartCmd(['splitProx -i 10000 -s 600 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,OneFN,BARTS_Aop.Others{:});
%%
end
%%
BARTCmd=strrep(BARTCmd.','-s 60 ','-s 600 ');
%%
BartCmdFN=[mainP filesep 'splitProxRun'];
fid=fopen(BartCmdFN,'w+');
fprintf(fid,'# splitProx run\r\n');
fprintf(fid,'echo "Starting!"\r\n');
for i=1:numel(BARTCmd)
    fprintf(fid,'%s\n',BARTCmd{i});
end
fclose(fid);
[status,msg,msgID] = fileattrib(BartCmdFN,'+x');
disp('Prepared script for splitProx runs');

% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSRPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigFN,BARTS_Aop.Others{:});
% QQ=bart(['splitProx -i 100 -s 60 -d 2 -g -F ' CurSPrefix ' ' LS_ScriptFN],BARTS_Aop.ImSz16,SigToUse,BARTS_Aop.Others{:});
%%
SliI=13;
for SliI=1:nSlices
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
SliI=15;

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
