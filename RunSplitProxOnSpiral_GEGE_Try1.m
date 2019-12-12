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
TimePoints_ms_GE2=TimePoints_Medb_ms(TimePoints_Medb_ms>TE_SE)-TE_SE;
TimePoints_ms_GEGE=[TimePoints_Med_ms TimePoints_ms_GE2];
TEs_ms=TimePoints_ms_GEGE;
TEs_ms3=perm32(TEs_ms);
nEchos_GEGE=numel(TEs_ms);
NTEs=Col(1:nEchos_GEGE);
EchoTimes_ms=TEs_ms;
EchoTimes_ms3=permute32(EchoTimes_ms);
%
nGE1=numel(TimePoints_Med_ms);
nGE2=numel(TimePoints_ms_GE2);
nT2Q=numel(TimePoints_Medb_ms)-nGE2;
GEGEIdxsB=[ones(1,nGE1) zeros(1,numel(TimePoints_Medb_ms)-nGE2) ones(1,nGE2)]>0;
% T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos_GEGE,TS_Dim);
T=cat(CS_Dim,perm71([ones(nGE1,1);zeros(nGE2,1)]),perm71([zeros(nGE1,1);ones(nGE2,1)]));
T2Sx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(Col(TimePoints_ms_GEGE),[TS_Dim 1]);
Tx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(Col(TimePoints_ms_GEGE),[TS_Dim 1]);
Tx=repmat(Tx,[1 1 1 1 2 1 1 1 1 1]);
T2Sx=repmat(T2Sx,[1 1 1 1 2 1 1 1 1 1]);
M = fmac(ones(nEchos_GEGE,1), T,[Ch_Dim M_Dim],[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos_GEGE), 1i*T,[Ch_Dim M_Dim],[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos_GEGE)/1000, -1i*2*pi*Tx/1000,Ch_Dim,[CS_Dim TS_Dim]);
T2S = fmac(NTEs, -T2Sx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%%
MaxB0=400;
M0MedianFac=20;
T2SRange=[4 200];

WhichRepsToUse=RepsToUse;

RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');

OneFN='/autofs/cluster/kawin/Gilad/One';
writecfl(OneFN,1);

ImSz16=FillOnesTo16(Sz);
ImSz16(TS_Dim)=nEchos_GEGE;
ImSz16(CS_Dim)=2;
ImSz16FN=[mainP filesep 'ImSz16'];
writecfl(ImSz16FN,ImSz16);
%%
clear BARTCmd
%%
c0=cat(CS_Dim,PDBase_GE,PDBase0_SES.*(TSCxPMedOnlyB0S_GESE(:,:,nGE1+nT2Q+1)));

m0=abs(c0);
% m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

% t0s=UpdatedT2SMap_ms_SES;
t2s0=(abs(UpdatedT2SMap_ms_GE)+abs(UpdatedT2SMap_ms_SES))/2;
t0s=t2s0;
t0s=min(T2SRange(2),max(T2SRange(1),t0s));

b0=B0VarAvg+CurB0;

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
Elem0={m0,p0,b0,t0s};
ElemL={M,P,B,T2S};
ElemsL={T,1i*T,-1i*2*pi*Tx/1000,-T2Sx};
disp('ok initials');

Mm0=ElemL{1}*Elem0{1};
expPp0 = exp(ElemL{2} * Elem0{2});
expBb0 = exp(ElemL{3} * Elem0{3});
expTs0 = exp(ElemL{4} * (1./Elem0{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTs0;
Rec0X=sum(Rec0,CS_Dim);
%
% SigToUse=DataCCP(:,:,WhichRepsToUse,:);
SigToUse=DataCCP_GESE(1,:,RepsToUse,1:nccToUse,1,1,1,1,:);

CurTSB=TSB_THLRP_GESE(:,:,:,:,:,:,GEGEIdxsB);
BARTS_Aopx.ImSz16=ImSz16FN;

SensFN=[mainP filesep 'Sens_S' num2str(SliI)];
STrajFN=[mainP filesep 'TrajFN_GEGE_' RepsStr];
TSBFN=[mainP filesep 'TSB_GEGE'];
KernsFN=[mainP filesep 'KernsGEGE_' RepsStr '_S' num2str(SliI)];

writecfl(SensFN,SensCC(:,:,1,1:nccToUse));
writecfl(STrajFN,CurSTraj_GESE);
writecfl(TSBFN,CurTSB);
writecfl(KernsFN,CurKerns_GESE(:,:,:,:,:,:,GEGEIdxsB));

BARTS_Aopx.Others={SensFN STrajFN TSBFN OneFN KernsFN};

ksp_adjFN=[mainP filesep 'sig_adj_' RepsStr '_S' num2str(SliI)];
if(~exist([ksp_adjFN '.cfl'],'file'))
    ksp_adj=bart(['linopScript -d 5 -A ' LS_ScriptFN],BARTS_Aopx.ImSz16,SigToUse,BARTS_Aopx.Others{:});
    disp('got ksp_adj');
    writecfl(ksp_adjFN,ksp_adj);
end

% AHA_Rec0=bart(['linopScript -d 5 -N ' LS_ScriptFN],BARTS_Aopx.ImSz16,Rec0,BARTS_Aopx.Others{:});
%%
for i=1:numel(ElemsL)
    writecfl([CurSRPrefix 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[gsize(Elem0{1},1:2) 1 1 1 1 1]));
end
%%
Elem0={m0*0,p0*0+1.5,b0,t0s*0+40};
%%
CurSRPrefix=[mainP filesep 'splitProx_' RepsStr '_S' num2str(SliI) filesep];
gmkdir(CurSRPrefix);

CurSuffix='Lm4c_';
CurSRPrefixOut=[CurSRPrefix CurSuffix];

ElementTypes=[1 2 3 4];
ElemsAlphas=[1e-4, 1e+0, 1e+3, 1e+4];
ElemsAlphas=[1e-4, 1e+0, 1e+3, 1e+4];
ElemsLambda=[1e-6,1e-7,1e-9,1e-10];
ninneriterBART=[2 2 2 2];

writecfl([CurSRPrefix 'ElemsAlpha'],ElemsAlphas.');
writecfl([CurSRPrefix 'ElemsLambda'],ElemsLambda.');
writecfl([CurSRPrefix 'ElementTypes'],ElementTypes.');
writecfl([CurSRPrefix 'ninneriter'],ninneriterBART);
disp(['saved hyperparams']);
%%
for i=1:numel(Elem0)
    writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Elem0{i});
end

disp(['saved warmstart ' CurSRPrefix]);

for i=1:4
    gdelete([CurSRPrefix 'Elem' num2str(i-1) '.hdr']);
    gdelete([CurSRPrefix 'Elem' num2str(i-1) '.cfl']);
end
% %% Continue?
ContinueRun=false;
if(ContinueRun)
    for i=1:numel(Elem0)
        writecfl([CurSRPrefix 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end

BARTCmd{SliI,rs}=bartCmd(['splitProx -i 10000 -I 300 -s 600 -d 2 -T 4000 -g -f -S ' ksp_adjFN ' -F ' CurSRPrefix ' -O ' CurSRPrefixOut ' ' LS_ScriptFN],BARTS_Aopx.ImSz16,OneFN,BARTS_Aopx.Others{:});
% system(BARTCmd{SliI,rs});
% end
% end
disp('Prepared BARTCmd');
%%
system(BARTCmd{SliI,rs});
%%
% CurSRPrefix=[mainP filesep 'splitProx_' RepsStr '_S' num2str(SliI) filesep];
% CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];

ErrVec=readcfl([CurSRPrefixOut 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);

DD=dir([CurSRPrefixOut 'Elem0*.hdr']);
[~,MI]=max([DD.datenum]);
FFN=[DD(MI).folder filesep DD(MI).name];
LastSuffix=FFN(numel(CurSRPrefixOut)+6:end-4);
for i=1:4
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter900']);
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) LastSuffix]);
end

MmM=ElemL{1}*Maps{1};
expPpM = exp(ElemL{2} * Maps{2});
expBbM = exp(ElemL{3} * Maps{3});
expTtM = exp(ElemL{4} * (1./Maps{4}));
RecM=MmM.*expPpM.*expBbM.*expTtM;
RecMX=squeeze(sum(RecM,CS_Dim));

disp(['Loaded maps ' LastSuffix]);
%%
for SliI=1:nSlices
CurSRPrefix=[mainP filesep 'gB0xSli' num2str(SliI) '_R' RepsStr '_'];
% CurSRPrefixOut=[CurSRPrefix 'Lm45b_'];
CurSRPrefixOut=[CurSRPrefix 'Lm4c_'];
disp(CurSRPrefixOut);

DD=dir([CurSRPrefixOut 'Elem*.hdr']);
for i=1:4
    Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1)]);
%     Maps{i}=readcfl([CurSRPrefixOut 'Elem' num2str(i-1) '_iter1000']);
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
