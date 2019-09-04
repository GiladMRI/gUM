THLRMultiShotM=permute(cat(4,THLRMultiShotS{:}),[1 2 7 4 6 5 3]);
fgmontage(grmss(THLRMultiShotM(:,:,:,ROrd),3))
%%
IdxsForCusack=[2 3];
CusackDeltaTE_us=diff(IdxsForCusack)*numel(TrajPartMed)*AcqDwellTime_us/(size(THLRMultiShotM,3)-1);
% Mag=mean(abs(DataForSlice{7}{6}(:,:,IdxsForCusack)),3);
Mag=mean(abs(THLRMultiShotM(:,:,IdxsForCusack,SliI)),3);
% Combined=Mag.*exp(1i*(angle(DataForSlice{7}{6}(:,:,IdxsForCusack(1)))-angle(DataForSlice{7}{6}(:,:,IdxsForCusack(2)))));
Combined=Mag.*exp(1i*(angle(THLRMultiShotM(:,:,IdxsForCusack(1),SliI))-angle(THLRMultiShotM(:,:,IdxsForCusack(2),SliI))));
dAngle=double(angle(Combined));
[unwrapped] = cusackUnwrap(repmat(dAngle,[1 1 2]), repmat(Mag,[1 1 2])/3000);
B0_HzUx=-unwrapped(:,:,1)/(2*pi*CusackDeltaTE_us/1e6);
fgmontage(B0_HzUx,[-400 400])

dAngleX=angle(THLRMultiShotM(:,:,2:end,SliI).*conj(THLRMultiShotM(:,:,1:end-1,SliI)));

Combined=Mag.*exp(1i*angle(R1_tmp));
%%
fgmontage(B0M2(:,:,SliRefMI(SliI)),[-300 300])

%%
[~,~,~,H_AllTS1_2]=ghankel(nTSMed,2,TrgSz);
for s=1:nSlices
    XX=permute(THLRMultiShotM(:,:,:,s),[1 2 7 6 5 4 3]);
    [ ~, s_tmp, V_tmp] = batch_svd(H_AllTS1_2*squeeze(XX));
    R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
    InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
    T2SMap_tmp_ms=-InnerTSDiff_ms./log(abs(R1_tmp));
    B0Map_tmp=-(angle(R1_tmp)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
    PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed-1)))));
    W=abs(R1_tmp).^permute32(0:(nTSMed-1));
    PDEstx=sum(PDMEx,3)./sum(W,3);
    
    MPBDS(:,:,1,s)=abs(PDEstx)/(15e-4);
    MPBDS(:,:,2,s)=(angle(PDEstx)+pi)/(2*pi);
    MPBDS(:,:,3,s)=abs(T2SMap_tmp_ms)/100;
    MPBDS(:,:,4,s)=(B0Map_tmp+100)/200;
end
MPBDS(:,:,1,:)=MPBDS(:,:,1,:)./grmss(MPBDS(:,:,1,:))/2;
MPBDS=MPBDS.*permute43(BNS).*permute43(MskS);
fgmontagex(MPBDS,[0 1]);
fgmontagex(MPBDS(:,:,:,ROrd),[0 1]);
%% Med2 slice independent stuff
AbsB0Bins=-0.5:1:500.5;
AbsB0BinCenters=(AbsB0Bins(1:end-1)+AbsB0Bins(2:end))/2;
    
for s=1:nSlices
    absB0WHist=gWeightedHist2D(abs(B0M2(:,:,SliRefMI(s))),AbsB0Bins,RefMg(:,:,SliRefMI(s)));
    cNabsB0WHist=cumsum(absB0WHist)./sum(absB0WHist);
    MaxAbsB0=AbsB0BinCenters(find(cNabsB0WHist>0.97,1));
    MaxAngleForTS=45;
    TSDiffNeededS(s)=(1000/MaxAbsB0)*(MaxAngleForTS/360);
end
%%
TSDiffNeeded=1;
dTS_planned4_ms=TSDiffNeeded;
nTSMed4=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned4_ms);
InnerMedTS4Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed4-1);
TimePointsMed4_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed4);
TimePointsMed4_ms3=permute32(TimePointsMed4_ms);

TSBMed4=GetTSCoeffsByLinear(nPointsMed,nTSMed4);
TSBPMed4=permute(TSBMed4,[3 1 4 5 6 7 2]);
%
for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    KernsByRepMed4{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRepC{CurRep},TSBMed4);
    KernsPMed4C{CurRep}=permute(KernsByRepMed4{CurRep},[1 2 7 6 5 4 3]);
end
KernsPMMed4=cat(3,KernsPMed4C{:});
clear KernsPMed4C
disp('Prepared Med4');
%% 2L-THLR
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
BaseB0ToUse=tt;
% BaseB0ToUse=B0M2(:,:,SliRefMI(SliI));

MaxdB0ToHandle=20;
dTSNeededFordB0=1000/(MaxdB0ToHandle*(360/MaxAngleForTS));
nTSNeededFordB0=TimePointsMed4_ms3(end)/dTSNeededFordB0;

nTS_THLR=15;
TSCxPMed4OnlyB0=permute73(exp(1i*2*pi*(1e-3)*BaseB0ToUse.*TimePointsMed4_ms3));
Comps=GetTSCoeffsByLinear(nTSMed4,nTS_THLR);
Sz16_2LTHLR=FillOnesTo16([Sz 1 1 1 nTS_THLR]);
CompsP=permute(Comps,[7:-1:3 2 1]);
%
THLR2L_lambda=10^1;
% RhoStr=[' -u ' num2str(0.01) ' '];
RhoStr=[' -u ' num2str(1e-1) ' '];
% RhoStr=[' -u ' num2str(1e-3) ' ']; % in THLRMultiShot

disp('Prepared THLR-2L');
%% allreps ncc3: 113sec ncc31 518sec almost same result; Now 58sec
% nccToUse=31;
nccToUse=3;
CurRep=1:(nRepsHdr-1);
%     -w 1
Rec_THLR2L=bart(['picsS -m ' RhoStr '  -R K:32:3:' num2str(THLR2L_lambda) ':2:1:0:5 ' ScriptFN_CompgBo],...
    Sz16_2LTHLR,DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
    SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed4,TSCxPMed4OnlyB0,...
    sum(KernsPMMed4(:,:,CurRep,:,:,:,:),3),CompsP);
Rec_THLR2LX=sum(Rec_THLR2L.*CompsP,6).*TSCxPMed4OnlyB0;
%%
[MPBD, B0_THLR2L]=getMPBDBySVDTH(Rec_THLR2LX,InnerMedTS4Diff_ms,FirstEchoTE_ms);
MPBD(:,:,1)=(MPBD(:,:,1)-0.5)/2+0.5;
fgmontagex(MPBD,[0 1]);
%% end 2L-THLR
%% LLR B0-Var with new B0
Mag_THLRMultiShot=squeeze(abs(Rec_THLR2LX));
Mag_THLRMultiShot2D=reshape(Mag_THLRMultiShot,prod(gsize(Mag_THLRMultiShot,1:2)),size(Mag_THLRMultiShot,3));
[MUd,MSd,MVd]=svd(Mag_THLRMultiShot2D(:,:),'econ');
%%
[~,~,~,H_AllTS]=ghankel(nTSMed4,2,TrgSz);
[ ~, s_LLR_AllTS, V_LLR_AllTS] = batch_svd(H_AllTS*squeeze(Rec_THLR2LX));
R1ts=V_LLR_AllTS(:,:,2,1)./V_LLR_AllTS(:,:,1,1); % R1 is simply the decay
UpdatedT2SMap_ms1=-InnerMedTS4Diff_ms./log(abs(R1ts));
UpdatedB0Map1=-(angle(R1ts)/(2*pi))/(InnerMedTS4Diff_ms/1e3); % in Hz
figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR-2L multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
%%
TempRMS=grmss(Rec_THLR2LX,3:30);
WHT2s=HT2s*0;
for i=1:numel(HT2s)
    WHT2s(i)=sum(TempRMS(T2sbin==i));
end
SWHT2s=max(WHT2s,sum(WHT2s)*0.03/numel(WHT2s));

clear WDecays4
for i=2:numel(HT2s)-1
    WDecays4(i,:)=SWHT2s(i)*exp(-TimePointsMed4_ms./T2sCenters(i-1));
end
[WUd4,WSd4,WVd4]=svd(WDecays4,'econ');
%%
UpdatedB0Map=UpdatedB0Map1;
UpdatedT2SMap_ms=UpdatedT2SMap_ms1;
R1=R1ts;
R1x=min(max(abs(R1),0.3),1).*exp(-1i*angle(R1));
%
BaseComponents=WVd4(:,1:nComponentsToUse);

nComponentsToUsewB0Var=3;
B0VarEst=2;
        
ComponentsToUse=BaseComponents(:,1:nComponentsToUsewB0Var);
B0VarVec=exp(1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed4_ms.'));
B0VarVecm=exp(-1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed4_ms.'));

Comps=[ComponentsToUse ComponentsToUse.*B0VarVec ComponentsToUse.*B0VarVecm];

nComponentsIncB0Var=size(Comps,2);

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsIncB0Var]);
CompsP=permute(Comps,[7:-1:3 2 1]);

LLR_lambda=10^1;
RhoStr=[' -u ' num2str(0.01) ' '];
BlkSz=3;
disp('Prepared LLR Med4');
%%
nccToUse=31;
RepsToCalc=[1 2 3 4 5 6 37 38 39];
for CurRep=RepsToCalc
    Rec_CompgB0_Med4_WB0Var=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
        Sz16CompgB0,DataCCP(:,TrajPartMed,CurRep,1:nccToUse),SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),...
        TSBPMed4,TSCxPMed4OnlyB0,KernsPMMed4(:,:,CurRep,:,:,:,:),CompsP);
    Rec_CompgB0_Med4_C_WB0VarCS{CurRep,SliI}=Rec_CompgB0_Med4_WB0Var;
end
%%
Rec_CompgB0_Med4_C_WB0VarCSM=cat(4,Rec_CompgB0_Med4_C_WB0VarCS{:,SliI});
Rec_CompgB0_Med4_C_WB0VarCSMX=permute73(sum(Rec_CompgB0_Med4_C_WB0VarCSM.*CompsP,6).*TSCxPMed4OnlyB0);
for i=1:size(Rec_CompgB0_Med4_C_WB0VarCSMX,4)
    [MPBDR(:,:,:,i), UpdatedB0Map_LLR_B0VarR4(:,:,i), s_tmpR(:,:,:,i)]=getMPBDBySVDTH(Rec_CompgB0_Med4_C_WB0VarCSMX(:,:,:,i),InnerMedTS4Diff_ms);
end
B0sR=(MPBDR(:,:,1,:)-0.5)*200;
UpdatedB0Map_LLR_B0VarRX(:,:,RepsToCalc)=UpdatedB0Map_LLR_B0VarR4;
%%
[MPBD, B0_THLRLLR_SingleS]=getMPBDBySVDTH(Rec_CompgB0_Med4_C_WB0VarCSMX(:,:,:,1),InnerMedTS4Diff_ms,FirstEchoTE_ms);
MPBD(:,:,1)=(MPBD(:,:,1)-0.5)/2+0.5;
fgmontagex(MPBD,[0 1]);
title(TtlX);
%% end LLR B0-Var with new B0
%% Now LLR some-shots based on B0Var
nComponentsToUseAfterB0Var=2;
% Comps=BaseComponents(:,1:nComponentsToUseAfterB0Var);
Comps=WVd4(:,1:nComponentsToUseAfterB0Var);
% MVd(:,1:nComponentsToUse);

CompsP=permute(Comps,[7:-1:3 2 1]);
Sz16CompgB0VarSomeShots=FillOnesTo16([Sz 1 1 1 nComponentsToUseAfterB0Var]);

WhichTSToUse_LLRSomeShots=7:nTSMed4;
WhichTSToUse_LLRSomeShotsB=zeros(1,nTSMed4);
WhichTSToUse_LLRSomeShotsB(WhichTSToUse_LLRSomeShots)=1;
WhichTSToUse_LLRSomeShotsB7=perm72(WhichTSToUse_LLRSomeShotsB);

%% Prepare updated B0 TSC
CurTSCxPMed4OnlyB0AfterPerSlice=perm93(perm73(exp(1i.*2*pi*Perm93(UpdatedB0Map_LLR_B0VarR4)*(1e-3).*TimePointsMed4_ms3)));
%% 300/123? sec per 1.
nShotsTogether=3;
RIdxs=1:nShotsTogether:6;
for r=1:numel(RIdxs)
    CurRepa=RIdxs(r);
    CurRep=CurRepa+(0:(nShotsTogether-1));
    Rec_CompgB0_C_gB0Var4_SomeShotsS{r,SliI}=bart(['picsS -d 3 -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgB0VarSomeShots],...
        Sz16CompgB0VarSomeShots,perm93(DataCCP(:,TrajPartMed,CurRep,1:nccToUse)),...
        SensCC(:,:,:,1:nccToUse),perm93(STraj3MMed(:,:,CurRep)),TSBPMed4.*WhichTSToUse_LLRSomeShotsB7,perm93(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,CurRep,:,:,:,:,:,:)),...
        perm93(KernsPMMed4(:,:,CurRep,:,:,:,:).*WhichTSToUse_LLRSomeShotsB7),CompsP);
end
%%
Rec_CompgB0_C_gB0Var4_SomeShotsSM=cat(4,Rec_CompgB0_C_gB0Var4_SomeShotsS{:,SliI});
Rec_CompgB0_C_gB0Var4_SomeShotsSMX=sum(Rec_CompgB0_C_gB0Var4_SomeShotsSM.*CompsP,6);
%% End LLR some-shots based on B0var
%%
nTSMed2=nTSMed*4+1;
TimePointsMed2_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed2);
TimePointsMed2_ms3=permute32(TimePointsMed2_ms);

TSBMed2=GetTSCoeffsByLinear(nPointsMed,nTSMed2);
TSBPMed2=permute(TSBMed2,[3 1 4 5 6 7 2]);
%
for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
    KernsByRepMed2{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRep,TSBMed2);
    KernsPMed2C{CurRep}=permute(KernsByRepMed2{CurRep},[1 2 7 6 5 4 3]);
end
KernsPMMed2=cat(3,KernsPMed2C{:});
clear KernsPMed2C
disp('Prepared Med2');
%%
SliI=7;
%%
SelfSens1=DataForSlice{SliI}{9};
sccmtx=DataForSlice{SliI}{10};

ADataIsL=ADatax.image(:,:,:,:,SliI,3,:,:,RepsToRead,:,:,:,:,:,:,:,:);
ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);

ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
for i=1:ncc
    ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
end
DataC=permute(ADataIsLCC,[1 3 2]);
disp('ok cc');
DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
DataCCP=permute(DataPC,[1:3 8 4:7]);

SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
%%
for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    SnufftStruct_CurRepC{CurRep} = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),Sz), Sz, [6 6], Sz*2); % st.om
end
%% Redo THLR
dTS_planned3_ms=2.2;
nTSMed3=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned3_ms);

TSBMed3=GetTSCoeffsByLinear(nPointsMed,nTSMed3);
TSBPMed3=permute(TSBMed3,[3 1 4 5 6 7 2]);

TimePointsMed3_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed3);
TimePointsMed3_ms3=permute(TimePointsMed3_ms,[1 3 2]);

for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    KernsByRepMed3{CurRep}=NUFFT_to_Toep_2blocks(SnufftStruct_CurRepC{CurRep},TSBMed3);
    KernsPMed3C{CurRep}=permute(KernsByRepMed3{CurRep},[1 2 7 6 5 4 3]);
end
KernsPMMed3=cat(3,KernsPMed3C{:});
clear KernsPMedC
disp('Prepared Med3');
%% Try to get all TSC, multi-shot, for B0,T2*
nccToUse=3;
CurReps=1:39;

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC3=FillOnesTo16(TrgSz);
Sz16AllTSC3(7)=nTSMed3;

THLR_lambda=10;
% RhoStr=[' -u ' num2str(1e-3) ' '];
RhoStr=[' -u ' num2str(1e-1) ' '];

THLRMultiShot3=bart(['picsS -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16AllTSC3,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed3,1,sum(KernsPMMed3(:,:,CurReps,:,:,:,:),3));
%%
[~,~,~,H_MedTS3_2]=ghankel(nTSMed3,2,TrgSz);
XX=squeeze(THLRMultiShot3);
[ ~, s_tmp, V_tmp] = batch_svd(H_MedTS3_2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
InnerMedTS3Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed3-1);
T2S_Med3_tmp_ms=-InnerMedTS3Diff_ms./log(abs(R1_tmp));
B0_Med3_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS3Diff_ms/1e3); % in Hz
PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed3-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed3-1));
PDEstx=sum(PDMEx,3)./sum(W,3);

figure;subplot(1,2,1);gmontage(B0_Med3_tmp,[-100 100]*2);title('B_0,T_2^* from THLR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(T2S_Med3_tmp_ms,[0 100]);removeTicks;
%%
[~,~,~,H_MedTS3_3]=ghankel(nTSMed3,3,TrgSz);
% DUse2ndSVMap=imdilate(Use2ndSVMap,strel('disk', 2),'same');
DUse2ndSVMap=SmoothBySlices(Use2ndSVMap,[3 3],1.5)>(1/9);

THLRFitted=FitToModel_THLR2Mskf(THLRMultiShot3,H_MedTS3_3,DUse2ndSVMap);
%%
THLRFittedFN=[mainP filesep 'THLRFitted'];
writecfl(THLRFittedFN,permute73(THLRFitted));
THLR_lambda=10;
THLRMultiShot3x=bart(['picsS -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 -W ' THLRFittedFN ' ' ScriptFN_AllTS],Sz16AllTSC3,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed3,1,sum(KernsPMMed3(:,:,CurReps,:,:,:,:),3));
%%
sTwixRef=load('/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/meas_MID01104_FID09968_gre_te7_40/sTwixX.mat');
sTwixRef=sTwixRef.sTwixX;
RefFOV=sTwixRef.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
QQ=load('/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/meas_MID01104_FID09968_gre_te7_40/B0M2.mat');
RefLocs=load('/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/meas_MID01104_FID09968_gre_te7_40/Locs.mat');
RefLocs=RefLocs.RotatedLocs;
QQF=load('/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/meas_MID01104_FID09968_gre_te7_40/FirstEcho.mat');
RefMg=grmss(QQF.FirstEcho,3);
RefMg=padarray(RefMg,[0 2 0]','Both');
RefMg=rot90(RefMg);
RefMg=imresizeBySlices(RefMg,TrgSz);
B0M2=padarray(QQ.B0M2,[0 2 0]','Both');
B0M2=rot90(B0M2);
B0M2=imresizeBySlices(B0M2,TrgSz);

for i=1:nSlices
    [SliRefdMI(i),SliRefMI(i)]=min(abs(RefLocs(3,:)-RotatedLocs(3,Ord(i))));
end
%%
tt=B0_Med3_tmp;
tt(40:80,70:110)=B0M2(44:84,74:114,17);
fgmontage(tt,[-250 250])
%% End THLR
%%
THLRMultiShot=THLRMultiShotS{SliI};

Mag_THLRMultiShot=squeeze(abs(THLRMultiShot));
Mag_THLRMultiShot2D=reshape(Mag_THLRMultiShot,prod(gsize(Mag_THLRMultiShot,1:2)),size(Mag_THLRMultiShot,3));
[MUd,MSd,MVd]=svd(Mag_THLRMultiShot2D(:,:),'econ');
%%
[~,~,~,H_AllTS]=ghankel(nTSMed,2,TrgSz);
[ ~, s_LLR_AllTS, V_LLR_AllTS] = batch_svd(H_AllTS*squeeze(THLRMultiShot));
R1ts=V_LLR_AllTS(:,:,2,1)./V_LLR_AllTS(:,:,1,1); % R1 is simply the decay
InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
UpdatedT2SMap_ms1=-InnerTSDiff_ms./log(abs(R1ts));
UpdatedB0Map1=-(angle(R1ts)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
% figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR multi-shot');removeTicks;colorbar
% subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
%%
% UpdatedB0Map=UpdatedB0Map1;
% UpdatedT2SMap_ms=UpdatedT2SMap_ms1;
% R1=R1ts;
% R1x=min(max(abs(R1),0.3),1).*exp(-1i*angle(R1));
% InnerShotDiff_ms=InnerTSDiff_ms;

% TSCxMed3=R1x.^(TimePointsMed3_ms3/InnerShotDiff_ms);
% TSCxPMed3=permute73(TSCxMed3);

ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
% TSCxPMed3OnlyB0=exp(1i.*angle(TSCxPMed3));

TSCxPMed3OnlyB0=permute73(exp(1i*2*pi*(1e-3)*tt.*TimePointsMed3_ms3));


%
BaseComponents=MVd(:,1:nComponentsToUse);
% BaseComponents2=interp1(TimePointsMed_ms,BaseComponents,TimePointsMed2_ms);
BaseComponents3=interp1(TimePointsMed_ms,BaseComponents,TimePointsMed3_ms);

% nComponentsToUsewB0Var=3;
% B0VarEst=1e3/InnerTSDiff_ms; % 2;
        
% ComponentsToUse=BaseComponents2(:,1:nComponentsToUsewB0Var);
% B0VarVec=exp(1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed2_ms.'));
% B0VarVecm=exp(-1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed2_ms.'));

Comps=BaseComponents3;

% B0VarVec2=exp(1i.*2*pi*2*B0VarEst*(1e-3).*(TimePointsMed2_ms.'));

% Comps=[ComponentsToUse ComponentsToUse.*B0VarVec ComponentsToUse.*B0VarVecm];
% Comps=[ComponentsToUse ComponentsToUse.*B0VarVec];
% Comps=[ComponentsToUse ComponentsToUse.*B0VarVec ComponentsToUse.*B0VarVec2];

nComponentsIncB0Var=size(Comps,2);

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsIncB0Var]);
CompsP=permute(Comps,[7:-1:3 2 1]);
%
LLR_lambda=10^1;
% LLR_lambda=10^0;
RhoStr=[' -u ' num2str(0.01) ' '];
BlkSz=3;
disp('Prepared LLR');
%% 180 sec/300 sec/ 77 sec
nccToUse=31;
nccToUse=3;
CurRep=1:(nRepsHdr-1);
% for CurRep=1:(nRepsHdr-1)
% for CurRep=1:nRepsToCalc
    Rec_LLR_B0Wraps_MS=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
        Sz16CompgB0,DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed3,TSCxPMed3OnlyB0,...
        sum(KernsPMMed3(:,:,CurRep,:,:,:,:),3),CompsP);
    Rec_LLR_B0Wraps_MSX=sum(Rec_LLR_B0Wraps_MS.*CompsP,6).*TSCxPMed3OnlyB0;
%     Rec_LLR_B0Wraps_MS=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
%         Sz16CompgB0,DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
%         SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed2,TSCxPMed2OnlyB0,...
%         sum(KernsPMMed2(:,:,CurRep,:,:,:,:),3),CompsP);
%     Rec_LLR_B0Wraps_MSX=sum(Rec_LLR_B0Wraps_MS.*CompsP,6).*TSCxPMed2OnlyB0;
%     Rec_CompgB0_C_WB0VarCS{CurRep,SliI}=Rec_CompgB0_C_WB0Var;
% end
%%
[~,~,~,H_MedTS3_2]=ghankel(nTSMed3,2,TrgSz);
XX=squeeze(Rec_LLR_B0Wraps_MSX);
[ ~, s_tmp, V_tmp] = batch_svd(H_MedTS3_2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
InnerMedTS3Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed3-1);
T2S_Med3_tmp_ms=-InnerMedTS3Diff_ms./log(abs(R1_tmp));
B0_Med3_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS3Diff_ms/1e3); % in Hz
PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed3-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed3-1));
PDEstx=sum(PDMEx,3)./sum(W,3);

figure;subplot(1,2,1);gmontage(B0_Med3_tmp,[-100 100]*2);title('B_0,T_2^* from THLR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(T2S_Med3_tmp_ms,[0 100]);removeTicks;

%%
[~,~,~,H_MedTS2_2]=ghankel(nTSMed2,2,TrgSz);
XX=squeeze(Rec_LLR_B0Wraps_MSX);
[ ~, s_tmp, V_tmp] = batch_svd(H_MedTS2_2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
InnerMedTS2Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed2-1);
T2S_Med2_tmp_ms=-InnerMedTS2Diff_ms./log(abs(R1_tmp));
B0_Med2_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS2Diff_ms/1e3); % in Hz
PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed2-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed2-1));
PDEstx=sum(PDMEx,3)./sum(W,3);

figure;subplot(1,2,1);gmontage(B0_Med2_tmp,[-100 100]*2);title('B_0,T_2^* from LLR multi-shot w B0Wrap components');removeTicks;colorbar
subplot(1,2,2);gmontage(T2S_Med2_tmp_ms,[0 100]);removeTicks;
%%
figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
%%
dAngle=angle(R1_tmp);
Mg=s_tmp(:,:,1);
[unwrapped] = cusackUnwrap(repmat(dAngle,[1 1 2]), repmat(Mg,[1 1 2]));
unwrapped=unwrapped(:,:,1);
%%
MskToUse=DMskS(:,:,SliI);

sRatio=(s_tmp(:,:,2)./s_tmp(:,:,1)).*MskToUse.*BNS(:,:,SliI);
% Use2ndSVMap=sRatio>0.1;
Use2ndSVMap=imfillholesBySlices( sRatio>0.1);
disp('Got Use2ndSVMap');
%%
nccToUse=3;
CurReps=1:(nRepsHdr-1);
%%
nccToUse=31;
CurReps=1;
%%
nccToUse=31;
CurReps=2;
%%
Sz16AllTSC2=Sz16AllTSC;
Sz16AllTSC2(7)=nTSMed2;

HankelOrder=3;
%%
TSBPMed2=TSBPMed3;
KernsPMMed2=KernsPMMed3;
Sz16AllTSC2=Sz16AllTSC3;
TSCxPMed2OnlyB0=TSCxPMed3OnlyB0;
nTSMed2=nTSMed3;
H_MedTS2_2=H_MedTS3_2;
H_MedTS2_3=H_MedTS3_3;
%% Initialize X0
for CurReps=1:(nRepsHdr-1)
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}(:,:,:,:,:,1:3).*CompsP(:,:,:,:,:,1:3,:),6);
% % Using THLR-Multishot warmstart
% X0NoB=THLRMultiShot.*conj(TSCxPMedOnlyB0);
% X0NoB=Rec_LLR_B0Wraps_MSX.*conj(TSCxPMed2OnlyB0);

% X0NoB=X0N.*conj(TSCxPMed2OnlyB0);
X0NoB=X0N_MS.*conj(TSCxPMed2OnlyB0);
%
disp('Running AHy');
AHy=bart(['linopScript -A ' ScriptFN_AllTS],Sz16AllTSC2,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
    SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed2,1,...
    sum(KernsPMMed2(:,:,CurReps,:,:,:,:),3));

X0=X0NoB.*TSCxPMed2OnlyB0;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC2,X0,...
    SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed2,1,...
    sum(KernsPMMed2(:,:,CurReps,:,:,:,:),3));
Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*Fac;
X0NoBN=X0NoB.*Fac;
disp('Initialized X0');

% AHA_Fac=grmss(AHAX0)./grmss(X0); % Poor effective Lipshitz factor
rcur=AHy-AHAX0N;
AHAr=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC2,rcur,...
    SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed2,1,...
    sum(KernsPMMed2(:,:,CurReps,:,:,:,:),3));
AHA_Fac=grmss(AHAr)./grmss(rcur); % better Poor effective Lipshitz factor

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;

[~,~,~,H_AllTS2]=ghankel(nTSMed2,HankelOrder,TrgSz);

iter=1;

X0new=squeeze(X0NoBN).*MskToUse;
alpha=0.01/AHA_Fac;

nIters=30;
MinDiffT=3e-4;
% Wlambda_base=0.1;
% Wlambda=0.0003;
% WT=wave_thresh('db4', 3, Wlambda_base);

disp('Ready to iterate');
%% Iterate
for jj=1:nIters
    disp(jj);
    AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC2,X0N,...
        SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed2,1,...
        sum(KernsPMMed2(:,:,CurReps,:,:,:,:),3));
    Fac=grmss(AHy)./grmss(AHAX0);
    AHAX0N=AHAX0*Fac;
    rcur=AHy-AHAX0N;
    rcur_rms(iter)=grmss(rcur);
    disp(rcur_rms(iter));
    CurRatioToDiff(iter)=AHy_rms./rcur_rms(iter);
    rcurNoB=rcur.*conj(TSCxPMed2OnlyB0);
    UpdatedX0=X0new+alpha*squeeze(rcurNoB);
    
    Fitted=FitToModel_THLR2Mskf(UpdatedX0,H_MedTS2_3,DUse2ndSVMap);
%     In=UpdatedX0;
%     FitToModel_THLR2Msk
%     X0new=Out.*MskToUse;
    X0new=Fitted.*MskToUse;
    LastX0N=X0N;
    X0N=permute(X0new,[1 2 7 6 5 4 3]).*TSCxPMed2OnlyB0;
        
    dX0N(iter)=grmss(LastX0N-X0N)./grmss(LastX0N);
    disp(['dX0N(iter) : ' num2str(dX0N(iter))]);
    if(dX0N(iter)<MinDiffT)
        disp('Less the MinDiffT, breaking');
        break;
    end
    
    iter=iter+1;
end
disp('Finished iterations');

X0NRC{min(CurReps)}=X0N;
end
%%
X0NRM=permute73(cat(4,X0NRC{:}));
XX=grmss(X0NRM,3);

XX=squeeze(X0NRM(:,:,20,:));

sXX=std(XX,[],3);
mXX=mean(XX,3);
tXX=mXX./sXX;
%%
% XX=squeeze(X0N);TtlX='THLR Multishot';
% XX=squeeze(X0N);TtlX='THLR Singleshot warmstart MS';
% XX=squeeze(X0N.*exp(-1i*angle(X0N_MS)));TtlX='B0 diff';
XX=squeeze(X0NRC{1}.*exp(-1i*angle(X0NRC{2})));TtlX='B0 diff';
% XX=squeeze(X0NoB.*TSCxPMed2OnlyB0);
[ ~, s_tmp, V_tmp] = batch_svd(H_MedTS2_2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
InnerMedTS2Diff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed2-1);
T2S_Med2_tmp_ms=-InnerMedTS2Diff_ms./log(abs(R1_tmp));
B0_Med2_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS2Diff_ms/1e3); % in Hz
PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed2-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed2-1));
PDEstx=sum(PDMEx,3)./sum(W,3);

MPBD=permute43(PartitionDim(cat(3,(B0_Med2_tmp+100)/200,abs(T2S_Med2_tmp_ms)/100,abs(PDEstx)/grmss(abs(PDEstx))/2.5,(angle(PDEstx)+pi)/2/pi),3,2));
fgmontagex(MPBD,[0 1]);title(TtlX);
%%
figure;
subplot(2,2,1);gmontage(B0_Med2_tmp,[-100 100]);title(['B_0,T_2^* from ' TtlX]);removeTicks;colorbar
subplot(2,2,2);gmontage(abs(T2S_Med2_tmp_ms),[0 100]);removeTicks;
subplot(2,2,3);gmontage(abs(PDEstx),[0 2.4e-3]);title('PD');removeTicks
subplot(2,2,4);gmontage(angle(PDEstx),[-pi pi]);removeTicks;
%%
X0N_MS=X0N;
% figure;subplot(1,2,1);gmontage(B0_Med2_tmp,[-100 100]*2);title('B_0,T_2^* from THLR single-shot cycled');removeTicks;colorbar
% subplot(1,2,2);gmontage(abs(T2S_Med2_tmp_ms),[0 100]);removeTicks;
