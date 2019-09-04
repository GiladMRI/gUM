nRepsToCalc=6;
SliI=2;
for SliI=1:nSlices
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
UpdatedB0Map=UpdatedB0Map1;
UpdatedT2SMap_ms=UpdatedT2SMap_ms1;
R1=R1ts;
R1x=min(max(abs(R1),0.3),1).*exp(-1i*angle(R1));
InnerShotDiff_ms=InnerTSDiff_ms;

TSCxMed=R1x.^(TimePointsMed_ms3/InnerShotDiff_ms);
TSCxPMed=permute(TSCxMed,[1 2 7 6 5 4 3]);

ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
TSCxPMedOnlyB0=exp(1i.*angle(TSCxPMed));
%
BaseComponents=MVd(:,1:nComponentsToUse);

nComponentsToUsewB0Var=3;
B0VarEst=2;
        
ComponentsToUse=BaseComponents(:,1:nComponentsToUsewB0Var);
B0VarVec=exp(1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed_ms.'));
B0VarVecm=exp(-1i.*2*pi*B0VarEst*(1e-3).*(TimePointsMed_ms.'));

Comps=[ComponentsToUse ComponentsToUse.*B0VarVec ComponentsToUse.*B0VarVecm];

nComponentsIncB0Var=size(Comps,2);

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsIncB0Var]);
CompsP=permute(Comps,[7:-1:3 2 1]);

LLR_lambda=10^1;
RhoStr=[' -u ' num2str(0.01) ' '];
BlkSz=3;
disp('Prepared LLR');
%%
% for CurRep=1:(nRepsHdr-1)
for CurRep=1:nRepsToCalc
    Rec_CompgB0_C_WB0Var=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
        Sz16CompgB0,DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed,TSCxPMedOnlyB0,...
        sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP);
    Rec_CompgB0_C_WB0VarCS{CurRep,SliI}=Rec_CompgB0_C_WB0Var;
end
%%
% Rec_CompgB0_C_WB0VarCS(:,SliI)=Rec_CompgB0_C_WB0VarC;
% Rec_CompgB0_C_WB0VarX1=squeeze(sum(Rec_CompgB0_C_WB0Var.*CompsP,6));
%%
Rec_CompgB0_C_WB0VarM=cat(4,Rec_CompgB0_C_WB0VarCS{:,SliI});
Rec_CompgB0_C_WB0VarM1=permute(sum(Rec_CompgB0_C_WB0VarM.*CompsP,6),[1 2 7 4 6 5 3]);
%%
WhichInnerShotsToUse=11:20;
[~,~,~,H_tmp]=ghankel(numel(WhichInnerShotsToUse),2,TrgSz);

for CurRep=1:nRepsToCalc
    CurX=squeeze((Rec_CompgB0_C_WB0VarM1(:,:,WhichInnerShotsToUse,CurRep)));
    [ ~, s_tmp, V_tmp] = batch_svd(H_tmp*CurX);
    R1a=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
    UpdatedT2SMap_LLR_B0VarR(:,:,CurRep,SliI)=-InnerTSDiff_PerShot_ms./log(abs(R1a));
    UpdatedB0Map_LLR_B0VarR(:,:,CurRep,SliI)=-(angle(R1a)/(2*pi))/(InnerTSDiff_PerShot_ms/1e3); % in Hz    
    PD_LLR_B0VarR(:,:,CurRep,SliI)=mean(CurX.*(R1a.^(-permute32(WhichInnerShotsToUse)+1)),3);
end
%% Prepare updated B0 TSC
for CurRep=1:nRepsToCalc
    CurB0VarUpdate=exp(1i.*2*pi*UpdatedB0Map_LLR_B0VarR(:,:,CurRep,SliI)*(1e-3).*TimePointsMed_ms3);
    CurB0VarUpdateP=permute(CurB0VarUpdate,[1 2 7 6 5 4 3]);
    CurTSCxPMedOnlyB0MS{SliI}(:,:,CurRep,1,1,1,:)=TSCxPMedOnlyB0.*CurB0VarUpdateP;
end
%% Rec several slices together, given updated B0
ScriptFN_CompgB0VarSomeShots=[BaseSP 'nuftCompgB0VarSomeShots_N.txt'];

Perm93=@(X) permute(X,[1 2 9 4 5 6 7 8 3]);
CurRep=[1 2 3];

Sz16CompgB0VarSomeShots=Sz16CompgB0;
%%
nShotsTogether=2;
RIdxs=1:nShotsTogether:nRepsToCalc;
for r=1:numel(RIdxs)
    CurRepa=RIdxs(r);
    CurRep=CurRepa+(0:(nShotsTogether-1));
    Rec_CompgB0_C_gB0Var_SomeShotsS{r,SliI}=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgB0VarSomeShots],...
        Sz16CompgB0VarSomeShots,Perm93(DataCCP(:,TrajPartMed,CurRep,1:nccToUse)),...
        SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed,Perm93(CurTSCxPMedOnlyB0MS{SliI}(:,:,CurRep,1,:,:,:,:,:)),...
        Perm93(KernsPMMed(:,:,CurRep,:,:,:,:)),CompsP);
end
%%
end % Slices loop
%%
for SliI=1:nSlices
    [Out B BN]=CalcSlicesSNR(max(grmss(THLRMultiShotS{SliI},3:30),eps),true,7);
    Msk=imfillholesBySlices(~BN);
    Msk=getLargestComponent(Msk);
    se = strel('disk', 3);
    DMsk=imdilate(Msk,se,'same');
    BS(:,:,SliI)=B;
    BNS(:,:,SliI)=BN;
    MskS(:,:,SliI)=Msk;
    DMskS(:,:,SliI)=DMsk;
end
BNS=(~BNS);
%%
for i=1:nSlices
    Rec_CompgB0_C_gB0Var_SomeShotsMS{i}=cat(5,Rec_CompgB0_C_gB0Var_SomeShotsS{:,i});
end
Rec_CompgB0_C_gB0Var_SomeShotsM=cat(4,Rec_CompgB0_C_gB0Var_SomeShotsMS{:});
Rec_CompgB0_C_gB0Var_SomeShotsM=squeeze(sum(Rec_CompgB0_C_gB0Var_SomeShotsM.*CompsP,6));

fgmontagex((squeeze(Rec_CompgB0_C_gB0Var_SomeShotsM(:,:,:,[1 2 5 10 15 20]))),[0 0.002]);
title('2shot, with per-shot B_0 variation, 3 different (2-)shots, several echos')

% ZZ=grmss(PartitionDim(Rec_CompgB0_C_gB0Var_SomeShotsM,5,4),5);
ZZ=grmss(PartitionDim(Rec_CompgB0_C_gB0Var_SomeShotsM,4,4),4);
fgmontagex(permute43(squeeze(ZZ(:,:,:,2:4))),[0 0.002]);
title('2-shot, 3 representing shots, 3 combined echos');

%%
XX=permute(Rec_CompgB0_C_gB0Var_SomeShotsM,[1 2 5 3 4]);
% XX=permute43(Rec_CompgB0_C_gB0Var_SomeShotsM);
WhichInnerShotsToUse=7:20;
[~,~,~,H_tmp]=ghankel(numel(WhichInnerShotsToUse),2,TrgSz);
for i=1:size(XX,4)
    for j=1:size(XX,5)
        CurX=XX(:,:,WhichInnerShotsToUse,i,j);
        [ ~, s_tmp, V_tmp] = batch_svd(H_tmp*CurX);
        R1a=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
        UpdatedT2SMap_LLRB0Var_SomsShots(:,:,i,j)=-InnerTSDiff_PerShot_ms./log(abs(R1a));
        UpdatedB0Map_LLRB0Var_SomsShots(:,:,i,j)=-(angle(R1a)/(2*pi))/(InnerTSDiff_PerShot_ms/1e3); % in Hz
        PD_LLRB0Var_SomsShots(:,:,i,j)=mean(CurX.*(R1a.^(-permute32(WhichInnerShotsToUse)+1)),3);
    end
end

fgmontagex(permute43(UpdatedT2SMap_LLRB0Var_SomsShots(:,:,ROrd,:)).*BNS(:,:,ROrd),[0 200]);
colormap hot;title('T_2^* map, 3 different 2-shots');

% fgmontagex(permute43(PartitionDim(UpdatedT2SMap_LLRB0Var_SomsShots(:,:,ROrd(1:9),2).*BNS(:,:,ROrd(1:9)),3,3)),[0 200]);colormap hot;title('T_2^* map, 2-shots');
fgmontagex(abs(permute43(PartitionDim(UpdatedT2SMap_LLRB0Var_SomsShots(:,:,ROrd(1:9),2).*BNS(:,:,ROrd(1:9)),3,3))),[0 200]);colormap hot;title('T_2^* map, 2-shots');
fgmontagex(abs(permute43(PartitionDim(PD_LLRB0Var_SomsShots(:,:,ROrd(1:9),2).*BNS(:,:,ROrd(1:9)),3,3))),[0 0.003]);title('PD map, 2-shots');

ZZ=grmss(PartitionDim(XX,3,4),3);
fgmontagex(permute43(squeeze(ZZ(:,:,ROrd(1:2:end),2,2:4).*BNS(:,:,ROrd(1:2:end)))),[0 0.0021]);title('2-shot, different echos');

fgmontagex(permute43(squeeze(ZZ(:,:,ROrd(1:2:end),:,2).*BNS(:,:,ROrd(1:2:end)))),[0 0.0021]);title('3 different 2-shot, same echo');

% fgmontagex(ZZ(:,:,ROrd([2 5 8]),2,2:4).*BNS(:,:,ROrd([2 5 8])),[0 0.0021]);title('3 different 2-shots, different echos');

fgmontagex(UpdatedT2SMap_LLRB0Var_SomsShots(:,:,:).*B,[0 200]);
colormap hot;title('T_2^* map, 3 different 2-shots');

fgmontagex(PD_LLRB0Var_SomsShots.*B,[0 0.003])
title('PD, 3 different 2-shots');

