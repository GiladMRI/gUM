FirstEchoTE_ms=2.38;

WhichTSToUse=8:nTSMed4;
WhichTSToUseB=zeros(1,nTSMed4);
WhichTSToUseB(WhichTSToUse)=1;
WhichTSToUseB7=perm72(WhichTSToUseB);
% for CurReps=1:(nRepsHdr-1)
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}(:,:,:,:,:,1:3).*CompsP(:,:,:,:,:,1:3,:),6);
% % Using THLR-Multishot warmstart
% X0NoB=THLRMultiShot.*conj(TSCxPMedOnlyB0);
% X0NoB=Rec_LLR_B0Wraps_MSX.*conj(TSCxPMed2OnlyB0);

% CurRep=[1 2];
CurRep=[1];
nShotsTogether=numel(CurRep);
PerShotSmoothRegion=[40 40];
PerShotSmoothSig=5;
PerShotSmoothRegion=[20 20];
PerShotSmoothSig=5;

if(nShotsTogether==1)
%     X0a=squeeze(Rec_THLR2LX);
%     TSCEst=exp(1i*angle(Rec_THLR2LX));
    X0a=squeeze(Rec_CompgB0_C_gB0Var4_SomeShotsSMX(:,:,CurRep,:,:,:,:,:,:)).*squeeze(perm93(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,CurRep,:,:,:,:,:,:)));
    TSCEst=perm93(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,CurRep,:,:,:,:,:,:));
    X0aNoB=X0a.*conj(perm73(TSCEst));
    TSCPerShotCorr=1;
    X0aNoBxNoPerShot=perm73(X0aNoB);
else
    TSCEst=exp(1i.*angle(sum(CurTSCxPMed4OnlyB0AfterPerSlice,3)));
    X0a=Rec_CompgB0_Med4_C_WB0VarCSMX(:,:,:,CurRep);
    X0aNoB=X0a.*conj(perm73(TSCEst));
    for i=1:nShotsTogether
        [~, B0PerShota(:,:,i)]=getMPBDBySVDTH(X0aNoB(:,:,:,i),InnerMedTS4Diff_ms);
    end
    for i=1:nShotsTogether
        W=grmss(X0aNoB(:,:,:,i),3);
        B0W=B0PerShota(:,:,i).*W;
        B0WS=SmoothBySlices(B0W,PerShotSmoothRegion,PerShotSmoothSig);
        WSm=SmoothBySlices(grmss(tmp2,3),PerShotSmoothRegion,PerShotSmoothSig);
        B0Sm=B0WS./max(eps,WSm);
        B0PerShotaSm(:,:,i)=B0Sm;
    end
    TSCPerShotCorr=perm94(perm73(exp(1i.*2*pi*permute43(B0PerShotaSm)*(1e-3).*TimePointsMed4_ms3)));

    B0PerShotaBase=B0PerShota;

    % X0aNoBx=perm73(mean(X0aNoB.*conj(perm94(perm73(TSCPerShotCorr))),4));
    % X0aNoBx=perm94(perm73(X0aNoB.*conj(perm94(perm73(TSCPerShotCorr))))); % no main B, no per-shot B
    X0aNoBxNoPerShot=mean(perm94(perm73(X0aNoB)).*conj(TSCPerShotCorr),9);
end
X0aNoBx=repmat(X0aNoBxNoPerShot,[1 1 1 1 1 1 1 1 nShotsTogether]);

% X0NoB=permute73(mean(Rec_CompgB0_Med4_C_WB0VarCSMX(:,:,:,1:2).*conj(permute43(squeeze(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,1:2,:)))),4));
% dTSC21=CurTSCxPMed4OnlyB0AfterPerSlice(:,:,2,:,:,:,:,:,:,:,:,:).*conj(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,1,:,:,:,:,:,:,:,:,:));
% dTSC21x=cat(9,dTSC21*0+1,dTSC21);

% X0NoB=X0N.*conj(TSCxPMed2OnlyB0);
% X0NoB=X0N_MS.*conj(TSCxPMed2OnlyB0);
%
Sz16THLRSomeShot=FillOnesTo16(TrgSz);
Sz16THLRSomeShot(7)=nTSMed4;

Sz16THLRSomeShotx=Sz16THLRSomeShot;
Sz16THLRSomeShotx(9)=nShotsTogether;

ScriptFN_AllTS_gB0=[BaseSP 'nuftAllTS_gB0_N.txt'];

disp('Running AHy');
AHy=bart(['linopScript -d 0 -A ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,Perm93(DataCCP(:,TrajPartMed,CurRep,1:nccToUse)),...
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4.*WhichTSToUseB7,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:).*WhichTSToUseB7));

X0=X0aNoBx.*TSCEst;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,X0,...
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4.*WhichTSToUseB7,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:).*WhichTSToUseB7));

Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*Fac;
% X0N=X0aNoBx1.*Fac.*MskToUse;
% X0NoBN=X0NoB.*Fac;
X0NoBN=X0aNoBx.*Fac;
disp('Initialized X0');

% AHA_Fac=grmss(AHAX0)./grmss(X0); % Poor effective Lipshitz factor
rcur=AHy-AHAX0N;
AHAr=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,rcur,...
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4.*WhichTSToUseB7,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:).*WhichTSToUseB7));
AHA_Fac=grmss(AHAr)./grmss(rcur); % better Poor effective Lipshitz factor

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;

[~,~,~,H_AllTS4x]=ghankel(nTSMed4,HankelOrder,TrgSz);
[~,~,~,H_AllTS4_4]=ghankel(nTSMed4,4,TrgSz);
[~,~,~,H_AllTS4_3]=ghankel(nTSMed4,3,TrgSz);
[~,~,~,H_AllTS4_2]=ghankel(nTSMed4,2,TrgSz);

[~,~,~,H_AllTS4x_2]=ghankel(numel(WhichTSToUse),2,TrgSz);
[~,~,~,H_AllTS4x_3]=ghankel(numel(WhichTSToUse),3,TrgSz);

alpha=0.1/AHA_Fac;
alphaPerShotB=0.1/AHA_Fac;

% CadzowK=1;
nIters=30;
MinDiffT=3e-4;


DUse2ndSVMap=SmoothBySlices(Use2ndSVMap,[3 3],1.5)>(1/9);
% DUse2ndSVMap=DUse2ndSVMap*0>1;

UseOnly1CS=true;
iter=1;

if(nShotsTogether==1)
    X0N=perm73(FitToModel_1CSfx(squeeze(X0N(:,:,:,:,:,:,WhichTSToUse)),H_AllTS4x_2,nTSMed4));
end

disp('Ready to iterate');
%% Iterate
for jj=1:nIters
    disp(jj);
    
    X0Nr=repmat(X0N,[1 1 1 1 1 1 1 1 nShotsTogether]);
    AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,X0Nr,...
        SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4.*WhichTSToUseB7,TSCPerShotCorr,...
        Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:).*WhichTSToUseB7));
    CurFac=grmss(AHy)./grmss(AHAX0);
    AHAX0N=AHAX0*CurFac;
    rcur=AHy-AHAX0N;
    rcur_rms(iter)=grmss(rcur);
    disp(rcur_rms(iter));
    CurRatioToDiff(iter)=AHy_rms./rcur_rms(iter);
    
    UpdatedX0WPerShotWMainB=X0Nr+alpha*rcur; % with perShotB
    if(nShotsTogether==1)
        UpdatedX0=UpdatedX0WPerShotWMainB;
    else
        CommonPhaseAllEchos=exp(1i*angle(sum(UpdatedX0WPerShotWMainB,9)));
        UpdatedX0WPerShot=UpdatedX0WPerShotWMainB.*conj(CommonPhaseAllEchos);
    %     rcurNoB=rcur.*conj(TSCEst);
    %     UpdatedX0=X0new+alpha*squeeze(rcurNoB); % with perShotB, no mainB
    %     UpdatedX0WPerShot=X0new+alpha*rcurNoB; % with perShotB, no mainB

        for i=1:nShotsTogether
            [~, B0PerShota(:,:,i)]=getMPBDBySVDTH(squeeze(UpdatedX0WPerShot(:,:,:,:,:,:,:,:,i)),InnerMedTS4Diff_ms);
        end
        for i=1:nShotsTogether
            W=grmss(UpdatedX0WPerShot(:,:,:,:,:,:,:,:,i),7);
            B0W=B0PerShota(:,:,i).*W;
            B0WS=SmoothBySlices(B0W,PerShotSmoothRegion,PerShotSmoothSig);
            WSm=SmoothBySlices(grmss(tmp2,3),PerShotSmoothRegion,PerShotSmoothSig);
            B0Sm=B0WS./max(eps,WSm);
            B0PerShotaSm(:,:,i)=B0Sm;
        end
        TSCPerShotCorr=perm94(perm73(exp(1i.*2*pi*permute43(B0PerShotaSm)*(1e-3).*TimePointsMed4_ms3)));
    %     fgmontage(B0PerShotaSm,[-4 4])
        UpdatedX0=mean(UpdatedX0WPerShot.*conj(TSCPerShotCorr),9).*CommonPhaseAllEchos;
    end
%     Fitted=FitToModel_1CSf(squeeze(UpdatedX0),H_AllTS4_2);
%     Fitted=FitToModel_1CSf(squeeze(UpdatedX0(:,:,:,:,:,:,WhichTSToUse)),H_AllTS4x_2);
    Fitted=FitToModel_1CSfx(squeeze(UpdatedX0(:,:,:,:,:,:,WhichTSToUse)),H_AllTS4x_2,nTSMed4);
    
%     Fitted=FitToModel_1CSf(squeeze(UpdatedX0),H_AllTS4_2);
    if(~UseOnly1CS)
        Fitted2=FitToModel_2CSf(squeeze(UpdatedX0(:,:,:,:,:,:,WhichTSToUse)),H_AllTS4x_3);
        Fitted=Fitted.*(1-DUse2ndSVMap)+Fitted2.*DUse2ndSVMap;
    end

    X0new=perm73(Fitted).*MskToUse;
    LastX0N=X0N;
    X0N=X0new;
        
    dX0N(iter)=grmss(LastX0N-X0N)./grmss(LastX0N);
    disp(['dX0N(iter) : ' num2str(dX0N(iter))]);
    if(dX0N(iter)<MinDiffT)
        disp('Less the MinDiffT, breaking');
        break;
    end
    
    iter=iter+1;
end
disp('Finished iterations');
%%
MPBD_2CS2Shot=getMPBDBySVDTH(X0N,InnerMedTS4Diff_ms,FirstEchoTE_ms);
fgmontagex(permute43(PartitionDim(MPBD_2CS2Shot,3,2)),[0 1])

% ShowAbsAngle(X0N(:,:,:,:,:,:,[1 23 50 75]))
ShowAbsAngle(X0N(:,:,:,:,:,:,[1 17 33 49]))
%%
MPBD_2CS2Shot=getMPBDBySVDTH(X0Nbase,InnerMedTS4Diff_ms,FirstEchoTE_ms);TtlX='base';
fgmontagex(permute43(PartitionDim(MPBD_2CS2Shot,3,2)),[0 1]);title(TtlX)

% ShowAbsAngle(X0N(:,:,:,:,:,:,[1 23 50 75]))
ShowAbsAngle(X0Nbase(:,:,:,:,:,:,[1 17 33 49]));title(TtlX)
%%
XX=PartitionDim(squeeze(X0N(:,:,:,:,:,:,1:80)),3,4);
fgmontagex(grmss(XX,3))

B0sR=(MPBDR(:,:,1,:)-0.5)*200;

X0NRC{min(CurReps)}=X0N;