% for CurReps=1:(nRepsHdr-1)
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
% X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}(:,:,:,:,:,1:3).*CompsP(:,:,:,:,:,1:3,:),6);
% % Using THLR-Multishot warmstart
% X0NoB=THLRMultiShot.*conj(TSCxPMedOnlyB0);
% X0NoB=Rec_LLR_B0Wraps_MSX.*conj(TSCxPMed2OnlyB0);

CurRep=[1 2];
nShotsTogether=numel(CurRep);
PerShotSmoothRegion=[40 40];
PerShotSmoothSig=5;

TSCEst=exp(1i.*angle(sum(CurTSCxPMed4OnlyB0AfterPerSlice,3)));
X0a=Rec_CompgB0_Med4_C_WB0VarCSMX(:,:,:,CurRep);
X0aNoB=X0a.*conj(permute73(TSCEst));
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
X0aNoBx=mean(perm94(perm73(X0aNoB)).*conj(TSCPerShotCorr),9);
X0aNoBx=repmat(X0aNoBx,[1 1 1 1 1 1 1 1 nShotsTogether]);

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
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));

% AHy=bart(['linopScript -d 0 -A ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,Perm93(DataCCP(:,TrajPartMed,CurRep,1:nccToUse)),...
%     SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,dTSC21x,...
%     Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));

% AHy=bart(['linopScript -d 0 -A ' ScriptFN_AllTS_gB0],Sz16THLRSomeShot,Perm93(DataCCP(:,TrajPartMed,CurRep,1:nccToUse)),...
%     SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,dTSC21x,...
%     Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));
    
% X0=X0NoB.*CurTSCxPMed4OnlyB0AfterPerSlice(:,:,1,:,:,:,:,:,:,:,:,:);
% X0=X0NoB.*TSCEst;
% AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShot,X0,...
%     SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,dTSC21x,...
%     Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));

X0=X0aNoBx.*TSCEst;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,X0,...
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));

Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*Fac;
% X0NoBN=X0NoB.*Fac;
X0NoBN=X0aNoBx.*Fac;
disp('Initialized X0');

% AHA_Fac=grmss(AHAX0)./grmss(X0); % Poor effective Lipshitz factor
rcur=AHy-AHAX0N;
AHAr=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,rcur,...
    SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,TSCPerShotCorr,...
    Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));
AHA_Fac=grmss(AHAr)./grmss(rcur); % better Poor effective Lipshitz factor

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;

UpdatedX0=mean(UpdatedX0WPerShot.*conj(TSCPerShotCorr),9);

[~,~,~,H_AllTS4x]=ghankel(nTSMed4,HankelOrder,TrgSz);
[~,~,~,H_AllTS4_4]=ghankel(nTSMed4,4,TrgSz);
[~,~,~,H_AllTS4_2]=ghankel(nTSMed4,2,TrgSz);

iter=1;

% X0new=squeeze(X0NoBN).*MskToUse;
X0new=X0NoBN.*MskToUse;
alpha=0.1/AHA_Fac;
alphaPerShotB=0.1/AHA_Fac;

nIters=30;
MinDiffT=3e-4;

CadzowK=1;

% DUse2ndSVMap=SmoothBySlices(Use2ndSVMap,[3 3],1.5)>(1/9);
DUse2ndSVMap=DUse2ndSVMap*0>1;

disp('Ready to iterate');
%% Iterate
for jj=1:nIters
    disp(jj);
    
    AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShotx,repmat(X0N,[1 1 1 1 1 1 1 1 2]),...
        SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,TSCPerShotCorr,...
        Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));
    CurFac=grmss(AHy)./grmss(AHAX0);
    AHAX0N=AHAX0*CurFac;
    rcur=AHy-AHAX0N;
    rcur_rms(iter)=grmss(rcur);
    disp(rcur_rms(iter));
    CurRatioToDiff(iter)=AHy_rms./rcur_rms(iter);
    rcurNoB=rcur.*conj(TSCEst);
%     UpdatedX0=X0new+alpha*squeeze(rcurNoB); % with perShotB, no mainB
    UpdatedX0WPerShot=X0new+alpha*rcurNoB; % with perShotB, no mainB
    
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

    UpdatedX0=mean(UpdatedX0WPerShot.*conj(TSCPerShotCorr),9);

%     tmpttt=repmat(X0new,[1 1 1 2])+alpha*squeeze(rcurNoB);
%     tmp2=mean(abs(tmpttt),4).*exp(1i*(angle(tmpttt(:,:,:,2))-angle(tmpttt(:,:,:,1))));
%     [ ~, s_tmp, V_tmp] = batch_svd(H_AllTS4_2*tmp2);
%     R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
%     B0_Med4_tmp=-(angle(R1_tmp)/(2*pi))/(InnerMedTS4Diff_ms/1e3); % in Hz
%     
%     W=grmss(tmp2,3);
%     B0W=B0_Med4_tmp.*W;
%     B0WS=SmoothBySlices(B0W,[40 40],10);
%     WSm=SmoothBySlices(grmss(tmp2,3),[20 20],5);
%     B0Sm=B0WS./max(eps,WSm);
%     fgmontage(B0Sm,[-5 5])

    
%     AHAX0=bart(['linopScript -N ' ScriptFN_AllTS_gB0],Sz16THLRSomeShot,X0N,...
%         SensCC(:,:,:,1:nccToUse),Perm93(STraj3MMed(:,:,CurRep)),TSBPMed4,dTSC21x,...
%         Perm93(KernsPMMed4(:,:,CurRep,:,:,:,:)));
%     
%     
%     
%     CurFac=grmss(AHy)./grmss(AHAX0);
%     AHAX0N=AHAX0*CurFac;
%     rcur=AHy-AHAX0N;
%     rcur_rms(iter)=grmss(rcur);
%     disp(rcur_rms(iter));
%     CurRatioToDiff(iter)=AHy_rms./rcur_rms(iter);

%     rcurNoB=rcur.*conj(CurTSCxPMed4OnlyB0AfterPerSlice(:,:,1,:,:,:,:,:,:,:,:,:));
%     UpdatedX0=X0new+alpha*squeeze(rcurNoB);
    
    Fitted=UpdatedX0;
    for k=1:CadzowK
        TH=H_AllTS4_2*squeeze(Fitted);
        [ U_tmp, s_tmp, V_tmp] = batch_svd(TH);
        s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);
        U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
        U1(:,:,2:end)=min(abs(U1(:,:,2:end)),abs(U1(:,:,1:end-1))).*exp(1i*angle(U1(:,:,2:end)));
        VH1=permute43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
        VH1(:,:,1,2)=min(abs(VH1(:,:,1,2)),abs(VH1(:,:,1,1))).*exp(1i*angle(VH1(:,:,1,2)));
        LR1Est=permute(sum(U1.*permute(VH1,[1 2 5 3 4]),4),[1 2 3 5 4]);
        Fitted=(H_AllTS4_2'*LR1Est);
    end    
    
%     Fitted=FitToModel_THLR2Mskf(UpdatedX0,H_AllTS4x,DUse2ndSVMap);
    
%     [ U_tmp, s_tmp0, V_tmp] = batch_svd(H_AllTS4x*UpdatedX0);
%     Fitted=FitToModel_THLR2Mskf(UpdatedX0,H_AllTS4x,DUse2ndSVMap);
%     [ U_tmp, s_tmp, V_tmp] = batch_svd(H_AllTS4x*Fitted);
%     for k=2:CadzowK
%         Fitted=FitToModel_THLR2Mskf(Fitted,H_AllTS4x,DUse2ndSVMap);
%     end
%     [ U_tmp, s_tmp2, V_tmp] = batch_svd(H_AllTS4x*Fitted);
%     Fitted=FitToModel_THLR2Mskf(Fitted,H_AllTS4x,DUse2ndSVMap);
%     [ U_tmp, s_tmp3, V_tmp] = batch_svd(H_AllTS4x*Fitted);
%     s_tmpx=cat(4,s_tmp0,s_tmp,s_tmp2,s_tmp3);
    
%     In=UpdatedX0;
%     H=H_AllTS4x;
%     X0NoBN_TH=H*squeeze(In);
% 
% [ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
% s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);
% 
% U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
% VH1=permute43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
% 
% U2=U_tmp(:,:,:,2).*(s_tmpP(:,:,1,2));
% VH2=permute43(conj(V_tmp(:,:,:,2))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
% 
% % v1: [1 a ab]
% % v2: [1 c cd]
% % have a balanced linear combination, w, 1-w:
% % aw+c(1-w), abw+cd(1-w)
% % (aw+c(1-w))^2=abw+cd(1-w)
% 
% % syms a b c d w
% % eqn=(a*w+c*(1-w))^2==a*b*w+c*d*(1-w);
% % sol=solve(eqn,w);
% % qs=a*sol+c*(1-sol);
% % simplify(qs)
% % w1=  ((a^2*b^2 - 4*a^2*b*c + 4*a^2*c*d + 4*a*b*c^2 - 2*a*b*c*d - 4*a*c^2*d + c^2*d^2)^(1/2) + a*b - 2*a*c - c*d + 2*c^2)/(2*(a - c)^2)
% % w2= -((a^2*b^2 - 4*a^2*b*c + 4*a^2*c*d + 4*a*b*c^2 - 2*a*b*c*d - 4*a*c^2*d + c^2*d^2)^(1/2) - a*b + 2*a*c + c*d - 2*c^2)/(2*(a - c)^2)
% 
% LR1Est=permute(sum(U1.*permute(VH1,[1 2 5 3 4]),4),[1 2 3 5 4]);
% LR2Est=permute(sum(U2.*permute(VH2,[1 2 5 3 4]),4),[1 2 3 5 4]);
% LREst=LR1Est+LR2Est.*Use2ndSVMap;
% 
% Db=X0NoBN_TH-LREst;
% Db_rms=grmss(Db(~isnan(Db)));
%     
% Out=(H'*LREst);
% 
% Fitted=Out;


    X0new=perm73(Fitted).*MskToUse;
    LastX0N=X0N;
%     X0N=permute73(X0new).*CurTSCxPMed4OnlyB0AfterPerSlice(:,:,1,:,:,:,:,:,:,:,:,:);
%     X0N=permute73(X0new).*TSCEst;
    X0N=X0new.*TSCEst;
        
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
MPBD_2CS2Shot=getMPBDBySVDTH(X0N,InnerMedTS4Diff_ms);
fgmontagex(permute43(PartitionDim(MPBD_2CS2Shot,3,2)),[0 1])

XX=PartitionDim(squeeze(X0N(:,:,:,:,:,:,1:80)),3,4);
fgmontagex(grmss(XX,3))

B0sR=(MPBDR(:,:,1,:)-0.5)*200;

X0NRC{min(CurReps)}=X0N;