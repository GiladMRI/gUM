%% Initialize X0
MskToUse=DMskS(:,:,SliI);

X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
X0=X0NoB.*TSCxPMedOnlyB0;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC,X0,...
        SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed,1,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*MskToUse.*Fac;
X0NoBN=X0NoB.*MskToUse.*Fac;
disp('Initialized X0');

AHA_Fac=grmss(AHAX0)./grmss(X0); % Poor effective Lipshitz factor

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;
%% THLR1-UV : initialize U,V
HankelOrder=3;
[~,~,~,H_AllTS]=ghankel(nTSMed,HankelOrder,TrgSz);
X0NoBN_TH=H_AllTS*squeeze(X0NoBN);
[ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
U=U_tmp(:,:,:,1).*sqrt(s_tmp(:,:,1));
V=V_tmp(:,:,:,1).*sqrt(s_tmp(:,:,1));
disp('Initialized U,V');
%% Iterate
X0new=squeeze(X0NoBN).*MskToUse;
alpha=0.1/AHA_Fac;
K=20;

iter=1;

for jj=1:5
    disp(jj);
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC,X0N,...
        SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed,1,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
rcur=AHy-AHAX0N;
rcur_rms(iter)=grmss(rcur);
disp(rcur_rms(iter));
CurRatioToDiff(iter)=AHy_rms./rcur_rms(iter);
rcurNoB=rcur.*conj(TSCxPMedOnlyB0);
UpdatedX0=X0new+alpha*squeeze(rcurNoB);
X0NoBN_TH=H_AllTS*squeeze(UpdatedX0);

for i=1:K
%     U=Mb*V/gsss(V);
    U=sum(X0NoBN_TH.*permute43(V),4)./gss(V,3);
%     V=(U'*Mb)'/gsss(U);
    V=permute43(sum(conj(U).*X0NoBN_TH,3)./gss(U,3));
end
U(isnan(U))=0;
V(isnan(V))=0;

% [ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
% s_tmpP=permute43(s_tmp);
% nHComps=3;
% U=U_tmp(:,:,:,1:nHComps).*(s_tmpP(:,:,1,1:nHComps));
% VH=permute43(conj(V_tmp(:,:,:,1:nHComps))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
% LREst=permute(sum(U.*permute(VH,[1 2 5 3 4]),4),[1 2 3 5 4]);
% LR1Est=U.*VH;
LR1Est=U.*conj(permute43(V));
Db=X0NoBN_TH-LR1Est;
disp(grmss(Db(~isnan(Db))))
Db_rms(iter)=grmss(Db(~isnan(Db)));
X0new=H_AllTS'*LR1Est;
X0N=permute(X0new,[1 2 7 6 5 4 3]).*TSCxPMedOnlyB0;

iter=iter+1;

R1var=V(:,:,2)./V(:,:,1);
% InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
UpdatedT2Sx_Map_ms=-InnerTSDiff_ms./log(abs(R1var));
UpdatedB0var=-(angle(R1var)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
end

figure;subplot(1,2,1);gmontage(UpdatedB0var,[-100 100]);title('\Delta B_0,T_2^* from TH-LR-1');removeTicks;colorbar
subplot(1,2,2);gmontage(abs(UpdatedT2Sx_Map_ms),[0 100]);removeTicks;