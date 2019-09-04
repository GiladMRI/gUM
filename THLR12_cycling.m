CurReps=1;

AHy=bart(['linopScript -A ' ScriptFN_AllTS],Sz16AllTSC,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
    SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed,1,...
    sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
%%
MskToUse=DMskS(:,:,SliI);

HankelOrder=3;
%%
THLRMultiShot=THLRMultiShotS{SliI};

[~,~,~,H_AllTS]=ghankel(nTSMed,HankelOrder,TrgSz);
[ ~, s_THLR_Multishot, V_THLR_Multishot] = batch_svd(H_AllTS*squeeze(THLRMultiShot));
% R1ts=V_THLR_Multishot(:,:,2,1)./V_THLR_Multishot(:,:,1,1); % R1 is simply the decay
% InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
% UpdatedT2SMap_ms1=-InnerTSDiff_ms./log(abs(R1ts));
% UpdatedB0Map1=-(angle(R1ts)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz

sRatio=(s_THLR_Multishot(:,:,2)./s_THLR_Multishot(:,:,1)).*MskToUse.*BNS(:,:,SliI);
Use2ndSVMap=sRatio>0.1;
disp('Got Use2ndSVMap');
%%
figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
%% Initialize X0
X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
% Using THLR-Multishot warmstart
X0NoB=THLRMultiShot.*conj(TSCxPMedOnlyB0);

X0=X0NoB.*TSCxPMedOnlyB0;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC,X0,...
    SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed,1,...
    sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*Fac;
X0NoBN=X0NoB.*Fac;
disp('Initialized X0');

% AHA_Fac=grmss(AHAX0)./grmss(X0); % Poor effective Lipshitz factor
rcur=AHy-AHAX0N;
AHAr=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC,rcur,...
    SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed,1,...
    sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
AHA_Fac=grmss(AHAr)./grmss(rcur); % better Poor effective Lipshitz factor

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;

[~,~,~,H_AllTS]=ghankel(nTSMed,HankelOrder,TrgSz);

iter=1;

X0new=squeeze(X0NoBN).*MskToUse;
alpha=0.01/AHA_Fac;

Wlambda_base=0.1;
Wlambda=0.0003;
WT=wave_thresh('db4', 3, Wlambda_base);

disp('Ready to iterate');
%% Iterate
for jj=1:30
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
    
    [ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
    s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);
%     [~,~,~,H_AllTS2]=ghankel(nTSMed,2,TrgSz);
%     [ ~, s_tmp, V_tmp] = batch_svd(H_AllTS2*squeeze(XX));
%     R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
%     PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed-1)))));
%     W=abs(R1_tmp).^permute32(0:(nTSMed-1));
%     PDEstx=sum(PDMEx,3)./sum(W,3);
%     PDEstx_W=WT(PDEstx,Wlambda);
%     XXWT=PDEstx_W.*(conj(R1_tmp).^permute32(0:(nTSMed-1)));
%     X0N=permute(XXWT,[1 2 7 6 5 4 3]);

    % U=U_tmp(:,:,:,1:2).*sqrt(s_tmpP(:,:,1,1:2));
    % V=V_tmp(:,:,:,1:2).*sqrt(s_tmpP(:,:,1,1:2));
    U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
    VH1=permute43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
    
    U2=U_tmp(:,:,:,2).*(s_tmpP(:,:,1,2));
    VH2=permute43(conj(V_tmp(:,:,:,2))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
    
    LR1Est=permute(sum(U1.*permute(VH1,[1 2 5 3 4]),4),[1 2 3 5 4]);
    LR2Est=permute(sum(U2.*permute(VH2,[1 2 5 3 4]),4),[1 2 3 5 4]);
    LREst=LR1Est+LR2Est.*Use2ndSVMap;
    
    Db=X0NoBN_TH-LREst;
    Db_rms(iter)=grmss(Db(~isnan(Db)));
    X0new=(H_AllTS'*LREst).*MskToUse;
    LastX0N=X0N;
    X0N=permute(X0new,[1 2 7 6 5 4 3]).*TSCxPMedOnlyB0;
    
    % Do wavelet thresholding to PD
%     XX=permute(X0N,[1 2 7 6 5 4 3]);
%     [~,~,~,H_AllTS2]=ghankel(nTSMed,2,TrgSz);
%     [ ~, s_tmp, V_tmp] = batch_svd(H_AllTS2*squeeze(XX));
%     R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
%     PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed-1)))));
%     W=abs(R1_tmp).^permute32(0:(nTSMed-1));
%     PDEstx=sum(PDMEx,3)./sum(W,3);
%     PDEstx_W=WT(PDEstx,Wlambda);
%     XXWT=PDEstx_W.*(conj(R1_tmp).^permute32(0:(nTSMed-1)));
%     X0N=permute(XXWT,[1 2 7 6 5 4 3]);
    %
    
    dX0N(iter)=grmss(LastX0N-X0N)./grmss(LastX0N);
    disp(['dX0N(iter) : ' num2str(dX0N(iter))]);
    
    iter=iter+1;
    
%     R1var=V(:,:,2)./V(:,:,1);
%     % InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
%     UpdatedT2Sx_Map_ms=-InnerTSDiff_ms./log(abs(R1var));
%     UpdatedB0var=-(angle(R1var)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
%     
%     sRatiox=(s_tmp(:,:,2)./s_tmp(:,:,1)).*MskToUse;
end

% fgmontagex(permute43(PartitionDim(squeeze(X0N),3,4)),[0 2.4e-3])
disp('Finished iterations');
%%
X0N_NoWS=X0N;
X0N_WS=X0N;
%%
X0N_NoWSC{SliI}=X0N_NoWS;
X0N_WSC{SliI}=X0N_WS;
THLRMultiShot=THLRMultiShotS{SliI};
THLRMultiShotN=THLRMultiShot.*grmss(X0N_NoWSC{SliI})./grmss(THLRMultiShot);
YY=cat(4,squeeze(THLRMultiShotN),squeeze(X0N_WSC{SliI}),squeeze(X0N_NoWSC{SliI}));
fgmontagex(permute43(YY(:,:,1:5:end,:)),[0 2.4e-3])
ylabel('Single-shot given B_0   |    single-shot with multi-shot warmstart   |   39-shot','FontSize',16);

% [~,~,~,H_AllTS2]=ghankel(nTSMed,2,TrgSz);
for i=1:3
    XX=permute(YY(:,:,:,i),[1 2 7 6 5 4 3]);
    [ ~, s_tmp, V_tmp] = batch_svd(H_AllTS2*squeeze(XX));
    R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
    InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
    T2SMap_tmp_ms=-InnerTSDiff_ms./log(abs(R1_tmp));
    B0Map_tmp=-(angle(R1_tmp)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz
    PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed-1)))));
    W=abs(R1_tmp).^permute32(0:(nTSMed-1));
    PDEstx=sum(PDMEx,3)./sum(W,3);
    
    MPBD(:,:,1,i)=abs(PDEstx)/(20e-4);
    MPBD(:,:,2,i)=(angle(PDEstx)+pi)/(2*pi);
    MPBD(:,:,3,i)=abs(T2SMap_tmp_ms)/100;
    MPBD(:,:,4,i)=(B0Map_tmp+100)/200;
end
MPBD(:,:,1,:)=MPBD(:,:,1,:)*1.3;
fgmontagex(permute43(MPBD),[0 1]);
title('PD magnitude, phase                                         T_2^*                        B_0','FontSize',16);
ylabel('Single-shot given B_0   |    single-shot with multi-shot warmstart   |   39-shot','FontSize',16);
%%
% XX=squeeze(THLRMultiShot);TtlX='THLRMultiShot';XX=XX.*grmss(X0N_WS)./grmss(XX);
% XX=X0N_WS;TtlX='THLR12 + warmstart';
% XX=X0N_NoWS;TtlX='THLR12';
XX=X0N;TtlX='THLR12 + WT';
[~,~,~,H_AllTS2]=ghankel(nTSMed,2,TrgSz);
[ ~, s_tmp, V_tmp] = batch_svd(H_AllTS2*squeeze(XX));
R1_tmp=V_tmp(:,:,2,1)./V_tmp(:,:,1,1); % R1 is simply the decay
InnerTSDiff_ms=numel(TrajPartMed)*AcqDwellTime_us/1e3/(nTSMed-1);
T2SMap_tmp_ms=-InnerTSDiff_ms./log(abs(R1_tmp));
B0Map_tmp=-(angle(R1_tmp)/(2*pi))/(InnerTSDiff_ms/1e3); % in Hz

PDMEx=squeeze(XX).*(exp(1i*angle(R1_tmp).*permute32((0:(nTSMed-1)))));
W=abs(R1_tmp).^permute32(0:(nTSMed-1));
PDEstx=sum(PDMEx,3)./sum(W,3);

figure;
subplot(2,2,1);gmontage(B0Map_tmp,[-100 100]);title(['B_0,T_2^* from ' TtlX]);removeTicks;colorbar
subplot(2,2,2);gmontage(abs(T2SMap_tmp_ms),[0 100]);removeTicks;
subplot(2,2,3);gmontage(abs(PDEstx),[0 2.4e-3]);title('PD');removeTicks
subplot(2,2,4);gmontage(angle(PDEstx),[-pi pi]);removeTicks;
%%
% WT=wave_thresh('db4', 3, 0.1);
% PDEstx_W=WT(PDEstx,0.001);

% XXWT=PDEstx_W.*(conj(R1_tmp).^permute32(0:(nTSMed-1)));