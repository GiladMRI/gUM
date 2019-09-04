SliI=10;

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
%%
figure;subplot(1,2,1);gmontage(UpdatedB0Map1,[-100 100]);title('B_0,T_2^* from TH-LR multi-shot');removeTicks;colorbar
subplot(1,2,2);gmontage(UpdatedT2SMap_ms1,[0 100]);removeTicks;
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

%%
ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];

CurReps=1;

AHy=bart(['linopScript -A ' ScriptFN_AllTS],Sz16AllTSC,DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSBPMed,1,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));

AHy_rms=grmss(AHy);
%% Initialize X0
MskToUse=DMskS(:,:,SliI);

X0NoB=sum(Rec_CompgB0_C_WB0VarCS{1,SliI}.*CompsP,6);
X0=X0NoB.*TSCxPMedOnlyB0;
AHAX0=bart(['linopScript -N ' ScriptFN_AllTS],Sz16AllTSC,X0,...
        SensCC(:,:,:,1:nccToUse).*MskToUse,STraj3MMed(:,:,CurReps),TSBPMed,1,...
        sum(KernsPMMed(:,:,CurReps,:,:,:,:),3));
Fac=grmss(AHy)./grmss(AHAX0);
AHAX0N=AHAX0*Fac;
X0N=X0.*Fac;
X0NoBN=X0NoB.*Fac;
disp('Initialized X0');

X0Nbase=X0N.*MskToUse;
X0NoBNbase=X0NoBN.*MskToUse;
%%
rcur=AHy-AHAX0N;
CurRatioToDiff=grmss(AHy)./grmss(AHy-AHAX0N);

rcurNoB=rcur.*conj(TSCxPMedOnlyB0);
rcurNoBM=rcurNoB.*abs(X0N); % ??? or use PD?
%%
PMat=[ones(1,nTSMed); 0:(nTSMed-1)];
PMat34=permute(PMat,[3 4 2 1]);
%%
PHrcurNoBM=permute43(sum(PMat34.*squeeze(rcurNoBM),3));
fgmontage(cat(3,real(PHrcurNoBM),imag(PHrcurNoBM)))
fgmontage(cat(3,real(PHrcurNoBM(:,:,2)),imag(PHrcurNoBM(:,:,2))))
UpdatedB0Map_TSUnits=-(angle(R1ts));
% Now real, or imag, or -real,imag?
% alpha
% add to Pcur
NewB0Map_TSUnits=UpdatedB0Map_TSUnits-10000*imag(PHrcurNoBM(:,:,2));
% Proximal
% Use now B0
%% THLR1-UV : initialize U,V
[~,~,~,H_AllTS]=ghankel(nTSMed,2,TrgSz);
X0NoBN_TH=H_AllTS*squeeze(X0NoBN.*MskToUse);
[ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
U=U_tmp(:,:,:,1).*sqrt(s_tmp(:,:,1));
V=V_tmp(:,:,:,1).*sqrt(s_tmp(:,:,1));
disp('Initialized U,V');
%% small test
% D=X0NoBN_TH-U.*conj(permute43(V));
% K=20;
% for i=1:K
% %     U=Mb*V/gsss(V);
%     U=sum(X0NoBN_TH.*permute43(V),4)./gss(V,3);
% %     V=(U'*Mb)'/gsss(U);
%     V=permute43(sum(conj(U).*X0NoBN_TH,3)./gss(U,3));
% end
% LR1Est=U.*conj(permute43(V));
% Db=X0NoBN_TH-LR1Est;
% X0new=H_AllTS'*LR1Est;
% X0N=permute(X0new,[1 2 7 6 5 4 3]).*TSCxPMedOnlyB0;
%%
X0new=squeeze(X0NoBN).*MskToUse;
alpha=0.03;
K=20;

iter=1;

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

for i=1:K
%     U=Mb*V/gsss(V);
    U=sum(X0NoBN_TH.*permute43(V),4)./gss(V,3);
%     V=(U'*Mb)'/gsss(U);
    V=permute43(sum(conj(U).*X0NoBN_TH,3)./gss(U,3));
end
U(isnan(U))=0;
V(isnan(V))=0;
LR1Est=U.*conj(permute43(V));
Db=X0NoBN_TH-LR1Est;
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
subplot(1,2,2);gmontage(UpdatedT2Sx_Map_ms,[0 100]);removeTicks;

%% Simulation
GT_Components
GT_X0
GT_B0
Delta_B0
Start_B0
cycle