% TEs8: 5.57 ms
% 8.67 ms
% 11.77 ms
% 14.87 ms
% 17.97 ms
% 21.07 ms
% 24.17 ms
% 27.27 ms


% RO is 1280
% TR 34ms
% TEs 5.65,11.95,18.25,24.55
% 576 slices
% Thickness 0.150mm
% FOV RO 192mm
% Phase 81.3%
% Phase resolution 100%
% Dist factor 20%?

% RO BW 180 Hz/Px
% EchoSpacing 6.3ms
%%
addpath('/autofs/cluster/exvivo/rock/4/users/hires/scripts/offline_recon/wrapper_scripts/matlab_code');
[prot, evp] = read_meas_prot('/autofs/space/minar_001/users/divya/forGilad/meas_MID200_GRE_4e_150um_FA10_seg0of7_[start+stop]_FID398.hdr');

mr_parms = mrir_measdat_parms(prot);
M0_vox2ras = mrir_measdat_vox2ras(prot, evp);

BaseP='/autofs/cluster/kawin/Gilad/FromDivya/FA10_ROneg/';
BaseP='/autofs/space/oribi_001/users/streaming/cases/I53_lh_32chexvivo_20190316/mri/FA10_ROneg/';
BaseP='/autofs/space/oribi_001/users/streaming/cases/I53_lh_32chexvivo_20190316/mri/FA20_ROneg/';
BaseP='/autofs/space/oribi_001/users/streaming/cases/I53_lh_32chexvivo_20190316/mri/FA10/';
BaseP='/autofs/space/oribi_001/users/streaming/cases/I53_lh_32chexvivo_20190316/mri/FA40/';
echoI=0;
chI=7;
BasePX='/autofs/cluster/kawin/Gilad/FromDivya/';
%%
XSliC=zeros(1280,1040,31,4);
%%
for echoI=0:3
    for chI=0:30
        disp([echoI chI]);
        disp(datestr(now));
        FN=[BaseP 'echo_' num2str(echoI) '_channel_' num2str(chI) '.complex_image'];
        
        fp = fopen(FN);
%         fseek(fid,2*,'bof');
        Xvec = single(fread(fp, 'float=>float'));
        fclose(fp);
        
        Xvol=reshape(Xvec, 2, evp.NImageCols, evp.NImagePar, evp.NImageLins);
        
        SliIToUse=250;
        XSli=squeeze(Xvol(1,:,SliIToUse,:)+1i*Xvol(2,:,SliIToUse,:));
        
        XSliC(:,:,chI+1,echoI+1)=XSli;
    end
end
disp('Saving...');
save('/autofs/cluster/kawin/Gilad/FromDivya/XSliC_FA40.mat','XSliC','-v7.3');
% save('XSliC_FA20_ROneg.mat','XSliC','-v7.3');
disp('Saved');
%%
% save('XSliC_FA10_ROneg.mat','XSliC','-v7.3');
load('/autofs/cluster/kawin/Gilad/FromDivya/XSliC_FA40.mat','XSliC');
%%
% XSliC02=XSliC(:,:,:,[1 3]);
% XSliC13=XSliC(:,:,:,[2 4]);
% % save('XSliC02.mat','XSliC02');
% % save('XSliC13.mat','XSliC13');
% %%
% load('XSliC13.mat','XSliC13');
% XSliC(:,:,:,[2 4])=XSliC13;
%% Renormalize
XSliCN=XSliC./grms(XSliC(1:50,1:50,:,1:4),[1:2 4]);
rXSliCN=grmss(XSliCN,3);
rXSliCNx=rXSliCN/grms(rXSliCN,1:3);
%% Good region
fgmontagex(grmss(XSliCN(600:end-100,200:end-200,:,:),3))
%%
RGB=grmss(XSliCN(800:end-100,200:end-200,:,1:3),3);

RGB=grmss(XSliCN(1:end,200:end-200,:,1:3),3);
%%
E3=rXSliCN(:,:,3);
E3x=sqrt(rXSliCN(:,:,2).*rXSliCN(:,:,4));
%%
load([BasePX 'XSliC_FA10.mat']);
XSliCp=XSliC;
%%
load([BasePX 'XSliC_FA10_ROneg.mat']);
XSliCn=XSliC;
%%
XSliCpN=XSliCp./grms(XSliCp(1:50,1:50,:,1:4),[1:2 4]);
rXSliCpN=grmss(XSliCpN,3);

XSliCnN=XSliCn./grms(XSliCp(1:50,1:50,:,1:4),[1:2 4]);
rXSliCnN=grmss(XSliCnN,3);


fgmontagex(rXSliCpN(600:end-100,200:end-200,3))
fgmontagex(rXSliCnN(600:end-100,200:end-200,3))
%%
fgmontagex(rXSliCpN(600:end-100,200:end-200,1))
fgmontagex(rXSliCpN(600:end-100,200:end-200,2))
fgmontagex(rXSliCpN(600:end-100,200:end-200,3))
fgmontagex(rXSliCpN(600:end-100,200:end-200,4))
%%
for i=1:4
    [~,threshold] = edge(rXSliCpN(600:end-100,200:end-200,i),'sobel');
    BWs(:,:,i) = edge(rXSliCpN(600:end-100,200:end-200,i),'sobel',threshold * 1);
end
%%
BWsX=sum(BWs.*perm32([1 2 4 8]),3);
figure;imagesc(BWsX);colormap colorcube
%%
E3=rXSliCpN(:,:,3);
E3x=sqrt(rXSliCpN(:,:,2).*rXSliCpN(:,:,4));
dE3=E3-E3x;
dE3x=cat(3,E3,E3x,dE3*10);

fgmontagex(abs(dE3x(800:end-100,500:900,:)),'Size',[1 3]);title('3rd Echo,   Simulated 3rd echo,    difference*10');
%%
CurSens=RunESPIRiTForSensMapsMultiMap((XSliCpN(:,:,:,1)),0);
CurSens1=CurSens(:,:,:,1);
CurSens1=single(CurSens1);
save([BasePX 'Sens.mat'],'CurSens1');
% save([BasePX 'Sens40.mat'],'CurSens1');
%%
load([BasePX 'Sens.mat'],'CurSens1');
%%
SliCpNC=squeeze(sum(XSliCpN.*conj(CurSens1),3));
%%
ES_ms=6.3;
FirstTE_ms=5.65;
WhichEchosToUseA=1:4;
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(SliCpNC,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
E3c=SliCpNC(:,:,3);
sE3c=sqrt(SliCpNC(:,:,2).*SliCpNC(:,:,4));
BE3c=cat(3,E3c,sE3c);
%%
fgmontagex(UpdatedB0Map_Hz,[-1 1]*(1000/(2*ES_ms)));colorbar
%%
[PDBasea, UpdatedB0Map_Hza, UpdatedT2SMap_msa, s_valsa, Fitted0a, PDBase0a]=...
    FitToModel_MPBD1CSf(abs(SliCpNC),WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
P1=angle(SliCpNC(:,:,3).*conj(SliCpNC(:,:,1)));
P2=angle(SliCpNC(:,:,4).*conj(SliCpNC(:,:,2)));

dP=angle(exp(1i*(P2-P1)));
dPx=dP;
dPx((dP>1) | (dP<0.3))=NaN;
%%
dAngle=UpdatedB0Map_Hz*2*pi*ES_ms/1000;
P1u = cusackUnwrap(repmat(P1,[1 1 3]), repmat(grmss(SliCpNC(:,:,[1 3]),3),[1 1 3]));
P1u=P1u(:,:,2);

P2u = cusackUnwrap(repmat(P2,[1 1 3]), repmat(grmss(SliCpNC(:,:,[2 4]),3),[1 1 3]));
P2u=P2u(:,:,2);
%%
C1=SliCpNC(:,:,3).*conj(SliCpNC(:,:,1));
C2=SliCpNC(:,:,4).*conj(SliCpNC(:,:,2));

C1s=SmoothBySlices(C1,[20 20],4);
C2s=SmoothBySlices(C2,[20 20],4);

P1su = cusackUnwrap(repmat(angle(C1s),[1 1 3]), repmat(abs(C1s),[1 1 3]));
P1su=P1su(:,:,2);

P2su = cusackUnwrap(repmat(angle(C2s),[1 1 3]), repmat(abs(C2s),[1 1 3]));
P2su=P2su(:,:,2);

fgmontagex([P1su P2su],[-pi pi]*4)
%%
E3c=SliCpNC(:,:,3);
E3cS=SmoothBySlices(E3c,[20 20],4);

sE3c=SliCpNC(:,:,2).*exp(1i*P2su/2);
sE3cS=SmoothBySlices(sE3c,[20 20],4);
fgmontagex([angle(E3cS) angle(sE3cS)])

dpE3S=sE3cS.*conj(E3cS);
dpE3S2=SmoothBySlices(dpE3S,[40 40],12);
fgmontage(angle(dpE3S2))

dP_EvenOdd=angle(dpE3S2);
%%
tmp=SliCpNC(800:1150,520:850,:);
[PDBaset, UpdatedB0Map_Hzt, UpdatedT2SMap_mst, s_valst, Fitted0t, PDBase0t]=...
    FitToModel_MPBD1CSf(tmp,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
tmp=SliCpNC;
tmp(:,:,[2 4])=tmp(:,:,[2 4]).*exp(-1i*dP_EvenOdd);
tmp=tmp(800:1150,520:850,:);

[PDBaset, UpdatedB0Map_Hzt, UpdatedT2SMap_mst, s_valst, Fitted0t, PDBase0t]=...
    FitToModel_MPBD1CSf(tmp,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
SliCpNC_EOp=SliCpNC;
SliCpNC_EOp(:,:,[2 4])=SliCpNC_EOp(:,:,[2 4]).*exp(-1i*dP_EvenOdd);

[PDBase_EOc, UpdatedB0Map_Hz_EOc, UpdatedT2SMap_ms_EOc, s_vals_EOc, Fitted0_EOc, PDBase0_EOc]=...
    FitToModel_MPBD1CSf(SliCpNC_EOp,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
B0_EOc_angle=UpdatedB0Map_Hz_EOc*pi/(1000/(2*ES_ms));
B0_EOc_W=s_vals_EOc(:,:,1);

B0_EOc_C=B0_EOc_W.*exp(1i*B0_EOc_angle);
B0_EOc_CS=SmoothBySlices(B0_EOc_C,[20 20],4);

B0_EOc_angle_u = cusackUnwrap(repmat(angle(B0_EOc_CS),[1 1 3]), repmat(abs(B0_EOc_CS),[1 1 3]));
B0_EOc_angle_u=B0_EOc_angle_u(:,:,2);

B0_EOc_u_Hz=B0_EOc_angle_u*(1000/(2*ES_ms))/pi;
%%
SmallVarT=50e6;
getSmallVars;
% save([BasePX 'SmallVars.mat'],SmallVars{:});
save([BasePX 'SmallVars40.mat'],SmallVars{:});
%%
nRO=size(B0_EOc_u_Hz,1);
%%
FPE=fft1cg(XSliCpN,2);
%%
InnerVoxDephasing=linspace(-pi,pi,size(FPE,2));
SPE=sinc(InnerVoxDephasing/(2*pi));
FPEC=FPE./SPE;
IPECorrected=ifft1cg(FPEC,2);
IPECorrectedCombined=squeeze(sum(IPECorrected.*conj(CurSens1),3));
% SliCpNC=squeeze(sum(XSliCpN.*conj(CurSens1),3));

%%
dB0dx=symD(B0_EOc_u_Hz,1);
dB0dy=symD(B0_EOc_u_Hz,2);
% FXSliCpN=fft1cg(XSliCpN,1);
FXSliCpN=fft1cg(IPECorrected,1);
%%
WhichPEIdxs=520:855;
nWhichPEIdxs=numel(WhichPEIdxs);

nTS=10;
% TimePointsCenters_ms=linspacecenters(0,ES_ms,nTS);
TimePointsCenters_ms=linspacecenters(-ES_ms/2,ES_ms/2,nTS);
TimePointsCenters_ms3=perm32(TimePointsCenters_ms);
TSC=exp(1i*2*pi*B0_EOc_u_Hz.*TimePointsCenters_ms3/1000);
%%
TSCSens_BeforeRODephasing=TSC.*perm43(CurSens1);
ROIdxs=floor(linspace(1,nRO+1,nTS+1));
SampMask=zeros(nRO,nTS);
for i=1:nTS
    SampMask(ROIdxs(i):(ROIdxs(i+1)-1),i)=1;
end
ROIdxsCenters=(ROIdxs(1:end-1)+ROIdxs(2:end))/2;
ROIdxsCentersN=ROIdxsCenters/nRO;
%%
PEPhaseDiffS=2*pi.*TimePointsCenters_ms3.*perm43(dB0dx)/1000;
PhaseDiffDueToImaging= (ROIdxsCentersN-0.5)*2*pi;
PhaseDiffDueToImaging3=perm32(PhaseDiffDueToImaging);

PEPhaseDiffTotal=PEPhaseDiffS+PhaseDiffDueToImaging3;
PEPhaseDiffTotalS=sinc(PEPhaseDiffTotal/(2*pi));

TSCSens=TSCSens_BeforeRODephasing.*PEPhaseDiffTotalS;

Ops={'fmac 0 0','fftc 1','fmac 1 0'}; % just multiply with TSCSens, FT along RO, multiply by sampling mask
ScriptFN=[BasePX 'Ops.txt'];
WriteLinopToFile(ScriptFN,Ops);
SampMaskP=repmat(perm42(SampMask),[1 nWhichPEIdxs 1 1]);
WhichEcho=1;
KspP=repmat(FXSliCpN(:,WhichPEIdxs,:,WhichEcho),[1 1 1 nTS]);
KspP=KspP.*SampMaskP;

KspPAllEchos=repmat(perm54(FXSliCpN),[1 1 1 nTS]);
KspPAllEchos=KspPAllEchos.*SampMaskP(:,1,:,:);

TSCSensP=perm43(TSCSens);

ImSz16=FillOnesTo16([nRO nWhichPEIdxs]);

ImSz16FN=[BasePX 'ImSz16'];
writecfl(ImSz16FN,ImSz16);

% KspPFN=[BasePX 'KspP'];
% writecfl(KspPFN,KspP);

TSCSensPFN=[BasePX 'TSCSensP'];
writecfl(TSCSensPFN,TSCSensP(:,WhichPEIdxs,:,:));

SampMaskPFN=[BasePX 'SampMaskP'];
writecfl(SampMaskPFN,SampMaskP);

for EchoI=1:4
    disp(EchoI);
    KspPFNE{EchoI}=[BasePX 'KspP_echo' num2str(EchoI)];
    writecfl(KspPFNE{EchoI},KspPAllEchos(:,WhichPEIdxs,:,:,EchoI));
end

mTSC=exp(-1i*2*pi*B0_EOc_u_Hz.*TimePointsCenters_ms3/1000);
mTSCSens_BeforeRODephasing=mTSC.*perm43(CurSens1);

mPEPhaseDiffTotal=-PEPhaseDiffS+PhaseDiffDueToImaging3;
mPEPhaseDiffTotalS=sinc(mPEPhaseDiffTotal/(2*pi));

mTSCSens=mTSCSens_BeforeRODephasing.*mPEPhaseDiffTotalS;
mTSCSensP=perm43(mTSCSens);

mTSCSensPFN=[BasePX 'mTSCSensP'];
writecfl(mTSCSensPFN,mTSCSensP(:,WhichPEIdxs,:,:));
%%
for EchoI=1:4
    disp(EchoI);
    Recx(:,:,EchoI)=bart(['picsS -m -R T:3:0:0.0001 ' ScriptFN],ImSz16FN,KspPFNE{EchoI},TSCSensPFN,SampMaskPFN);
    mRecx(:,:,EchoI)=bart(['picsS -m -R T:3:0:0.0001 ' ScriptFN],ImSz16FN,KspPFNE{EchoI},mTSCSensPFN,SampMaskPFN);
end
%%
RecxN=Recx.*grms(SliCpNC(:,WhichPEIdxs,:),1:2)./grms(Recx,1:2);
mRecxN=mRecx.*grms(SliCpNC(:,WhichPEIdxs,:),1:2)./grms(mRecx,1:2);
RecA=RecxN;
RecA(:,:,[2 4])=mRecxN(:,:,[2 4]);

RecB=mRecxN;
RecB(:,:,[2 4])=RecxN(:,:,[2 4]);

E3a=abs(SliCpNC(:,WhichPEIdxs,3));
sE3a=abs(sqrt(SliCpNC(:,WhichPEIdxs,2).*SliCpNC(:,WhichPEIdxs,4)));
dE3a=E3a-sE3a;
dE3ac=cat(3,E3a,sE3a,dE3a*5);

E3b=abs(Recx(:,:,3));
sE3b=abs(sqrt(Recx(:,:,2).*Recx(:,:,4)));
dE3b=E3b-sE3b;
dE3bc=cat(3,E3b,sE3b,dE3b*5);

E3c=abs(mRecx(:,:,3));
sE3c=abs(sqrt(mRecx(:,:,2).*mRecx(:,:,4)));
dE3c=E3c-sE3c;
dE3cc=cat(3,E3c,sE3c,dE3c*5);

E3d=abs(RecA(:,:,3));
sE3d=abs(sqrt(RecA(:,:,2).*RecA(:,:,4)));
dE3d=E3d-sE3d;
dE3dc=cat(3,E3d,sE3d,dE3d*5);

E3e=abs(RecB(:,:,3));
sE3e=abs(sqrt(RecB(:,:,2).*RecB(:,:,4)));
dE3e=E3e-sE3e;
dE3ec=cat(3,E3e,sE3e,dE3e*5);

E3All=cat(4,dE3ac,dE3bc,dE3cc,dE3dc,dE3ec);

E3All=E3All./grms(E3All(:,:,1,:),1:2);

WhichROIdxs=800:1150;

WhichROIdxs2=950:1050;
PEIdxs2=100:200;
fgmontagex(abs((E3All(WhichROIdxs2,PEIdxs2,:,:))))
%%
% tmp=RecA;
tmp=RecB;
tmp(:,:,[2 4])=tmp(:,:,[2 4]).*exp(-1i*dP_EvenOdd(:,WhichPEIdxs));
tmp=tmp(WhichROIdxs,:,:);
tmp=abs(tmp);

[PDBase_EOB0cm, UpdatedB0Map_Hz_EOB0cm, UpdatedT2SMap_ms_EOB0cm, s_vals_EOB0cm, Fitted0_EOB0cm, PDBase0_EOB0cm]=...
    FitToModel_MPBD1CSf(tmp,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
tmp=abs(SliCpNC(WhichROIdxs,WhichPEIdxs,:));
[PDBasetm, UpdatedB0Map_Hztm, UpdatedT2SMap_mstm, s_valstm, Fitted0tm, PDBase0tm]=...
    FitToModel_MPBD1CSf(tmp,WhichEchosToUseA,ES_ms,FirstTE_ms);
disp('ok');
%%
Before=SliCpNC(WhichROIdxs,WhichPEIdxs,:);
After=Fitted0_EOB0cm;
BeforeAfter=cat(4,Before,After);


%%
BeforeAfter=cat(4,Before,After);
%%
WhichEchoToShow=4;
fgmontagex(BeforeAfter(:,:,WhichEchoToShow,:))

fgmontagex(BeforeAfter(200:330,100:250,WhichEchoToShow,:))
%%
SmallVarT=50e6;
getSmallVars;
save([BasePX 'SmallVars2_40.mat'],SmallVars{:});
%%
BeforeAfterD=cat(4,Before,After_NoDephase,After);

fgmontagex(BeforeAfterD(:,:,WhichEchoToShow,:),'Size',[1 3])

fgmontagex(BeforeAfterD(200:330,100:250,WhichEchoToShow,:),'Size',[1 3])

%%
DPhase=abs(After_NoDephase-After);

fgmontagex(DPhase)

fgmontagex(DPhase(200:330,100:250,:))

%%
AA=bart(['linopScript ' ScriptFN],ImSz16,SliCpNC(:,:,WhichEcho),TSCSensP(:,WhichPEIdxs,:,:),SampMaskP);
AAX=bart(['linopScript -A ' ScriptFN],ImSz16,KspP,TSCSensP(:,WhichPEIdxs,:,:),SampMaskP);
%%
Rec1=bart(['picsS -d 5 -R W:3:0:0.0001 ' ScriptFN],ImSz16,KspP,TSCSensP(:,WhichPEIdxs,:,:),SampMaskP);

Rec1g=bart(['picsS -g -d 5 -R T:3:0:0.0001 ' ScriptFN],ImSz16,KspP,TSCSensP(:,WhichPEIdxs,:,:),SampMaskP);

Rec1=bart(['picsS -m -d 5 -R T:3:0:0.0001 ' ScriptFN],ImSz16FN,KspPFN,TSCSensPFN,SampMaskPFN);
Rec1g=bart(['picsS -m -g -d 5 -R T:3:0:0.0001 ' ScriptFN],ImSz16FN,KspPFN,TSCSensPFN,SampMaskPFN);