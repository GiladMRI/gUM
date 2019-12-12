ScanP='/autofs/cluster/kawin/Gilad/Ma_5x5/';
BaseFN='meas_MID134_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_5x5_FID8424_1mm_CD110';
RefFldMapP=[ScanP 'meas_MID135_BP_fieldmap_5echosX_FID8425_1mm' filesep];

BaseFN='meas_MID127_gBP_ep2d_bold_multiecho_ASL_SMS_Spi_5x5_FID8417_2mm_CD110';
RefFldMapP=[ScanP 'meas_MID128_BP_fieldmap_5echosX_FID8418_2mm' filesep];

MIDStr=BaseFN(6:11);
FN=[ScanP BaseFN '.dat'];
disp('ok');
%%
mainP=[ScanP BaseFN];
mkdir(mainP);

system(['chmod +777 -R ' mainP]);
disp([mainP ' Created']);

setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/bart-0.4.04b')

BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';
nccToUse=31;

Rows2Complex=@(X) X(1,:)+1i*X(2,:);

CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];
disp('ok 2');
%% Read raw
AData = mapVBVD(FN);
if(iscell(AData))
    ADatax=AData{end};
else
    ADatax=AData;
end
asSlice=ADatax.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[ADatax.hdr.Phoenix.sSliceArray.asSlice{:}];
end

nSlices=numel(ADatax.hdr.Phoenix.sSliceArray.asSlice);
try
    WipMemBlock=ADatax.hdr.MeasYaps.sWiPMemBlock;
catch
    WipMemBlock=ADatax.hdr.MeasYaps.sWipMemBlock;
end
MB=1; %WipMemBlock.alFree{9};
% MB
nSlicesNoMB=nSlices/MB;

% nRepsHdr=1+ADatax.hdr.Meas.lRepetitions;
nRepsHdr=1+ADatax.hdr.MeasYaps.lRepetitions;
disp('Read raw');
%%
TE0_ms=2.38;

TrajType=WipMemBlock.adFree{13};
% ResType=floor((TrajType-10)/2)+1; % 1.9,1.3
ResType=1;
% TimingType=mod(TrajType-10,4)+1; % 6ms, 8ms
% load('GAll68x.mat');

LL=getLines('SpiGradVec.cpp');
LLL=LL{4}(2:end-3);
LLL2=strrep(LLL,'f','');
ZZ=eval(['[' LLL2 ']']);
ZZ=ZZ*ZZ(end);
ZZ(end)=0;
ZZZ=ZZ(1:2500)+1i*ZZ(2501:end);
GAll=ZZZ.';

GTrajaCBase=GAll;
% GNav=load('GNav1ms.mat');
% GNav=GNav.GNav;
% GTrajaCBase=GAll(:,TimingType,ResType);

% if(TimingType==1 || TimingType==3)
%     nInnerShots=8;
% else
nInnerShots=5;
% end
%%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;
%% Do something with the noise!
% noise=AData{1}.noise();
%%
for s=1:nSlices
    try
        SlbLoc(1,s)=asSlice(s).sPosition.dSag;
    catch
        SlbLoc(1,s)=0;
    end
    try
        SlbLoc(2,s)=asSlice(s).sPosition.dCor;
    catch
        SlbLoc(2,s)=0;
    end
    try
        SlbLoc(3,s)=asSlice(s).sPosition.dTra;
    catch
        SlbLoc(3,s)=0;
    end
end

RotMat = transpose(Quat2RotMat(ADatax.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;

[U IA IB]=unique(ADatax.image.slicePos(3,:));
Qoffset=ADatax.image.slicePos(1:3,IA);

Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

try
    FOVx=ADatax.hdr.Meas.ReadFOV;
catch
    FOVx=ADatax.hdr.Config.ReadFoV;
end

FOVx=220;

dFOV=FOVx/1000;

FOV_mm=FOVx;

disp('ok');

GradReduceFac=WipMemBlock.adFree{4};
PhiRotPerRep=WipMemBlock.adFree{8};
AcqDwellTime_us=WipMemBlock.adFree{13}/1000;

dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;
%%
RLocs=load([RefFldMapP 'Locs.mat']);
RefLocs=RLocs.RotatedLocs;
CurzLocs=RotatedLocs(3,:);
RefzLocs=RefLocs(3,:);

RefMaps=load([RefFldMapP 'B0T2S.mat']);
RefSens=load([RefFldMapP 'Sens.mat']);
disp('Loaded ref');
%%
GTrajaC=GTrajaCBase/GradReduceFac;
GNavC=GNav/GradReduceFac;

for j=1:5
    GTrajaCM(:,j)=[zeros((j-1)*50,1); GTrajaC; zeros((5-j)*50,1); GNavC; -GNavC];
end

for j=1:5
    g=GTrajaCM(:,j);
    k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
    s=diff(g)/GradDwellTime_ms;
    
    kK=k*FOV_mm/1000/2/pi;
    kkM(:,j)=kK;
end

Kmax=ceil(max(abs(kK)));
NTrg=Kmax*2;
res_mm=FOVx/(max(abs(kK))*2);
%
figure;
plot(kK);
axis square;
axis equal;
ylabel('k');
%%
TrgSz=[NTrg NTrg];
OnesSens=ones(TrgSz);
%%
Extra=4;
kKA=[zeros(4,5); kkM; repmat(kkM(end,:),[10 1])]; % in us

% Delay_us=-2.5; % positive data rightward %-6; % clockwise
% Delays_us=-2:0.1:4;

Delays_us=1.6;
dlyI=1;
% for dlyI=1:numel(Delays_us)
    disp('-------');
    disp(dlyI);
    Delay_us=Delays_us(dlyI);
    AcqTimePoints_us=Extra*GradDwellTime_us+Delay_us+(0:AcqDwellTime_us:52000);

    nRepsToUse=nRepsHdr-1;
    
    clear TrajM
    for i=1:nRepsToUse
        CurJitterIdx=mod(i-1,5)+1;
        Traj=interp1((0:(size(kKA,1)-1))*GradDwellTime_us,kKA(:,CurJitterIdx),AcqTimePoints_us);
        TrajM(i,:)=Traj.*exp(-1i*2*pi*PhiRotPerRep/360*(i-1));
    end
    
    AcqLen_us=AcqTimePoints_us(end)-AcqTimePoints_us(1)+AcqDwellTime_us;
    
    BARTTraj=cat(3,real(TrajM),imag(TrajM),imag(TrajM)*0);
    BARTTrajP=permute(BARTTraj,[3 2 1]);

    kx=BARTTrajP(1,:,:)*2*pi;
    ky=BARTTrajP(2,:,:)*2*pi;

    modx=exp(-1i.*(dx*kx+dy*ky));
    disp('ok mod');
    
    nAcqPoints=numel(Traj);

Sz=TrgSz;
Sz16=FillOnesTo16(Sz);
%%
SlicesToRead=[5];
for SlicesToRead=1:nSlices
RepsToRead=1:(nRepsHdr-1);
% ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,ADCsToRead,:,:,:,:,:,:);
ADataIsL=ADatax.image(:,:,:,:,SlicesToRead,3,:,:,RepsToRead,:,:,:,:,:,:,:,:);

ADataIsL=permute(ADataIsL,[1 2 9 11 5 3:4 6:8 10]);
ADataIsL=CombineDims(ADataIsL,[4 1]);

disp('Read data');
%%
ChRMS=grmss(ADataIsL,[1 3]);
[SChRMS ChOrd]=sort(ChRMS,'descend');
% Data=squeeze(ADataIsL(:,ChOrd(1),1:nRepsToUse));

Ch2D=CombineDims(ADataIsL,[3 1]);
[~,S,sccmtx] = svd(Ch2D(1:10:end,:),'econ');
clear Ch2D 

ncc=31;
ADataIsLCC=single(zeros([size(ADataIsL,1) ncc size(ADataIsL,3)]));
for i=1:ncc
    ADataIsLCC(:,i,:)=sum(ADataIsL.*permute(sccmtx(:,i),[3 1 4 5 6 7 8 9 2]),2);
end
DataC=permute(ADataIsLCC,[1 3 2]);
disp('ok cc');

    DataC=permute(ADataIsLCC,[1 3 2]);
    %% Recon each channel separately, no CC
nTrajToUse=size(BARTTrajP,2);
TrajPartToUse=0+(1:2000);
RepsToUse=1:(nRepsHdr-1);

DataPC=permute(ADataIsL(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 3 5 6 7 8 2]).*modx;
OnesSensC=repmat(OnesSens,[1 1 1 1 1 1 1 size(DataPC,8)]);

Rec1p=bart('pics -t ',BARTTrajP(:,TrajPartToUse,RepsToUse),(DataPC(:,TrajPartToUse,RepsToUse,1,1,1,1,:)),OnesSensC);
% fgmontage(grmss(Rec1p,8));removeTicks;
% title([num2str(numel(RepsToUse)) ' shots data x ' num2str(numel(TrajPartToUse)) 'points, pre channel recon, RMS']);
%%
Rec1pRMS=grmss(Rec1p,8);
[Out B BN]=CalcSlicesSNR(Rec1pRMS,false,7);
Msk=imfillholesBySlices(~BN);
Msk=getLargestComponent(Msk);
% Msk=imfillholesBySlices(Rec1pRMS>0.001);
se = strel('disk', 3);
DMsk=imdilate(Msk,se,'same');
SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(Rec1p).*DMsk,0,TrgSz);
SelfSens1=SelfSens(:,:,:,1);
disp('ok SelfSens');
%% Recon each channel separately, CC
DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;

disp('ok rec pre channel cc');
%
SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
SensCC=permute43(SensCC);
disp('ok SensCC');

DataPC=permute(DataC(1:nTrajToUse,:,:,:,:,:,:,:,:),[4 1 2 5 6 7 8 3]).*modx;
DataCCP=permute(DataPC,[1:3 8 4:7]);
disp('ok 0');
%%
TightMask=DMsk;
%% Try to get all TSC, multi-shot, for B0,T2*
CurReps=1:39;

nTS_THLR=15;

TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

nTSMed=nTS_THLR;

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC=FillOnesTo16([TrgSz 1 1 1 1 nTS_THLR]);
%
TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

dTS_planned_ms=2.5;

nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);
TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

TimePointsMed_ms=linspace(0,AcqTimePoints_us(nPointsMed)/1000,nTSMed);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);
%
nTraj=numel(TrajPartMed);
TotalAcqTime_ms=AcqDwellTime_us*nTraj/1000;

nPointsNoNav=floor(50000/AcqDwellTime_us);
NoNavTime_ms=nPointsNoNav*AcqDwellTime_us/1000;
NoNavB=zeros(1,nTraj);
NoNavB(1:nPointsNoNav)=1;

% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
[TSBMed, dT_Med, TimePointsR_Med]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTSMed);
dT_Med_ms=dT_Med*NoNavTime_ms;
FirstT_Med_ms=TimePointsR_Med(1)*NoNavTime_ms;
TimePoints_Med_ms=TimePointsR_Med*NoNavTime_ms;
TimePoints_Med_ms3=permute(TimePoints_Med_ms,[1 3 2]);
TSBMed(nPointsMed,1)=0;
TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);
KernsPMMed=getKernsFromTrajM(TrajM(1:(nRepsHdr-1),TrajPartMed),Sz,TSBMed);

ScriptFN_AllTS=[BaseSP 'nuftAllTSC_N.txt'];
Sz16AllTSC=FillOnesTo16([Sz 1 1 1 1 nTSMed]);

disp('Prepared TSB,Kerns for Med');
%%
% TSB_THLR=GetTSCoeffsByLinear(nPointsMed,nTS_THLR);
[TSB_THLR, dT_THLR, TimePointsR_THLR]=GetTSCoeffsByLinearWithPlateau(nPointsNoNav,nTS_THLR);
dT_THLR_ms=dT_THLR*NoNavTime_ms;
FirstT_THLR_ms=TimePointsR_THLR(1)*NoNavTime_ms;
TSB_THLR(nPointsMed,1)=0;
TSB_THLRP=permute(TSB_THLR,[3 1 4 5 6 7 2]);
Sz16_THLR=FillOnesTo16([Sz 1 1 1 1 nTS_THLR]);
KernsP_TSTHLR=getKernsFromTrajM(TrajM(1:(nRepsHdr-1),TrajPartMed),Sz,TSB_THLR);
disp('Prepared TSB,Kerns for THLR');
%%
clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
    disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
end
%%
TightMask=DMsk;
THLR_lambda=10;
RhoStr=[' -u ' num2str(1e-3) ' '];
%%
RepsSets={1:(nRepsHdr-1), 1:3};
nRepsSets=numel(RepsSets);
for rs=1:nRepsSets
    CurReps=RepsSets{rs};
    disp(['RepSet ' num2str(rs) ' ! ' datestr(now)]);
    tmp=bart(['picsS -w 1 -m ' RhoStr ' -R K:64:3:' num2str(THLR_lambda) ':2:1:0:6 ' ScriptFN_AllTS],Sz16_THLR,NoNavB.*DataCCP(:,TrajPartMed,CurReps,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurReps),TSB_THLRP,1,...
        sum(KernsP_TSTHLR(:,:,CurReps,:,:,:,:),3));
    THLRMultiShot_RS(:,:,:,rs)=squeeze(tmp);
end
%%
WhichTSToUs=2:12;
clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
for rs=1:nRepsSets
    [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=...
        FitToModel_MPBD1CSf(THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
end
% ShowAbsAngle(THLRMultiShot_RS(:,:,2:4:end,:));ylabel('39 shots             3 shots','FontSize',16);title('THLR');
% ShowAbsAngle(Fitted0_RS(:,:,2:4:end,:));      ylabel('39 shots             3 shots','FontSize',16);title('1CS fitted THLR');
% ShowAbsAngle(THLRMultiShot_RS(:,:,2:4:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('THLR');
% ShowAbsAngle(Fitted0_RS(:,:,2:4:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('1CS fitted THLR');
% ShowAbsAngle(THLRMultiShot_RS(:,:,2:4:end,1));ylabel('3 shots','FontSize',16);title('THLR');
% ShowAbsAngle(Fitted0_RS(:,:,2:4:end,1));ylabel('3 shots','FontSize',16);title('1CS fitted THLR');
%%
T2sCenters=1:100;
T2Sbins=[-Inf T2sCenters-0.5 Inf];

nComponentsToUse=4;

clear TSCxPMed_RS
for rs=1:nRepsSets
    TSCMed_RS(:,:,:,rs)=exp(-TimePoints_Med_ms3./UpdatedT2SMap_ms_RS(:,:,rs)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,rs).*TimePoints_Med_ms3/1e3);
    TSCxPMed_RS(:,:,rs,1,1,1,:)=permute(TSCMed_RS(:,:,:,rs),[1 2 7 6 5 4 3]);
    
    TempRMS=grmss(THLRMultiShot_RS(:,:,:,rs),3:30);
    CurT2SMap=UpdatedT2SMap_ms_RS(:,:,rs);
    [HT2s,edges,T2sbin] = histcounts(CurT2SMap(:),T2Sbins);
    WHT2s=HT2s*0;
    for i=1:numel(HT2s)
        WHT2s(i)=sum(TempRMS(T2sbin==i));
    end
    SWHT2s=max(WHT2s,sum(WHT2s)*0.03/numel(WHT2s));

    clear WDecays WDecays4
    for i=2:numel(HT2s)-1
        WDecays(i,:)=SWHT2s(i)*exp(-TimePointsMed_ms./T2sCenters(i-1));
    end
    [WUd,WSd,WVd]=svd(WDecays,'econ');
    WVd_RS(:,:,rs)=WVd;
    
    CompsP_RS(rs,1,1,1,1,:,:)=permute(WVd_RS(:,1:nComponentsToUse,rs),[7:-1:3 2 1]);
end
disp('Prepared for LLR RepSets');

LLR_lambda=0.1;
RhoStr=[' -u ' num2str(1e-3) ' '];
% RhoStr=[' -u ' num2str(1e-1) ' '];
BlkSz=4;
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];
Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);

TSCxPMedOnlyB0_RS=exp(1i.*angle(TSCxPMed_RS));
%%
% [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=FitToModel_MPBD1CSf(THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,TotalAcqTime_ms,TE0_ms);
% save(['For_NU_MPBD2_S' num2str(SliI) '_.mat'],'THLRMultiShot_RS','PDBase_RS','UpdatedB0Map_RS','UpdatedT2SMap_ms_RS','Fitted0_RS');
%%
clear Rec_CompgB0_C_RS
for rs=2:nRepsSets
    CurRep=RepsSets{rs};
    Rec_CompgB0_C_RS{rs}=bart(['picsS -m -w 1 ' RhoStr ' -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
        Sz16CompgB0,NoNavB.*DataCCP(:,TrajPartMed,CurRep,1:nccToUse),...
        SensCC(:,:,:,1:nccToUse),STraj3MMed(:,:,CurRep),TSBPMed,TSCxPMedOnlyB0_RS(:,:,rs,1,1,1,:),...
        sum(KernsPMMed(:,:,CurRep,:,:,:,:),3),CompsP_RS(rs,1,1,1,1,:,:));
end
%%
Rec_CompgB0_RS_M=cat(4,Rec_CompgB0_C_RS{:});
Rec_CompgB0_RS_MX=perm43(squeeze(sum(Rec_CompgB0_RS_M.*gpermute(CompsP_RS,[4 1]),6)));
Rec_CompgB0_RS_MX=Rec_CompgB0_RS_MX.*perm43(squeeze(TSCxPMedOnlyB0_RS));

WhichTSToUse_LLR=2:20;

for rs=1:nRepsSets
    [PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs), Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
        FitToModel_MPBD1CSf(Rec_CompgB0_RS_MX(:,:,:,rs),WhichTSToUse_LLR,dT_Med_ms,TE0_ms+FirstT_Med_ms(1));
end
% ShowAbsAngle(Rec_CompgB0_RS_MX(:,:,2:5:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('LLR');
% ShowAbsAngle(Fitted0_RS_LLR(:,:,2:5:end,:));ylabel('20 shots      10 shots    5 shots     3 shots','FontSize',16);title('1CS Fitted LLR');
% ShowAbsAngle(Rec_CompgB0_RS_MX(:,:,2:5:end,2));ylabel('3 shots','FontSize',16);title('LLR');
% ShowAbsAngle(Fitted0_RS_LLR(:,:,2:5:end,2));ylabel('3 shots','FontSize',16);title('1CS Fitted LLR');
%%
CurReps=RepsSets{2};
CurSig=DataCCP(:,TrajPartMed,CurReps,1:nccToUse);
CurKernsA=KernsPMMed(:,:,CurReps,:,:,:,:);
CurSens=SensCC(:,:,:,1:nccToUse);
CurTraj=STraj3MMed(:,:,CurReps);
SliI=SlicesToRead;
save([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat'],'CurSig','CurKernsA','CurSens','CurTraj','Fitted0_RS_LLR',...
    'UpdatedB0Map_RS_LLR','UpdatedT2SMap_ms_RS_LLR','TotalAcqTime_ms','TE0_ms','THLRMultiShot_RS','Fitted0_RS','Rec_CompgB0_RS_MX');
%%
% DataForSlice{SlicesToRead}={UpdatedB0Map,UpdatedT2SMap_ms,THLRMultiShot,SelfSens1,sccmtx};

close all;

end % Slices loop
%%
% save([mainP filesep 'DataForSlice.mat'],'DataForSlice');
% % load([mainP filesep 'DataForSlice.mat']);
% %%
% clear KernsByRepMed              	1513758480    		1513.76MB
% clear KernsPMMed                 	1513754112    		1513.75MB
% clear ADataIsL                   	499200000     		499.20MB
% clear DataForSlice               	492372000     		492.37MB
% clear ADataIsLCC                 	483600000     		483.60MB
% clear DataC                      	483600000     		483.60MB
% clear DataCCP                    	439640760     		439.64MB
% clear DataPC                     	439640760     		439.64MB
% clear InnerShotData              	193440000     		193.44MB
% clear InnerShotDataIn            	193440000     		193.44MB
% save([mainP filesep 'CurStatus.mat']);
%%
SliI=7;
for SliI=1:nSlices
    disp(SliI);
    QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
    WhichTSToUse_LLR=2:20;
    Fitted0_RS_LLRS(:,:,:,SliI)=QQ.Fitted0_RS_LLR(:,:,:,2);
    Fitted0_RSS(:,:,:,:,SliI)=QQ.Fitted0_RS;
    UpdatedB0Map_RS_LLRS(:,:,:,SliI)=QQ.UpdatedB0Map_RS_LLR;
    UpdatedT2SMap_ms_RS_LLRS(:,:,:,SliI)=QQ.UpdatedT2SMap_ms_RS_LLR;
    rs=2;
    [PDBase_RS_LLR(:,:,rs), UpdatedB0Map_RS_LLR(:,:,rs), UpdatedT2SMap_ms_RS_LLR(:,:,rs), s_vals_RS_LLR(:,:,:,rs), Fitted0_RS_LLR(:,:,:,rs), PDBase0_RS_LLR(:,:,rs)]=...
            FitToModel_MPBD1CSf(QQ.Rec_CompgB0_RS_MX(:,:,:,rs),WhichTSToUse_LLR,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    PDBase0_RS_LLRS(:,:,SliI)=PDBase0_RS_LLR(:,:,rs);
    
    WhichTSToUs=2:12;
    clear PDBase_RS UpdatedB0Map_RS UpdatedT2SMap_ms_RS s_vals_RS Fitted0_RS PDBase0_RS
    for rs=1:nRepsSets
        [PDBase_RS(:,:,rs), UpdatedB0Map_RS(:,:,rs), UpdatedT2SMap_ms_RS(:,:,rs), s_vals_RS(:,:,:,rs), Fitted0_RS(:,:,:,rs), PDBase0_RS(:,:,rs)]=...
            FitToModel_MPBD1CSf(QQ.THLRMultiShot_RS(:,:,:,rs),WhichTSToUs,dT_THLR_ms,TE0_ms+FirstT_THLR_ms);
    end
    UpdatedB0Map_RSS(:,:,:,SliI)=UpdatedB0Map_RS;
    UpdatedT2SMap_ms_RSS(:,:,:,SliI)=UpdatedT2SMap_ms_RS;
end
%%
fgmontagex(Fitted0_RS_LLRS(:,:,12,ROrd))
fgmontagex(Fitted0_RS_LLRS(:,:,12,:))
%%
figure;
subplot(2,2,1);gmontage(abs(PDBase0_RS_LLRS(:,:,SliI)),[0 5.5e-4]);title('3-shot LLR');
subplot(2,2,2);gmontage(angle(PDBase0_RS_LLRS(:,:,SliI)),[-pi pi]);
subplot(2,2,3);gmontage(UpdatedB0Map_RS_LLRS(:,:,2,SliI),[-300 300]);
subplot(2,2,4);gmontage(UpdatedT2SMap_ms_RS_LLRS(:,:,2,SliI),[0 100]);
%%
for SliI=1:nSlices
    CurSlizLoc=CurzLocs(Ord(SliI));
    [MinRefD, CorrespondingRefSliIdx]=min(abs(RefzLocs-CurSlizLoc));
    RefPDBase0=rot90(padarray(RefMaps.PDBase0(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
    RefT2SM=rot90(padarray(RefMaps.UpdatedT2SMap_ms(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
    RefPDBase0=min(abs(RefPDBase0),grmss(RefPDBase0)*6).*exp(1i*angle(RefPDBase0));
    RefB0M=rot90(padarray(RefMaps.B0M_Hz(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));

    MidEcho_ms=20;
    RefMidEcho=RefPDBase0.*exp(-MidEcho_ms./RefT2SM);

    RefMidEchorS(:,:,SliI)=imresizeBySlices(RefMidEcho,Sz);
    RefB0MrS(:,:,SliI)=imresizeBySlices(RefB0M,Sz);
end

fgmontagex(UpdatedB0Map_RSS(:,:,1,ROrd),[-200 200]);title('From THLR multishot')
fgmontagex(RefB0MrS(:,:,ROrd),[-200 200]);title('From GRE-ME')
%%
BB=imread('batch016757_out_0.0174.png');
BBM=BB(1:116,116*4*5+(1:116*4),1);
BBP=BB(1:116,116*4*7+(1:116*4),1);
BBC=double(BBM).*exp(1i*2*pi*double(BBP)/256);
BBCX=PartitionDim(BBC,2,4);
%%
EPTIRecP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Recon/Human/';

S16=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_16_GE.mat']);
S6=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_6_GE.mat']);
EPTIRef6=imresize(gflip(S6.im_recon(:,:,1:21:end),2),Sz);
EPTIRef16=imresize(gflip(S16.im_recon(:,:,1:21:end),2),Sz);

fgmontagex(EPTIRef16Fr(:,:,1:21:end),[0 7]);title('EPTI fitted')
fgmontagex(EPTIRef6,[0 7]);
fgmontagex(EPTIRef6Fr(:,:,1:21:end),[0 7]);title('EPTI fitted')

%%
figure;
subplot(2,2,1);gmontage(abs(PDBase0_EPTI),[0 7]);title('Fitted EPTI');
subplot(2,2,2);gmontage(angle(PDBase0_EPTI),[-pi pi]);
subplot(2,2,3);gmontage(UpdatedB0Map_EPTI,[-300 300]);
subplot(2,2,4);gmontage(UpdatedT2SMap_ms_EPTI,[0 100]);
%%
WhichTSToUse_EPTI=10:60;
dT_EPTI_ms=0.78;

YY=gflip(S16.im_recon(:,:,1:1:end),2);
[PDBase_EPTI, UpdatedB0Map_EPTI, UpdatedT2SMap_ms_EPTI, s_vals_EPTI, Fitted0_EPTI, PDBase0_EPTI]=...
    FitToModel_MPBD1CSf(YY,WhichTSToUse_EPTI,dT_EPTI_ms,TE0_ms+dT_EPTI_ms);
EPTIRef16F=Fitted0_EPTI;
EPTIRef16Fr=imresize(EPTIRef16F,Sz);

YY=gflip(S6.im_recon(:,:,1:1:end),2);
[PDBase_EPTI, UpdatedB0Map_EPTI, UpdatedT2SMap_ms_EPTI, s_vals_EPTI, Fitted0_EPTI, PDBase0_EPTI]=...
    FitToModel_MPBD1CSf(YY,WhichTSToUse_EPTI,dT_EPTI_ms,TE0_ms+dT_EPTI_ms);
EPTIRef6F=Fitted0_EPTI;
EPTIRef6Fr=imresize(EPTIRef6F,Sz);
%%
SliI=4;
XX=cat(4,Fitted0_RSS(:,:,1:4:end,1,SliI)/5.5e-4,Fitted0_RS_LLRS(:,:,1:6:end,SliI)/5.5e-4,EPTIRef16/6.5);
fgmontagex(perm43(XX),[0 1])

fgmontagex(XX(:,:,3,:),[0 1],'Size',[2 2])
%%
SliI=10;
XX=cat(4,Fitted0_RSS(:,:,1:4:end,1,SliI)/5.5e-4,Fitted0_RS_LLRS(:,:,1:6:end,SliI)/5.5e-4,EPTIRef6/6.5);
fgmontagex(perm43(XX),[0 1])

fgmontagex(XX(:,:,3,:),[0 1],'Size',[1 3])
%%
fgmontagex(BBCX);title('SKEPTIC MLN, 1-shot');
fgmontagex(QQ.Fitted0_RS(:,:,1:4:end,1),[0 5.5e-4]);title('SKEPTIC THLR Fitted, 39-shot');
fgmontagex(QQ.Fitted0_RS_LLR(:,:,1:6:end,2),[0 5.5e-4]);title('SKEPTIC LLR Fitted, 3-shot');
% fgmontagex(Fitted0_RS_LLRS(:,:,12,3))
%%
% XX=cat(4,QQ.Fitted0_RS(:,:,1:4:end,1)/5.5e-4,QQ.Fitted0_RS_LLR(:,:,1:6:end,2)/5.5e-4,BBCX/255,EPTIRef/6.5);
XX=cat(4,QQ.Fitted0_RS(:,:,1:4:end,1)/5.5e-4,QQ.Fitted0_RS_LLR(:,:,1:6:end,2)/5.5e-4,Fitted0_MLN/255,EPTIRef/6.5);
fgmontagex(perm43(XX),[0 1]);
tmp=get(gca,'Children');
WSz=size(tmp.CData);
WSzX=WSz./[4 4];

StartX=10;
StartY=10;
text(StartX,StartY,'SKEPTIC THLR Fitted, 39-shot','FontSize',16,'Color','red')
text(StartX,StartY+WSzX(2),'SKEPTIC LLR Fitted, 3-shot','FontSize',16,'Color','red')
text(StartX,StartY+2*WSzX(2),'SKEPTIC MLN, 1-shot','FontSize',16,'Color','red')
text(StartX,StartY+3*WSzX(2),'EPTI subspace 3-shot','FontSize',16,'Color','red')
%%
fgmontagex(XX(:,:,2,:),[0 1])
StartX=10;
StartY=5;
text(StartX,StartY,'SKEPTIC THLR Fitted, 39-shot','FontSize',16,'Color','red')
text(StartX,StartY+WSzX(2),'SKEPTIC MLN, 1-shot','FontSize',16,'Color','red')
text(StartX+WSzX(1),StartY,'SKEPTIC LLR Fitted, 3-shot','FontSize',16,'Color','red')
text(StartX+WSzX(1),StartY+WSzX(2),'EPTI subspace 3-shot','FontSize',16,'Color','red')
%%
[PDBase_MLN, UpdatedB0Map_MLN, UpdatedT2SMap_ms_MLN, s_vals_MLN, Fitted0_MLN, PDBase0_MLN]=...
        FitToModel_MPBD1CSf(BBCX,1:4,15,3);
%%
% QQ=load([RefFldMapP 'sTwixX.mat']);
%%
CurSlizLoc=CurzLocs(Ord(SliI));
[MinRefD, CorrespondingRefSliIdx]=min(abs(RefzLocs-CurSlizLoc));
% CorrespondingRefSliIdx=28;
RefPDBase0=rot90(padarray(RefMaps.PDBase0(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
RefT2SM=rot90(padarray(RefMaps.UpdatedT2SMap_ms(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));
RefPDBase0=min(abs(RefPDBase0),grmss(RefPDBase0)*6).*exp(1i*angle(RefPDBase0));
RefB0M=rot90(padarray(RefMaps.B0M_Hz(:,:,CorrespondingRefSliIdx),[0 2 0],'both'));

MidEcho_ms=20;
RefMidEcho=RefPDBase0.*exp(-MidEcho_ms./RefT2SM);

RefMidEchor=imresizeBySlices(RefMidEcho,Sz);
RefB0Mr=imresizeBySlices(RefB0M,Sz);
%%
fgmontage(RefMidEchor);removeTicks
fgmontage(Fitted0_RS(:,:,7,1));removeTicks

fgmontage(RefB0Mr,[-300 300]);removeTicks
fgmontage(UpdatedB0Map_RS(:,:,1),[-300 300]);removeTicks
%%
nNeighbors=12;
lambda=0.1;
% lambda=1;
CurRep=1;
SliI=4;
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);
CurSig=QQ.CurSig(1,:,CurRep,:);
CurSens=QQ.CurSens;
CurTraj=TrajM(CurRep,:);
TrajPartToUse=1:45501;
CurTraj=(CurTraj(:,TrajPartToUse(1:3:end))+CurTraj(:,TrajPartToUse(2:3:end))+CurTraj(:,TrajPartToUse(3:3:end)))/3;
CurSig=squeeze(QQ.CurSig(1,:,CurRep,:));


CurSig=(CurSig(TrajPartToUse(1:3:end),:)+CurSig(TrajPartToUse(2:3:end),:)+CurSig(TrajPartToUse(3:3:end),:))/3;
CurB0M=UpdatedB0Map_RSS(:,:,1,SliI);

TrajTimePoints_ms=(0:size(CurTraj,2)-1)*AcqDwellTime_us*3/1000;
TrajTimePoints_ms3=perm32(TrajTimePoints_ms);

TrgnTS=7;
TrgTSTimePoints_ms=linspace(0,TrajTimePoints_ms(end),TrgnTS);
TrgTSTimePoints_ms3=perm32(TrgTSTimePoints_ms);

nCh=size(CurSens,4);
xx=3;
yy=4;
W=grmss(Fitted0_RSS(:,:,:,1,SliI),3);
W=W./grmss(W);
disp('ok');
%%
NeiM=zeros([Sz nNeighbors]);
xS=-Sz(1)/2:Sz(1)/2-1;
yS=-Sz(2)/2:Sz(2)/2-1;
VecX=(-Sz(1)/2:Sz(1)/2-1)/Sz(1);
VecY=(-Sz(1)/2:Sz(1)/2-1)/Sz(2);
% iAHAtikM=zeros(numel(xS),numel(yS),nCh*nNeighbors,nCh*nNeighbors);
% CoeffsM=zeros(numel(xS),numel(yS),nCh*nNeighbors);
CoeffsM=zeros(numel(xS),numel(yS),nCh*nNeighbors,TrgnTS);
ValsM=zeros(numel(xS),numel(yS),nNeighbors,nCh);
xSIToDo=1:numel(xS);
ySIToDo=1:numel(yS);
xSIToDo=floor(Sz(1)/2)+(-7:7);
ySIToDo=floor(Sz(2)/2)+(-7:7);
for xI=xSIToDo
    disp([num2str(xI,'%d ') ' ' datestr(now)]);
    for yI=ySIToDo
        xx=xS(xI);
        yy=yS(yI);
        cc=xx+1i*yy;
        
%         X1=exp(1i*linspace(-pi,pi,Sz(1))*xx);
%         Y1=exp(1i*linspace(-pi,pi,Sz(2))*yy).';
        X1=exp(1i*2*pi*VecX*xx);
        Y1=exp(1i*2*pi*VecY*yy).';
        TrgM=X1.*Y1;
        TrgMWB0=W.*TrgM.*exp(1i*2*pi*CurB0M.*TrgTSTimePoints_ms3/1000);
        TrgMWB02=reshape(TrgMWB0,Sz(1)*Sz(2),[]);

        dTraj=CurTraj-cc;
        [SdTraj, TOrd]=sort(abs(dTraj),'ascend');
        
        CurNeiLocs=CurTraj(TOrd(1:nNeighbors));
        for n=1:nNeighbors
            NeiM(:,:,n)=exp(1i*linspace(-pi,pi,Sz(1))*real(CurNeiLocs(n))).*(exp(1i*linspace(-pi,pi,Sz(2))*imag(CurNeiLocs(n))).');
        end
        CurB0MEffect=exp(1i*2*pi*CurB0M.*TrajTimePoints_ms3(1,1,TOrd(1:nNeighbors))/1000);
        NeiMWB0=W.*NeiM.*CurB0MEffect;
        NeiMWB0WSens=NeiMWB0.*CurSens;
        NeiMWB0WSens2=reshape(NeiMWB0WSens,Sz(1)*Sz(2),[]);
        
        AHy=NeiMWB0WSens2'*TrgMWB02;
        AHA=NeiMWB0WSens2'*NeiMWB0WSens2;
        AHAtik=AHA+eye(size(AHA,1))*lambda;
        
        Coeffs=AHAtik\AHy;
%         grmss(AHy)
%         grmss(TrgMWB02)
%         grmss(AHy-AHAtik*Coeffs)
% grmss(AHy-AHA*Coeffs)
% grmss(TrgMWB02-NeiMWB0WSens2*Coeffs)
%         Coeffs=NeiMWB0WSens2\TrgWB02;
% RR=NeiMWB0WSens2*Coeffs;
% ShowAbsAngle(reshape(RR,[Sz 7]))
% ShowAbsAngle(reshape(TrgMWB02,[Sz 7]))
        CoeffsM(xI,yI,:,:)=Coeffs;
        
        Vals=CurSig(TOrd(1:nNeighbors),:);
        ValsM(xI,yI,:,:)=Vals;
%         iAHAtik=inv(AHAtik);
%         Coeffs=iAHAtik*Aty;
%         CoeffsM(xI,yI,:)=Coeffs;
%         iAHAtikM(xI,yI,:,:)=iAHAtik;
    end
end
disp('Finished');
%%
OutF=zeros([Sz TrgnTS]);
for xI=xSIToDo
%     disp([num2str(xI,'%d ') ' ' datestr(now)]);
    for yI=ySIToDo
        Vals=squeeze(ValsM(xI,yI,:,:));
        Vals1=Vals(:);
        OutF(xI,yI,:)=squeeze(CoeffsM(xI,yI,:,:)).'*Vals1;
    end
end
disp('ok');
%
OutI=ifft2cg(OutF);
ShowAbsAngle(OutI)
%%
IValsM=ifft2cg(ValsM);
%%
% save('ForMLN.mat','TrajM','Sensr','sig','UpdatedB0Map_RS','UpdatedT2SMap_ms_RS');
BaseBaseOutP='/autofs/cluster/kawin/Gilad/TF/MEx/';
mkdir(BaseBaseOutP);
system(['chmod -R 777 ' BaseBaseOutP]);
%%
SliI=1;
for SliI=1:nSlices
    disp(SliI);
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);

Sensr=QQ.CurSens;
sig=QQ.CurSig;
%
BaseOutDir=[BaseBaseOutP 'Sli' num2str(SliI) filesep];
mkdir(BaseOutDir);
system(['chmod -R 777 ' BaseOutDir]);
nCh=13;
batchSize=16;
nTS=15;
SensCC=squeeze(Sensr(:,:,1:nCh)); % [X Y Ch]
SensMsk=grmss(SensCC,3)>0.01; % [X Y]
save([BaseOutDir 'SensCC1.mat'],'SensCC','SensMsk');

WhichRep=1;
% AcqDwellTime_us=2*1.1;
AcqDwellTime_us=3*1.1;
% TrajPartToUse=1:24000;
TrajPartToUse=1:45501;
% CurSig=     sig(1,TrajPartToUse(1:2:end),WhichRep,1:nCh)+...
%             sig(1,TrajPartToUse(2:2:end),WhichRep,1:nCh);
CurSig=     sig(1,TrajPartToUse(1:3:end),WhichRep,1:nCh)+...
            sig(1,TrajPartToUse(2:3:end),WhichRep,1:nCh)+...
            sig(1,TrajPartToUse(3:3:end),WhichRep,1:nCh);
CurSig=squeeze(CurSig);
tmp=Row(CurSig);
tmp2=[real(tmp) imag(tmp)]*600;
Data=tmp2;
Data(batchSize,end)=0;
save([BaseOutDir 'RealDataForNN.mat'],'Data');
nTrajA=size(CurSig,1);
TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TS_TimePoints3=permute(TS_TimePoints,[1 3 2]);
TSBF=GetTSCoeffsByLinear(nTrajA,nTS).';
WhichRSToUse=1;
% TSC=exp(-TS_TimePoints3./UpdatedT2SMap_ms_RS(:,:,WhichRSToUse)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,WhichRSToUse).*TS_TimePoints3/1e3);
% TSBF: [15×7162 double]
% TSC: [128×128×15 double]
% B0_Hz=UpdatedB0Map_RS(:,:,WhichRSToUse);
B0_Hz=RefB0MrS(:,:,SliI);
save([BaseOutDir 'B0TS.mat'],'TSBF','B0_Hz');
%
% Traj=(TrajM(WhichRep,TrajPartToUse(1:2:end))+TrajM(WhichRep,TrajPartToUse(2:2:end)))/2;
Traj=(TrajM(WhichRep,TrajPartToUse(1:3:end))+TrajM(WhichRep,TrajPartToUse(2:3:end))+TrajM(WhichRep,TrajPartToUse(3:3:end)))/3;
Sz128=gsize(SensCC,1:2);
clear Trajm2
Trajm2(1,:)=real(Traj);
Trajm2(2,:)=imag(Traj);
[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
save([BaseOutDir 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');

TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TSstr=strrep(num2str(TS_TimePoints,'%3.5f,'),' ','');
TSstr=TSstr(1:end-1);
% TSstr=['TimePoints_ms ' TSstr(1:end-1)];

St=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForS10.txt');
St.TimePoints_ms=TSstr;
TS_TimePointsForRec=linspace(0,TimePoints_ms(end),8);
TS_TimePointsForRec_Str=strrep(num2str(TS_TimePointsForRec,'%3.5f,'),' ','');
TS_TimePointsForRec_Str=TS_TimePointsForRec_Str(1:end-1);
St.TimePointsForRec_ms=TS_TimePointsForRec_Str;
St.SessionNameBase=['ME_S' num2str(SliI)];
St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
St.BaseTSDataP=BaseOutDir;
St.BaseNUFTDataP=BaseOutDir;
% disp(TSstr)
ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForMES' num2str(SliI) '.txt'];
Txt=gStruct2txt(St,ParamsOutFn);
disp('ok')
end
%%
BaseBaseOutP='/autofs/cluster/kawin/Gilad/TF/ME3/';
mkdir(BaseBaseOutP);
system(['chmod -R 777 ' BaseBaseOutP]);
disp(['Created ' BaseBaseOutP]);
%%
SliI=10;
% for SliI=1:nSlices
    disp(SliI);
QQ=load([mainP filesep 'For_NU_MPBD3_S' num2str(SliI) '.mat']);

Sensr=QQ.CurSens;
sig=QQ.CurSig;
%
BaseOutDir=[BaseBaseOutP 'Sli' num2str(SliI) filesep];
mkdir(BaseOutDir);
system(['chmod -R 777 ' BaseOutDir]);
nCh=13;
batchSize=16;
nTS=15;
SensCC=squeeze(Sensr(:,:,1:nCh)); % [X Y Ch]
SensMsk=grmss(SensCC,3)>0.01; % [X Y]
save([BaseOutDir 'SensCC1.mat'],'SensCC','SensMsk');

WhichReps=1:3;
% AcqDwellTime_us=2*1.1;
% AcqDwellTime_us=3*1.1;
AcqDwellTime_usz=3*1.1/3;
% TrajPartToUse=1:24000;
TrajPartToUse=1:45501;
% CurSig=     sig(1,TrajPartToUse(1:2:end),WhichRep,1:nCh)+...
%             sig(1,TrajPartToUse(2:2:end),WhichRep,1:nCh);
CurSig=     sig(1,TrajPartToUse(1:3:end),WhichReps,1:nCh)+...
            sig(1,TrajPartToUse(2:3:end),WhichReps,1:nCh)+...
            sig(1,TrajPartToUse(3:3:end),WhichReps,1:nCh);
CurSig=squeeze(CurSig);
CurSig=permute(CurSig,[2 1 3]);
CurSig=reshape(CurSig,prod(gsize(CurSig,1:2)),size(CurSig,3));
tmp=Row(CurSig);
tmp2=[real(tmp) imag(tmp)]*600;
Data=tmp2;
Data(batchSize,end)=0;
save([BaseOutDir 'RealDataForNN.mat'],'Data');
nTrajA=size(CurSig,1);
TimePoints_ms=(1:nTrajA)*AcqDwellTime_usz/1000;
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TS_TimePoints3=permute(TS_TimePoints,[1 3 2]);
TSBF=GetTSCoeffsByLinear(nTrajA,nTS).';
WhichRSToUse=1;
% TSC=exp(-TS_TimePoints3./UpdatedT2SMap_ms_RS(:,:,WhichRSToUse)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,WhichRSToUse).*TS_TimePoints3/1e3);
% TSBF: [15×7162 double]
% TSC: [128×128×15 double]
% B0_Hz=UpdatedB0Map_RS(:,:,WhichRSToUse);
B0_Hz=RefB0MrS(:,:,SliI);
save([BaseOutDir 'B0TS.mat'],'TSBF','B0_Hz');
%
% Traj=(TrajM(WhichRep,TrajPartToUse(1:2:end))+TrajM(WhichRep,TrajPartToUse(2:2:end)))/2;
TrajZ=(TrajM(WhichReps,TrajPartToUse(1:3:end))+TrajM(WhichReps,TrajPartToUse(2:3:end))+TrajM(WhichReps,TrajPartToUse(3:3:end)))/3;
TrajZ=TrajZ(:).';
Sz128=gsize(SensCC,1:2);
clear Trajm2
Trajm2(1,:)=real(TrajZ);
Trajm2(2,:)=imag(TrajZ);
[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
save([BaseOutDir 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');

TimePoints_ms=(1:nTrajA)*AcqDwellTime_usz/1000;
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TSstr=strrep(num2str(TS_TimePoints,'%3.5f,'),' ','');
TSstr=TSstr(1:end-1);
% TSstr=['TimePoints_ms ' TSstr(1:end-1)];

St=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForS10.txt');
St.TimePoints_ms=TSstr;
TS_TimePointsForRec=linspace(0,TimePoints_ms(end),8);
TS_TimePointsForRec_Str=strrep(num2str(TS_TimePointsForRec,'%3.5f,'),' ','');
TS_TimePointsForRec_Str=TS_TimePointsForRec_Str(1:end-1);
St.TimePointsForRec_ms=TS_TimePointsForRec_Str;
St.SessionNameBase=['ME3_S' num2str(SliI)];
St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
St.BaseTSDataP=BaseOutDir;
St.BaseNUFTDataP=BaseOutDir;
% disp(TSstr)
ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForME3S' num2str(SliI) '.txt'];
Txt=gStruct2txt(St,ParamsOutFn);
disp('ok')
% end
%%
StBaser=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForMES10r.txt');

srezNP='/autofs/cluster/kawin/Gilad/TF/srezN/';
DD=dir(srezNP);
DD=DD([DD.isdir]);
DD=DD(3:end);
DD=DD(strhas({DD.name},'checkpoint'));
%%
SliI=9;
for SliI=9:14
SliStr=num2str(SliI);
CurSP=DD(strhas({DD.name},['ME_S' SliStr '__'])).name;

CurSt=StBaser;
CurSt.LoadAndRunOnData_checkpointP=[srezNP CurSP];
CurSt.LoadAndRunOnData_Prefix= ['/autofs/cluster/kawin/Gilad/TF/ME/s' SliStr 'r'];
CurSt.LoadAndRunOnData_OutP= ['/autofs/cluster/kawin/Gilad/TF/ME/S' SliStr]; 
CurSt.BaseTSDataP=['/autofs/cluster/kawin/Gilad/TF/MEx/Sli' SliStr '/'];
CurSt.BaseNUFTDataP=['/autofs/cluster/kawin/Gilad/TF/MEx/Sli' SliStr '/'];

mkdir(CurSt.LoadAndRunOnData_OutP);
system(['chmod +777 -R ' CurSt.LoadAndRunOnData_OutP]);

BaseBaseOutP='/autofs/cluster/kawin/Gilad/TF/MEx/';
BaseOutDir=[BaseBaseOutP 'Sli' num2str(SliI) filesep];
RealDataFNa=[BaseOutDir 'RealDataForNN.mat'];
RealDataFNb=[CurSt.LoadAndRunOnData_Prefix '01.mat'];
copyfile(RealDataFNa,RealDataFNb);
                   
ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForME1S' SliStr 'r.txt'];
Txt=gStruct2txt(CurSt,ParamsOutFn);
end
disp('Wrote run params.txt');
%%
for SliI=9:14
    SliStr=num2str(SliI);
    RDataFN=['/autofs/cluster/kawin/Gilad/TF/ME/S' SliStr '/OnRealData01.mat'];
    RR=load(RDataFN);
    RRC=RR.x(:,:,:,1)+1i*RR.x(:,:,:,2);
    RRC1=squeeze(RRC(1,:,:));
    RRC1M=PartitionDim(RRC1,2,8);
    RRC1MS(:,:,:,SliI-8)=RRC1M;
end
%%
% Rec1pRMS=grmss(Rec1p,8);
% [Out B BN]=CalcSlicesSNR(Rec1pRMS,false,7);
% Msk=imfillholesBySlices(~BN);
% Msk=getLargestComponent(Msk);
% % Msk=imfillholesBySlices(Rec1pRMS>0.001);
% se = strel('disk', 3);
% DMsk=imdilate(Msk,se,'same');

for SliI=1:6
    WhichTSToUse_MLN=2:7;
dT_MLN_ms=7.15016;
[PDBase_MLNS(:,:,SliI), UpdatedB0Map_MLNS(:,:,SliI), UpdatedT2SMap_ms_MLNS(:,:,SliI), s_vals_MLNS(:,:,:,SliI), Fitted0_MLNS(:,:,:,SliI), PDBase0_MLNS(:,:,SliI)]=...
    FitToModel_MPBD1CSf(RRC1MS(:,:,:,SliI),WhichTSToUse_MLN,dT_MLN_ms,TE0_ms+dT_MLN_ms);
end
%%
fgmontagex(perm43(Fitted0_MLNS),[0 1.4])
%%
fgmontagex(UpdatedT2SMap_ms_MLNS.*(squeeze(s_vals_MLNS(:,:,1,:))>0.7),[0 200]);colormap hot
%%
fgmontagex(angle(PDBase0_MLNS.*(squeeze(s_vals_MLNS(:,:,1,:))>0.7)));

%%
OutGifFN=[mainP filesep 'MLN' '.gif'];
delete(OutGifFN);
Mx=1.3;
figure;pause(1)
fgmontagex(Fitted0_MLNS(:,:,1,:),[0 Mx]);
gmontage(Fitted0_MLNS(:,:,1,:),[0 Mx]);
removeTicks;daspect([1 1 1]);set(get(gcf,'Children'),'Position',[0.01 0.01 0.98 0.9]);pause(4)
gif(OutGifFN)
pause(1)
for i=2:size(Fitted0_MLNS,3)
    gmontage(Fitted0_MLNS(:,:,i,:),[0 Mx]);
    removeTicks;daspect([1 1 1]);set(get(gcf,'Children'),'Position',[0.01 0.01 0.98 0.9]);
    gif
end
disp(OutGifFN);

%%
% See ParamsForMES10r.txt
% HowManyToRun 2
RR=load('/autofs/cluster/kawin/Gilad/TF/ME/OnRealData01.mat');
RRC=RR.x(:,:,:,1)+1i*RR.x(:,:,:,2);
RRC1=squeeze(RRC(1,:,:));
RRC1M=PartitionDim(RRC1,2,8);
%%
WhichTSToUse_MLN=2:7;
dT_MLN_ms=7.15016;
[PDBase_MLN, UpdatedB0Map_MLN, UpdatedT2SMap_ms_MLN, s_vals_MLN, Fitted0_MLN, PDBase0_MLN]=...
    FitToModel_MPBD1CSf(RRC1M,WhichTSToUse_MLN,dT_MLN_ms,TE0_ms+dT_MLN_ms);
%%
figure;
subplot(2,2,1);gmontage(abs(PDBase0_MLN),[0 1]);title('Fitted MLN');
subplot(2,2,2);gmontage(angle(PDBase0_MLN),[-pi pi]);
subplot(2,2,3);gmontage(UpdatedB0Map_MLN,[-300 300]);
subplot(2,2,4);gmontage(UpdatedT2SMap_ms_MLN,[0 100]);

