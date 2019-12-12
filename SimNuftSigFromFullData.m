directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
FNBase='meas_MID00868_FID32103_ep2d_ge_sms1_EPTI_1p9_fully';
% FNBase='meas_MID00892_FID32127_ep2d_ge_sms1_EPTI_1p9_fully_Cor';

filename = [directory_rawdata,FNBase,'.dat'];
save_filename = 'Calib';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);
%%
slice=9;
load([OutP 'kdata_Slice' num2str(slice) '.mat'],'kdata');
%%
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;
CTo2Rows=@(X) [real(X);imag(X)];
CTo3Rows=@(X) [real(X);imag(X);imag(X)*0];

FOV=228; % meas.prot.dReadoutFOV
FOV_mm=FOV;
Sz=[120 120];
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));

pf_echo=0;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo);

rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/'))
%%
kdatax=kdata(1:120,1:120,:,:);
idatax=ifft2cg(kdatax);
%% Sensitivity map
% Sens=RunESPIRiTForSensMapsMultiMap(squeeze(idatax(:,:,1,:)),0,Sz);
Sens=RunESPIRiTForSensMapsMultiMap(squeeze(idatax(:,:,1,:)),14,Sz);
Sens=perm43(Sens(:,:,:,1));

CombinedAllEchos=sum(idatax.*conj(Sens),4);
%%
% ES_ms=meas.prot.iEffectiveEpiEchoSpacing/1000;
ES_ms=0.7273;
%% Fit
TE0_ms=9;
WhichTSToUs=5:50;
[~, UpdatedB0MapA, UpdatedT2SMap_msA, s_valsA, FittedA, PDBase0A]=...
    FitToModel_MPBD1CSf(CombinedAllEchos,WhichTSToUs,ES_ms,TE0_ms+ES_ms);
%%
G2mm=load('Grads2mmb.mat');
% QQ=load('GAll1p9mmVD1PAT3Pause.mat');
QQ=load('GAll1p9mmVD1PAT3Pauseb.mat');
AllGrads=cat(1,G2mm.Grads{:}).';
AllGrads=cat(2,AllGrads,QQ.GAll(:,7:8));
G11mm=load('Grads11mmb.mat');
AllGrads11=cat(1,G11mm.Grads{:}).';
AllGrads=cat(2,AllGrads,AllGrads11);

%%
GTraj16=tmp.';
%% EPTI R=4 jitter 2
KOrd=[-60:4:60 -58:4:60 -60:4:60];
NPerRO=floor(48000/numel(KOrd));
ROPoints=linspace(-60,60,NPerRO);
ROP=repmat(ROPoints,[numel(KOrd) 1]).';
PEP=repmat(KOrd.',[1 NPerRO]).';

% Traj=(ROP(:)+1i.*PEP(:)).';
Traj=(1i*ROP(:)+PEP(:)).';
STraj=Traj;

nTraj=numel(Traj);
nTrajBase=nTraj;
STrajd=STraj;
disp('EPTI');
%% Simulate sig from complete kData
load('GAll1p9mmVD1PAT3.mat');
% GAll 4800           6
GTraj=GAll(:,6)/1.02;
%%
TrajIdx=1;
for TrajIdx=1:12
    disp('TrajIdx');
GTraj=AllGrads(:,TrajIdx);

% GTraj=GAll(:,1);
% 
% GTraj=GAll(:,2);

% GTraj=GTraj16;

nEchosData=size(idatax,3);
TotalAcqTime_ms=ES_ms*nEchosData;
%
k=cumsum([0; GTraj])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kK=k*FOV_mm/1000/2/pi;

Kmax=ceil(max(abs(kK)));
% NTrg=Kmax*2;
res_mm=FOV_mm/(max(abs(kK))*2);
AcqDwellTime_us=1.1;
AcqTimePoints_us=(0:AcqDwellTime_us:TotalAcqTime_ms*1000);

Traj=interp1((0:(size(kK,1)-1))*GradDwellTime_us,kK,AcqTimePoints_us);
Traj=Traj(~isnan(Traj));
STraj=Traj;
nTraj=numel(Traj);
nTrajBase=nTraj;
AcqTimePoints_usx=AcqTimePoints_us;

STrajd=STraj;

Rs=interp1(1:nTraj,abs(Traj),linspace(1,nTraj,nEchosData));

disp('Prepared traj');
% %% Delay
% Trajd=interp1((0:(size(kK,1)-1))*GradDwellTime_us,kK,AcqTimePoints_us+2);
% Trajd=Trajd(~isnan(Trajd));
% STrajd=Trajd;
% STrajd(end+1:nTraj)=0;
% disp('Delay');
% %% distort
% kKd=kK;
% kKd(1600+(1:1600))=kKd(1600+(1:1600))+0.4;
% kKd(2*1600+(1:1600))=kKd(2*1600+(1:1600))+0.8;
% Trajd=interp1((0:(size(kK,1)-1))*GradDwellTime_us,kKd,AcqTimePoints_us);
% Trajd=Trajd(~isnan(Traj));
% STrajd=Trajd;
% disp('Distorted');
% %%
% I2Spiral=kK*0;
% I2Spiral(1600+(1:1600))=1;
% I3Spiral=kK*0;
% I3Spiral(1*1600+(1:1600))=1;
% 
% I2I=interp1((0:(size(kK,1)-1))*GradDwellTime_us,I2Spiral,AcqTimePoints_us);
% I2I=I2I(~isnan(Traj));
% I2IB=I2I>0;
% 
% I3I=interp1((0:(size(kK,1)-1))*GradDwellTime_us,I3Spiral,AcqTimePoints_us);
% I3I=I3I(~isnan(Traj));
% I3IB=I3I>0;
% %% Multi shot
% kKf=flip(kK,1);
% Trajf=-interp1((0:(size(kKf,1)-1))*GradDwellTime_us,kKf,AcqTimePoints_us);
% Trajf=Trajf(~isnan(Traj));
% STrajf=Trajf;
% 
% STrajf=Traj.*exp(1i*2*pi/3);
% STrajf2=Traj.*exp(1i*2*pi*2/3);
% 
% STrajf=Traj.*exp(1i*2*pi*(140/360));
% STrajf2=Traj.*exp(1i*2*pi*(280/360));
% 
% STrajX=[STraj; STrajf; STrajf2];
% 
% % STrajX=[STraj; STrajf];
% STraj=STrajX(:).';
% nTraj=numel(STraj);
% Traj=STraj;
% 
% STrajd=STraj;
% AcqTimePoints_usx=AcqTimePoints_us*nTrajBase/nTraj;
% disp('Applied multishot');
%%
disp('Caluclating sig...');
STraj3=CTo3Rows(STraj);
STraj3d=CTo3Rows(STrajd);
Sigx=bart('nufft ',STraj3d,perm73(idatax));
disp('Caluclated sig');
% Siga=Sigx;
% Sigx(1,I3IB,:,:,:,:,:,:)=Sigx(1,I3IB,:,:,:,:,:,:).*exp(1i*pi/5);
%%
nTSx=21;
TSBa=GetTSCoeffsByLinear(nTraj,nEchosData);
TSBb=GetTSCoeffsByLinear(nEchosData,nTSx);
Sigy=sum(squeeze(perm73(Sigx)).*TSBa,2); % Traj 1 Channels
% Sigy=sum(Sigx.*TSBa,2); % Traj 1 Channels

TSB=GetTSCoeffsByLinear(nTraj,nTSx);
KernsPMMed=getKernsFromTrajM(STraj,Sz,TSB);

TimePointsMed_ms=linspace(0,AcqTimePoints_us(end)/1000,nTSx);
TimePointsMed_ms3=permute(TimePointsMed_ms,[1 3 2]);

TimePointsAllEchos_ms=linspace(0,AcqTimePoints_us(end)/1000,nEchosData);
TimePointsAllEchos_ms3=permute(TimePointsAllEchos_ms,[1 3 2]);

BaseSP='/autofs/space/daisy_002/users/Gilad/gUM/';
disp('Prepared kerns');
% Now recon subspace according to perfect B0?
T2svalues_ms=linspace(5,200,200);
Decays=exp(-TimePointsMed_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');
%
nComponentsToUse=2;
ScriptFN_CompgBo=[BaseSP 'nuftCompgB0_N.txt'];

TSCxPMedOnlyB0=exp(-1i*2*pi.*UpdatedB0MapA.*perm72(TimePointsMed_ms)/1000);
% TSCxPMedOnlyB0=exp(-1i*2*pi.*UpdatedB0MapA*0.5.*perm72(TimePointsMed_ms)/1000);
% TSCxPMedOnlyB0=exp(-1i*2*pi.*(UpdatedB0MapA+10).*perm72(TimePointsMed_ms)/1000);

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);

TSBPMed=permute(TSB,[3 1 4 5 6 7 2]);

LLR_lambda=0.1;
RhoStr=[' -u ' num2str(1e-3) ' '];
BlkSz=4;
disp('Prepared LLR');
%%
CurData=permute(Sigy,[4 1 2 3]); % 1       nTraj           1          Channels
% RegStr='-R T:3:0:0.000001';
RegStr='-R W:3:0:0.001';
Rec_CompgB0_C=bart(['picsS -m ' RhoStr ' -d 5 -S -b ' num2str(BlkSz) ' ' RegStr ' ' ScriptFN_CompgBo],...
    Sz16CompgB0,CurData,Sens,STraj3,TSBPMed,TSCxPMedOnlyB0,KernsPMMed,CompsP);

Rec_CompgB0_MX=squeeze(sum(Rec_CompgB0_C.*CompsP,6).*TSCxPMedOnlyB0);
% Fit
WhichTSToUsA=7:18;
dT_msA=nTraj*AcqDwellTime_us/(1000*(nTSx-1));
[~, UpdatedB0MapA, UpdatedT2SMap_msA, s_valsA, FittedA, PDBase0A]=...
    FitToModel_MPBD1CSf(Rec_CompgB0_MX,WhichTSToUsA,dT_msA,TE0_ms);

UpdatedT2SMap_msAx=min(300,max(5,abs(UpdatedT2SMap_msA)));
m0=abs(PDBase0A);
m0=min(m0,median(m0)*56);
c0=m0.*exp(1i*angle(PDBase0A));
FitAFull=c0.*exp(-1i*2*pi*UpdatedB0MapA.*(TE0_ms+TimePointsAllEchos_ms3)/1000).*exp(-(TE0_ms+TimePointsAllEchos_ms3)./UpdatedT2SMap_msAx);
% FitAFull=FitAFull.*grmss(Rec_CompgB0_MX)/grmss(FitAFull);
%%
TSBbP=permute(TSBb,[3 4 2 1]);
Rec_CompgB0_MX_AllEchos=squeeze(sum(Rec_CompgB0_MX.*TSBbP,3));

FitAFull=squeeze(sum(FittedA.*TSBbP,3));

% Both=cat(4,CombinedAllEchos,Rec_CompgB0_MX_AllEchos);
Both=cat(4,CombinedAllEchos,Rec_CompgB0_MX_AllEchos,FitAFull);
%%
kCombinedAllEchos=fft2cg(CombinedAllEchos);
kRec_CompgB0_MX_AllEchos=fft2cg(Rec_CompgB0_MX_AllEchos);
kFitRec_CompgB0_MX_AllEchos=fft2cg(FitAFull);

dktRec=kRec_CompgB0_MX_AllEchos-kCombinedAllEchos;
adktRec=abs(dktRec);

dktRecA=kFitRec_CompgB0_MX_AllEchos-kCombinedAllEchos;
adktRecA=abs(dktRecA);

NadktRec=adktRec./abs(kCombinedAllEchos);
[X,Y]=ndgrid(linspace(-1,1,Sz(1)),linspace(-1,1,Sz(1)));
R=sqrt(X.^2+Y.^2);
Edges=linspace(0,1,61);
[H,B,B2]=histcounts(R(:),Edges);
Rx=reshape(B2,Sz);
for r=1:60
    CurB=   Rx==r;
    QQ=Reshape4d22d(perm43(adktRec),CurB);
    mkt(r,:)=mean(QQ,1);
    
    QQ=Reshape4d22d(perm43(adktRecA),CurB);
    mktA(r,:)=mean(QQ,1);
    
    QQ=Reshape4d22d(perm43(abs(kCombinedAllEchos)),CurB);
    Nfac(r,:)=mean(QQ,1);
    
    QQ=Reshape4d22d(perm43(NadktRec),CurB);
    Nmkt(r,:)=mean(QQ,1);
end
N2mkt=mkt./Nfac;
N2mktA=mktA./Nfac;
disp('got N2mkt');
%%
% fgmontagex(flip(N2mkt,1))
fgmontagex(flip(N2mkt,1),[0 1]);colormap jet;
xlabel('echos ->');
ylabel('k radius');
title('Normalized error');
set(get(gcf,'Children'),'Position',[0.1 0.1 0.8 0.8]);colorbar
hold on;
% plot(1:nEchosData,Sz(1)/2-Rs);
% plot(1:nEchosData,Sz(1)/2-Rs);
% plot(linspace(1,nEchosData,numel(kK)),Sz(1)/2-abs(kK));
N2mktC{TrajIdx}=N2mkt;
N2mktCX{TrajIdx}={mkt mktA Nfac};
end
%%
save('N2mktC.mat','N2mktC','N2mktCX');
%%
N2mkt_3Spi=N2mkt;
N2mkt_Spi6ms=N2mkt;
N2mkt_3SpiAndMinusFlip=N2mkt;

N2mkt_2Shot_Spi6msAnd120=N2mkt;
N2mkt_2Shot_Spi6msAnd120240=N2mkt;

N2mkt_Spi16ms=N2mkt;

N2mkt_3Shot_Spi16ms1480280=N2mkt;

save('TrajErrorAnalysis.mat','N2mkt_3Spi','N2mkt_Spi6ms','N2mkt_3SpiAndMinusFlip','N2mkt_2Shot_Spi6msAnd120',...
    'N2mkt_2Shot_Spi6msAnd120240','N2mkt_Spi16ms','N2mkt_3Shot_Spi16ms1480280','GTraj16');
%%
N2mkt36R=N2mkt3./N2mkt6;
fgmontage(flip(N2mkt36R,1),[0 2]);colorbar;colormap jet

title('3spirals / 8spirals')
%%
nTS_THLR=15;

nTSMed=nTS_THLR;

Sz16AllTSC=FillOnesTo16([TrgSz 1 1 1 1 nTS_THLR]);
%
TrajPartMed=1:nTrajToUse;

nPointsMed=numel(TrajPartMed);

dTS_planned_ms=2.5;

nTSMed=ceil((nPointsMed+1)*AcqDwellTime_us/1000/dTS_planned_ms);
% TSBMed=GetTSCoeffsByLinear(nPointsMed,nTSMed);
% TSBPMed=permute(TSBMed,[3 1 4 5 6 7 2]);

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
disp('Prepared TSB,Kerns for THLR');
%%
clear STraj3MMed
for CurRep=1:(nRepsHdr-1)
%     disp(CurRep);
    STraj=TrajM(CurRep,TrajPartMed);
    STraj3MMed(:,:,CurRep)=CTo3Rows(STraj);
    STraj3M(:,:,CurRep)=CTo3Rows(TrajM(CurRep,:));
end
%%
%%
RepsForKerns=[1 13 25];
% RepsForKerns=1:3;
% RepsForKerns=1:(nRepsHdr-1);
if(~exist('KernsPMMed','var'))
    KernsPMMed=getKernsFromTrajM(TrajM(RepsForKerns,TrajPartMed),Sz,TSBMed);
    KernsP_TSTHLR=getKernsFromTrajM(TrajM(RepsForKerns,TrajPartMed),Sz,TSB_THLR);
end

%%
QQ=imread('/autofs/cluster/kawin/Gilad/GIRF.png');
QQ=double(QQ);
NQQ=256-QQ;
AQQ=NQQ(50:270,55:380,:);
AQQ(1:50,200:end,:)=1;
AQQ12=AQQ(:,:,1).*AQQ(:,:,2);
AQQyy=AQQ12>20000;
for i=1:size(AQQyy,2)
    F=find(AQQyy(:,i));
    if(~isempty(F))
        GIRFyy(i)=mean(size(AQQyy,1)-F);
    else
        GIRFyy(i)=NaN;
    end
end

f=fit(find(~isnan(GIRFyy)).',GIRFyy(~isnan(GIRFyy)).','smoothingspline');
sGIRFyy=f(1:size(AQQyy,2));
NGIRFyy=(sGIRFyy-27)/7;
NGIRFyy=NGIRFyy./sum(NGIRFyy);
save('NGIRFyy.mat','NGIRFyy');
GXs=(1:326)-126;

NGIRFyySym=[zeros(75,1); NGIRFyy];
%%
GTraj_us=interp1(1:numel(GTraj),GTraj,1:0.1:numel(GTraj)).';

cGTraj_usr=conv(real(GTraj_us),NGIRFyySym,'same');
cGTraj_usi=conv(imag(GTraj_us),NGIRFyySym,'same');
cGTraj_us=cGTraj_usr+1i*cGTraj_usi;

k1=cumsum([0; GTraj_us])*0.001*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kK1=k1*FOV_mm/1000/2/pi;

kc=cumsum([0; cGTraj_us])*0.001*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kKc=kc*FOV_mm/1000/2/pi;

figure; plot(real(kK1(1000:19000)));hold on
plot(real(kKc(1000:19000)),'r');

k=cumsum([0; GTraj])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m  
kK=k*FOV_mm/1000/2/pi;

GTraj=GAll(:,1);


%%
save('GTraj16.mat','GTraj16');
TrajToWrite=GTraj16;
Ttl=PadStringWithBlanks('Traj16',50);
FN='/autofs/space/daisy_002/users/Gilad/gUM/Traj16.grd';

fid=fopen(FN,'wb');
fwrite(fid,numel(TrajToWrite),'int32');
fwrite(fid,uint8(Ttl),'uint8');
fwrite(fid,real(TrajToWrite),'float32');
fwrite(fid,imag(TrajToWrite),'float32');
fclose(fid);
%%