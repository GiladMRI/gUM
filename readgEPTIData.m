clear;
close all;
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%%
Pattern=[[31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90],
[32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91],
[33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92],
[34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93],
[35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94],
[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99,101,103,105,107,109,111,113,115,117,119],
[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98,100,102,104,106,108,110,112,114,116,118,120],
[119,117,115,113,111,109,107,105,103,101, 99, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1],
[120,118,116,114,112,110,108,106,104,102,100, 98, 96, 94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 74, 72, 70, 68, 66, 64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2],
[5, 8, 12, 29, 43, 45, 57, 59, 73, 87, 88, 92, 94,101,104,114,112,110,109, 89, 85, 71, 63, 54, 48, 41, 32, 22, 20, 6, 2, 4, 11, 24, 37, 52, 69, 76, 77, 80, 83, 97,116,117,120,106, 99, 81, 70, 67, 65, 61, 55, 51, 34, 33, 26, 19, 16, 15],
[13, 15, 29, 35, 39, 42, 43, 45, 50, 59, 68, 93,102,105,117,119,111,108, 99, 96, 87, 83, 81, 80, 78, 64, 63, 61, 55, 22, 5, 8, 18, 24, 25, 37, 47, 53, 62, 70, 89, 90,106,113,115,103,101, 92, 75, 72, 66, 58, 57, 56, 33, 28, 26, 20, 11, 6],
[24, 35, 37, 52, 54, 57, 69, 75, 85, 95,104,106,110,113,116,112, 93, 88, 83, 66, 64, 63, 61, 58, 33, 25, 22, 12, 11, 9, 1, 3, 4, 5, 15, 31, 67, 70, 77, 89, 96, 99,101,108,115,102, 91, 80, 73, 62, 59, 49, 47, 45, 40, 29, 27, 20, 18, 16],
[5, 8, 12, 29, 43, 45, 57, 59, 73, 87, 88, 92, 94,101,104,114,112,110,109, 89, 85, 71, 63, 54, 48, 41, 32, 22, 20, 6, 2, 4, 11, 24, 37, 52, 69, 76, 77, 80, 83, 97,116,117,120,106, 99, 81, 70, 67, 65, 61, 55, 51, 34, 33, 26, 19, 16, 15]];
%%
% directory = '/autofs/cluster/kawin/FuyixueWang/EPTI/EPTI_Rec_pack/';
% directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
% directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
% directory_rawdata = '/autofs/cluster/kawin/Gilad/gep_Phantom/';
% directory_rawdata = '/autofs/space/daisy_001/users/data/Gilad/gep_CL/';
directory_rawdata = '/autofs/space/daisy_001/users/data/Gilad/gEPTI_Try1/';
% directory_rawdata = '/autofs/space/daisy_001/users/data/Gilad/SCEPTI_AM_14Nov19/';
% filename_calib = 'Calib';
%% Get Parameters
FNBase='meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns';
FNBaseCalib='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
% filename_calib = 'Human_Calib_1p3_5shot';
% filename_calib = 'Human_Calib_1p9';

FNBase='meas_MID00872_FID32107_ep2d_ge_EPTI_1p9_1shot_4dyns';
FNBaseCalib='meas_MID00870_FID32105_ep2d_ge_EPTI_1p9_1shot_calib';

FNBase='meas_MID00932_FID42872_gepBase_T0_10rep';
% FNBase='meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns';
% % FNBase='meas_MID00878_FID32113_ep2d_ge_EPTI_1p9_5shot_4dyns';
% FNBaseCalib='meas_MID00874_FID32109_ep2d_ge_EPTI_1p9_3and5shot_calib';
% % 
% FNBase='meas_MID00903_FID32138_ep2d_ge_EPTI_1p9_3shot_4dyns_Cor';
% % FNBase='meas_MID00908_FID32143_ep2d_ge_EPTI_1p9_5shot_4dyns_Cor';
% FNBaseCalib='meas_MID00905_FID32140_ep2d_ge_EPTI_1p9_3and5shot_calib_Cor';

FNBase='meas_MID01946_FID43887_ep2d_ge_EPTI_1p9_3shot_16dyns';
FNBaseCalib='meas_MID01944_FID43885_ep2d_ge_EPTI_1p9_3and5shot_calib';

FNBase='meas_MID00540_FID48437_gEPTI_1Slice_Try1';
FNBase='meas_MID00131_FID48696_gEPTI_1Slice_Try2';
FNBase='meas_MID00470_FID49035_gEPTI_1Slice_Try2';

% FNBase='meas_MID00688_FID49241_gEPTI_1Slice_Try2';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);

% CalibP=[directory_rawdata FNBaseCalib filesep];
filename = [directory_rawdata,FNBase,'.dat'];
%%
apodization_para=0;

% read parameters
pf_echo=0;
MB_factor=1;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
%%
% [meas] = read_meas_dat_memmap_EPTI_GESE_SMSnoRef(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo,MB_factor);
% [meas] = read_meas_dat_memmap_EPTI_SMSnoRef_g(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo,MB_factor);

sTwixA = mapVBVD(filename,'removeOS','ignoreSeg','doAverage','rampSampRegrid');

sTwixA2 = mapVBVD(filename,'removeOS','doAverage','rampSampRegrid');

sTwixA3 = mapVBVD(filename,'removeOS','rampSampRegrid');

ADatax=sTwixA3{end};
asSlice=ADatax.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[ADatax.hdr.Phoenix.sSliceArray.asSlice{:}];
end
nSlices=numel(ADatax.hdr.Phoenix.sSliceArray.asSlice);
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

try
    FOVx=ADatax.hdr.Meas.ReadFOV;
catch
    FOVx=ADatax.hdr.Config.ReadFoV;
end

dx=RotatedLocs(2,1)/FOVx;
dy=RotatedLocs(1,1)/FOVx;

QQa=sTwixA2{2}.phasecor.unsorted();

QQb=sTwixA2{2}.image.unsorted();
%%
opt.ReturnStruct=1;
    opt.ReadMultipleRepetitions = 0;
    opt.SqueezeChannels = 1;
    opt.ReturnCellArray = 1; 
    
    % grab meas.data as a cell which make things a lot faster but some cells can be empty due to PAT and SMS

%     opt.SMSRefscan = 1;

    opt.SMSRefscan = 0;
    
    disp('read & save First Rep')
    tic
     
    %readin first rep and save prot,evp,(patrefscan,patrefscan_phascor)
    
    addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));

    meas_first = read_meas_dat(filename,opt);

    
PhaseCor3=sTwixA3{end}.phasecor();
PhaseCor3a=PhaseCor3(:,:,end,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

PhaseCor3b=cat(6,PhaseCor3a(:,:,:,:,:,1,:,:,:,:,1),PhaseCor3a(:,:,:,:,:,1,:,:,:,:,2),PhaseCor3a(:,:,:,:,:,2,:,:,:,:,1));
PhaseCor3b=perm32(perm63(PhaseCor3b));
PhaseCor3b(:,2,:,:,:,:,:,2,:,:,:,:)=PhaseCor3b(:,2,:,:,:,:,:,1,:,:,:,:);
PhaseCor3b(:,2,:,:,:,:,:,1,:,:,:,:)=PhaseCor3b(:,2,:,:,:,:,:,1,:,:,:,:)*0;
% size(meas.data)            240    85    32     1     1     1     1     2     1    32
% size(meas.data_phascor1d)  240     3    32     1     1     1     1     2     1    32
% 3×2 single matrix
% 
%     0.4579         0
%          0    0.4231
%     0.4018         0

    
PhaseCor=sTwixA{end}.phasecor();
Data=sTwixA{end}.image();
Data2=sTwixA2{end}.image();
Data2b=permute(Data2,[1 3 2 4:7 11 8:10 12:20]);

% linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d);
% 
% 
% linear_fit_coeff = mrir_artifact_ghost_compute(PhaseCor3b);

meas_Indiv=meas_first;
meas_Indiv.data=Data2b;
meas_Indiv.data_phascor1d=PhaseCor3b;
meas_Indiv.data=permute(Data2b,[1:6 9 8 7 10:20]);
meas_Indiv.data_phascor1d=permute(PhaseCor3b,[1:6 9 8 7 10:20]);
meas_Indiv.offline=[];
[meas_Indiv, linear_fit_coeff] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indiv,[],0);  % no meas_Collapsed
meas_Indiv.data = sum(meas_Indiv.data,8);

QQ=meas_Indiv.data;

QQc=QQ(:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

QQd=QQc.*exp(1i*linspaceWithHalfStep(-pi,pi,120)*120*dy);
IQQ=ifft2cg(QQd);

% QQa=grmss(QQ(:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[1 3]);

Sz=gsize(IQQ,1:2);
%%
try
    SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(IQQ(:,:,:,1)),0,Sz);
catch
    SelfSens=RunESPIRiTForSensMapsMultiMap(squeeze(IQQ(:,:,:,1)),28,Sz);
end
SelfSens1=SelfSens(:,:,:,1);
SelfSensP=perm43(SelfSens1);
%%
Rec1=bart('pics -m -S -R W:3:0:0.001 ',perm43(QQd(:,:,:,1)),SelfSensP);
Rec2=bart('pics -m -S -R W:3:0:0.001 ',perm43(QQd(:,:,:,2)),SelfSensP);
Rec3=bart('pics -m -S -R W:3:0:0.001 ',perm43(QQd(:,:,:,3)),SelfSensP);

Rec6=bart('pics -m -S -R W:3:0:0.00000001 ',perm43(Data6),SelfSensP);

%%
tmp=perm32(PartitionDim(QQa,3,25));
clear tmpa
tmpa(:,1,:,1,1,1,:,1)=tmp(:,1,:,:);
tmpa(:,3,:,1,1,1,:,1)=tmp(:,3,:,:);
tmpa(:,2,:,1,1,1,:,2)=tmp(:,2,:,:);
PhaseCorBase=tmpa;
%%
Idx=6;
for Idx=1:13
meas_Indivx=meas_Indiv;
CurPattern=Pattern(Idx,:);
meas_Indivx.data_phascor1d=PhaseCorBase(:,:,:,1,1,1,Idx,:);
CurData=QQb(:,:,60*(Idx-1)+(1:60));
tmp=zeros(120,120,32,1,1,1,1,2,1,1);
tmp(:,CurPattern(1:2:end),:,1,1,1,1,1)=perm32(CurData(:,:,1:2:end));
tmp(:,CurPattern(2:2:end),:,1,1,1,1,2)=perm32(CurData(:,:,2:2:end));
meas_Indivx.data=tmp;
[meas_Indivx, linear_fit_coeff] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indivx,[],0);  % no meas_Collapsed
CData = sum(meas_Indivx.data,8);
CData=CData.*exp(1i*linspaceWithHalfStep(-pi,pi,120)*120*dy);

CDataX(:,:,:,Idx)=CData;
end
% meas_Indiv.data_phascor1d 120     3    32     1     1     1    25     2
% data 120   122    32     1     1     1     1     2     1    25
%%
for Idx=1:13
    RecX(:,:,Idx)=bart('pics -m -S -R W:3:0:0.00000001 ',perm43(CDataX(:,:,:,Idx)),SelfSensP);
end
%%
FirstTE_ms=9;

ES_ms=ADatax.hdr.Config.EchoSpacing_us/1000;

ZZ=CDataX(:,60+(-24:25),:,1:5);
SelfSensPr=imresize(SelfSensP,[120 50]);

HW=hanning(size(ZZ,2)).';
Recr=squeeze(bart('pics -m -S ',HW.*perm83(perm43(ZZ)),SelfSensPr));

[~, UpdatedB0MapA, UpdatedT2SMap_msA, s_valsA, FittedA, PDBase0A]=...
    FitToModel_MPBD1CSf(Recr,1:5,ES_ms,FirstTE_ms);

B0x=imresize(UpdatedB0MapA,Sz);
% fgmontagex(B0x,[-100 100])
%%
save([OutP 'CurStatus1.mat']);
%%
nPatterns=size(Pattern,1);
nEchoes=size(Pattern,2);
PatternX=zeros([120,nEchoes,nPatterns]);
for i=1:nPatterns
    for j=1:nEchoes
        PatternX(Pattern(i,j),j,i)=1;
    end
end
%%
TEs_ms=(FirstTE_ms+(0:(nEchoes-1))*ES_ms);
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms3=TEs_ms3;
disp('ok');
%% T2*, TSC...
BaseFP=OutP;

T2svalues_ms=linspace(5,300,200);
Decays=exp(-TEs_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');

Ops={'fmac 0 1'};
FMScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots.txt';
WriteLinopToFile(FMScriptFN,Ops);

Ops_Subspace=['fmac 1 32',Ops];
FMScriptFN_Subspace='/autofs/space/daisy_002/users/Gilad/gUM/EPTI_Shots_Subspace.txt';
WriteLinopToFile(FMScriptFN_Subspace,Ops_Subspace);

ImSz16=FillOnesTo16([Sz 1 1 1 1 nEchoes]);

CompsPFN=[BaseFP 'CompsP'];

nComponentsToUse=2;
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);

writecfl(CompsPFN,CompsP);

ImSz16Comps=FillOnesTo16([Sz 1 1 1 nComponentsToUse 1]);
disp('Saved comps');
%%
nccToUse=size(SelfSensP,4);

PatternXP=permute(PatternX,[1 4 5 6 7 8 2 3]);

for i=1:nEchoes
    for s=1:nPatterns
        kLoc(i,s)=find(PatternXP(:,1,1,1,1,1,i,s));
    end
end
disp('got kLoc');
%%
for i=1:nPatterns
    for e=1:nEchoes
        sigA(:,Pattern(i,e),1,:,1,1,e,i)=CDataX(:,Pattern(i,e),:,i);
    end
end
CurSens=SelfSensP;

isigA=ifft1cg(sigA,1);
sisigPSA=sum(perm21(isigA),1);

% isigPS=zeros([Sz 1 nccToUse 1 1 nEchoes nShots]);
isigPSA=perm21(isigA);
MaskSamplesSA=abs(isigPSA)>0;
FMaskSamplesSA=fft1cg(single(MaskSamplesSA),1);

SensP=permute(CurSens(:,:,:,1:nccToUse),[2 1 3 4]);
%%
disp('a');
SigFN=[BaseFP 'sisigPS' ];
TSCPSensFMaskSamplesS_FN=[BaseFP 'TSCPSensFMaskSamplesS'];
SensFMaskSamplesS_FN=[BaseFP 'SensFMaskSamplesS'];

% WhichShots=1:5;
WhichShots=[4:7];
WhichShots=6:9;
% WhichShots=8:9;
WhichShots=10:13;
nShots=numel(WhichShots);
sig=sigA(:,:,:,:,:,:,:,WhichShots);

% sig kRO kPE 1 CC 1 1 Echos
isig=ifft1cg(sig,1);
sisigPS=sum(perm21(isig),1);

FMaskSamplesS=FMaskSamplesSA(:,:,:,:,:,:,:,WhichShots);

SensFMaskSamplesS=SensP.*FMaskSamplesS;

writecfl(SigFN,sisigPS);
writecfl(SensFMaskSamplesS_FN,SensFMaskSamplesS);

TSCP=exp(1i*2*pi*perm21(B0x).*TEs_ms3/1000);
TSCP7=perm73(TSCP);
TSCPSensFMaskSamplesS=TSCP7.*SensFMaskSamplesS;
writecfl(TSCPSensFMaskSamplesS_FN,TSCPSensFMaskSamplesS);

disp(['Saved for subspace']);
%
WhichEchosToUse=15:65;
    
RecStr=' -S ';
RecStr=' -S -R W:3:0:0.001 ';

RecSubspace=bart(['picsS  ' RecStr FMScriptFN_Subspace],ImSz16Comps,SigFN,TSCPSensFMaskSamplesS_FN,CompsPFN);
RecSubspaceX=sum(RecSubspace.*CompsP,6).*TSCP7;
disp(['done']);
%%
WhichEchoes=1:45;
WhichShots=10:13;
CDataX1=squeeze(sum(sigA(:,:,1,:,1,1,WhichEchoes,WhichShots),7));
% CDataX1=CDataX(:,:,:,10:13);
CDataX2=sum(CDataX1,4)./sum(abs(CDataX1)>0,4);
CDataX2(~isfinite(CDataX2))=0;

RecY=bart('pics -m -S -R W:3:0:0.00000001 ',perm43(CDataX2),SelfSensP);