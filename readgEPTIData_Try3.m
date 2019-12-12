clear;
close all;
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%%
directory_rawdata = '/autofs/space/daisy_001/users/data/Gilad/SEPTI_2Echo_Try1_1Dec19/';
%% Get Parameters

FNBase='meas_MID01247_FID53829_gep2d_327_P1_1Sli_16Rep';
FNBase='meas_MID01249_FID53831_gep2d_327_P2_1Sli_16Rep';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);

filename = [directory_rawdata,FNBase,'.dat'];
%
apodization_para=0;

pf_echo=0;
MB_factor=1;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
%%
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));

% sTwixA = mapVBVD(filename,'removeOS','ignoreSeg','doAverage','rampSampRegrid');
% 
% sTwixA2 = mapVBVD(filename,'removeOS','doAverage','rampSampRegrid');

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

% QQa=sTwixA2{2}.phasecor.unsorted();
% 
% QQb=sTwixA2{2}.image.unsorted();
disp('Read data');
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
disp('Prepared a');
    %%
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

    
% PhaseCor=sTwixA{end}.phasecor();
% Data=sTwixA{end}.image();
% Data2=sTwixA2{end}.image();
try
%     Data2=sTwixA2{end}.image(:,:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:);
    Data2=sTwixA3{end}.image(:,:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:);
catch
%     Data2=sTwixA2{end}.image(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
    Data2=sTwixA3{end}.image(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
end
Data2b=permute(Data2,[1 3 2 4:7 11 8:10 12:20]);

% linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d);
% 
% 
% linear_fit_coeff = mrir_artifact_ghost_compute(PhaseCor3b);

meas_Indiv=meas_first;
meas_Indiv.data=Data2b;
meas_Indiv.data_phascor1d=PhaseCor3b;
% meas_Indiv.data=permute(Data2b,[1:6 9 8 7 10:20]);
meas_Indiv.data=permute(Data2b,[1:6 10 8 7 9 11:20]);
meas_Indiv.data_phascor1d=permute(PhaseCor3b,[1:6 9 8 7 10:20]);
meas_Indiv.offline=[];
[meas_Indiv, linear_fit_coeff] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indiv,[],0);  % no meas_Collapsed
meas_Indiv.data = sum(meas_Indiv.data,8);

QQ=meas_Indiv.data;
if(size(QQ,2)<120)
    QQ(end,120,end,1)=0;
end

QQc=QQ(:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);

QQd=QQc.*exp(1i*linspaceWithHalfStep(-pi,pi,120)*120*dy);
IQQ=ifft2cg(QQd);
disp('Finished base phase cor');

% QQa=grmss(QQ(:,1:120,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),[1 3]);

Sz=gsize(IQQ,1:2);
%%
A=squeeze(Data2b); % 1 : 22:81
AA=squeeze(IQQ);

figure;plot(grmss(A(:,:,:,:,1:14),[1 3 4]))

figure;plot(grmss(QQc(:,:,:,:,1:14),[1 3:20]))

% ZZ=QQc(:,60-3+(-11:11),:,15);IZZ=ifft2cg(ZZ);ShowAbsAngle(IZZ)
ZZ=QQc(:,60-3+(-11:10),:,15);IZZ=ifft2cg(ZZ);ShowAbsAngle(IZZ)
%%
% ForSelfSens=squeeze(IQQ(:,:,:,12));
ForSelfSens=IZZ;
try
    SelfSens=RunESPIRiTForSensMapsMultiMap(ForSelfSens,0,Sz);
catch
    SelfSens=RunESPIRiTForSensMapsMultiMap(ForSelfSens,20,Sz);
end
SelfSens1=SelfSens(:,:,:,1);
SelfSensP=perm43(SelfSens1);
%%
tmp=squeeze(Data2b(:,:,:,1,1,1,1,:,1,1));
seg1=tmp(:,:,:,1);
seg2=tmp(:,:,:,2);
Idx1=find(grmss(seg1,[1 3])>0);
Idx2=find(grmss(seg2,[1 3])>0);

clear tmpData
tmpData(:,(1:numel(Idx1))*2-1,:)=seg1(:,Idx1,:);
tmpData(:,(1:numel(Idx2))*2-1,:,1,1,1,1,2)=seg2(:,Idx2,:);
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