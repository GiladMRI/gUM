EPTIRecP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Recon/Human/';

G6=load([EPTIRecP 'Recon_EPTI_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_6_GE.mat']);
G16=load([EPTIRecP 'Recon_EPTI_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_16_GE.mat']);
S6=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_6_GE.mat']);
S16=load([EPTIRecP 'Recon_EPTI_Subspace_Human_SMS1_1p9_3SHOT_GE_Dyn2_Slice_16_GE.mat']);
%%
ShowAbsAngle(G6.im_EPTI_correct(:,:,1,:))
%%
fgmontage(G6.im_EPTI_correct(:,:,1,11:11:end))
fgmontage(G16.im_EPTI_correct(:,:,1,11:11:end))

fgmontage(gflip(G6.im_EPTI_correct(:,:,11:11:end),2))
%%
fgmontage(S6.im_recon(:,:,11:11:end))
fgmontage(S16.im_recon(:,:,11:11:end))
%%
fgmontage(gflip(S6.im_recon(:,:,11:11:end),2))
fgmontagex(gflip(S16.im_recon(:,:,1:21:end),2),[0 6.5]);title('EPTI subspace 3-shot');
%%
CalibP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/data/Data_acq/';
C6=load([CalibP 'Human_Calib_1p9_Slice6.mat']);
C6g=load([CalibP 'Human_Calib_1p9_Slice16_generated.mat']);
%%
EBaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
R=load([EBaseP 'meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib_FirstRep_Raw.mat']);
%%
ScanP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
BaseFN='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
MIDStr=BaseFN(6:11);
FN=[ScanP BaseFN '.dat'];
disp('ok');
%%
mainP=[ScanP BaseFN];
mkdir(mainP);

system(['chmod +777 -R ' mainP]);
disp([mainP ' Created']);
%%
AData = mapVBVD(FN);
if(iscell(AData))
    ADatax=AData{end};
else
    ADatax=AData;
end
%%
ADataI=ADatax.image(); % 240    32    27     1    32     1     1     1    39     1     2
%%
QQ=squeeze(ADataI(:,12,3:end-2,1,23,1,1,1,23,1,:));
QQA=sum(QQ,3);
%%
%% demo for SMS PSF processing
% Loop by slice since the dataset could be very large
clear;
close all;
%%
BaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
FN='meas_MID04688_FID21650_ep2d_ge_EPTI_1p3_7shot_calib';
% FN='meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns';

mkdir([BaseP FN]);
%%
directory = '/autofs/cluster/kawin/FuyixueWang/EPTI/EPTI_Rec_pack/';
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
%% Get Parameters
filename = [directory_rawdata,FN,'.dat'];
save_filename = 'Human_Calib_1p3_7shot';
% save_filename = 'Human_1p9_3shot';

% read parameters
pf_echo=0;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo);
Nreps = meas.evp.RawRep;
Necho=meas.evp.NEcoMeas;
parameters=meas.prot;
clear meas.data;
nslice_group=meas.evp.NSlcMeas;
generate_calib = 1;


for slice= 1:nslice_group
% for slice= 6
    %% Process of Refscan
    RepsToRead=1:Nreps;
%     kdata  = EPTI_SMS_Preprocess_Imgscan_memorySave_GESE(filename,1,RepsToRead,SMS_data,slice,pf_echo);
    kdata  = EPTI_SMS_Preprocess_Imgscan_memorySave_GESE(filename,0,RepsToRead,SMS_data,slice,pf_echo);
    kdataS(:,:,:,:,slice)=kdata;
end
%%
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'))
%%
nSlices=nslice_group;
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);
disp('ok');
%%
% save([BaseP FN filesep 'Locs.mat'],'RotatedLocs');
% disp(['Saved ' BaseP FN filesep 'Locs.mat']);
%%
% D=gpermute(Data,[1 3 4 2 8 5 6 7]);
D=permute(kdataS(3:end-2,42+(1:84),:,:,:),[2 3 4 1 6 5]);
% PD=padLeft(D,24,2);
% D=D(:,:,:,:,:,ROrd);
% PD=RepDotMult(PD,permute( mod(1:nSlices,2)*2-1,[1 6 3 4 5 2]));
I=squeeze(ifft2cg(D));

EchoSpace_us=parameters.iEffectiveEpiEchoSpacing;
% TEs_us=[sTwixX.hdr.Phoenix.alTE{:}];

nEchos=size(I,4);
TEs_us=2e3+(0:nEchos-1)*EchoSpace_us;
TEs_us=TEs_us(1:nEchos);
TEs_ms=TEs_us/1000;
dTEs_ms=diff(TEs_ms,1,2);

FirstEcho=squeeze(I(:,:,:,1,:));

SnsSzB=gsize(I,1:2);

disp('ok');
%%
for SliI=1:nSlices
    disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
    SensBMM{SliI}=RunESPIRiTForSensMapsMultiMap(FirstEcho(:,:,:,SliI),0,SnsSzB);
%     SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);1
end
SensB=permute(cat(5,SensBMM{:}),[1 2 3 5 4]);

SensMsk=grmss(SensB(:,:,:,:,1),3)>0.01;
%%
save([BaseP FN filesep 'Sens.mat'],'SensB');
disp(['Saved ' BaseP FN filesep 'Sens.mat']);
%%
for SliI=1:nSlices
    figure;
    subplot(2,2,1);gmontage(abs(SensB(:,:,:,SliI,1)),[0 0.7]);
    subplot(2,2,2);gmontage(angle(SensB(:,:,:,SliI,1)),[-pi pi]);
    subplot(2,2,3);gmontage(abs(SensB(:,:,:,SliI,2)),[0 0.7]);
    subplot(2,2,4);gmontage(angle(SensB(:,:,:,SliI,2)),[-pi pi]);
%     ShowAbsAngle(SensB(:,:,:,SliI,:))
    YLbl=['Sli' num2str(SliI,'%02d')];
    ylabel(YLbl);
    gprint(get(gcf,'Number'),[BaseP FN filesep 'Sens_' YLbl],[]) 
    close(gcf);
end
disp('printed sens images');
%%
% SensB [X Y Ch Sli Maps]
% I [X Y Ch Echo Slices]
Combined=gpermute(sum(I.*gpermute(conj(SensB(:,:,:,:,1)),[5 4]),3),[5 3]);
%%
WhichEchosToUse=1:nEchos;
for i=1:nSlices
    disp(i);
    [PDBase(:,:,i), UpdatedB0Map_Hz(:,:,i), UpdatedT2SMap_ms(:,:,i), s_vals(:,:,:,i), Fitted0(:,:,:,i), PDBase0(:,:,i)]=FitToModel_MPBD1CSf(Combined(:,:,i,:),WhichEchosToUse,dTEs_ms(1),TEs_ms(1));
end
%%
fgmontagex(abs(UpdatedT2SMap_ms(:,:,3:8:end)),[0 100]);title('EPTI-Calib T_2^*');
fgmontagex(UpdatedB0Map_Hz(:,:,3:8:end),[-300 300]);title(['EPTI-Calib B_0, TEs: ' num2str(TEs_ms,'%.1f ')]);
ShowAbsAngle(PDBase0(:,:,3:8:end),1,'Size',[2 2])
fgmontagex(s_vals(:,:,:,3:8:end));title('EPTI-Calib SV maps');
%%
PDBase0x=min(PDBase0,6*grmss(PDBase0));
[Out B1 BN1]=CalcSlicesSNR(abs(PDBase0x(:,:,:)),false,5);
B2=(~BN1).*SensMsk;
B2D=imdilate(B2,strel('disk',3,8));
B3=imfillholesBySlices( B2D );
for i=1:nSlices
    B4(:,:,i)=getLargestComponent(B3(:,:,i));
end
B4=B4.*SensMsk;
%%
dAngle=UpdatedB0Map_Hz*2*pi*dTEs_ms(1)/1000;
[unwrapped] = cusackUnwrap(dAngle, grmss(Combined,4));
unwrapped=unwrapped.*B4;
B0M_Hz=unwrapped*1000/2/pi/dTEs_ms(1);
fgmontagex(B0M_Hz(:,:,1:16),[-300 300]);title(['EPTI-Calib B_0 unwrapped, ES: ' num2str(EchoSpace_us,'%.1f ')]);
fgmontagex(UpdatedB0Map_Hz(:,:,1:16),[-300 300]);title(['EPTI-Calib B_0, ES: ' num2str(EchoSpace_us,'%.1f ')]);
%%
save([BaseP FN filesep 'B0T2S.mat'],'SensMsk','B1','BN1','B2','B2D','B3','B4',...
    'B0M_Hz','UpdatedB0Map_Hz','UpdatedT2SMap_ms','s_vals','PDBase0','TEs_ms');