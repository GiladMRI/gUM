clear;
close all;
%%
% directory = '/autofs/cluster/kawin/FuyixueWang/EPTI/EPTI_Rec_pack/';
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';

addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
%% Get Parameters
FNBase='meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns';
FNBaseCalib='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
filename_calib = 'Human_Calib_1p3_5shot';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);

CalibP=[directory_rawdata FNBaseCalib filesep];
filename = [directory_rawdata,FN,'.dat'];

apodization_para=0;

% read parameters
pf_echo=0;
MB_factor=1;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_GESE_SMSnoRef(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo,MB_factor);

Nseg = meas.prot.sWipMemBlock_alFree(34);
Nreps = Nseg;
Ndyn = meas.evp.RawRep/Nreps;
Necho=meas.evp.NEcoMeas;
parameters=meas.prot;
clear meas.data;
Rseg=meas.prot.sWipMemBlock_alFree(33);
Rpe=meas.prot.sWipMemBlock_alFree(35);
Ngroup=floor(parameters.lPhaseEncodingLines/(Rseg/Rpe));
MB=1;
param.Nseg=Nseg;
param.Rseg=Rseg;
param.Rpe=Rpe;
param.Ngroup=Ngroup;
parameters.nechoGE=meas.evp.NLinMeas;
%%
% parameters.nechoGE=size(kdata,1);
save([OutP,'parameters.mat'],'parameters');
%% Slice Loop
nslice_group=meas.evp.NSlcMeas/MB;
slice_sep=nslice_group;

for dyn=2:2

    RepsToRead=1:Nreps;
    RepsToRead=RepsToRead+(dyn-1)*Nreps;
    [kdata_GE_all] = EPTI_SMS_Preprocess_Imgscan_SMS1_GESE(filename,1,RepsToRead,SMS_data,1:nslice_group,pf_echo);
%% reconstruction loop of all slices
%     k_recon=zeros(meas.evp.NImageLins,Nseg*Rseg,nslice_group*MB,Ngroup*(Rseg/Rpe),meas.evp.NChaMeas);
%     kdata_correct=zeros(size(k_recon));
    for slice=1:nslice_group
        disp(['---slice ---- ',num2str(slice)]);
        tic;
        kdata_GE=kdata_GE_all(:,:,:,:,slice);
        [kdata_GE] = putRawtoPattern(kdata_GE,Nseg,Rseg,Rpe,Ngroup,'GE',MB);
        kdata_GES(:,:,:,:,slice)=kdata_GE;
    end
end
%%
save([OutP 'kdata_GES.mat'],'kdata_GES','-v7.3');
%%
nRO=size(kdata_GE,2);
nPE=size(kdata_GE,3);
nEchos=size(kdata_GE,1);
nChannels=size(kdata_GE,4);
nSlices=nslice_group;
nShots=Nseg;
%% Fuyixue GRAPPA "Script1"
for slice=1:nSlices
    disp(slice);
    disp(datestr(now));
    kdata_GE=kdata_GES(:,:,:,:,slice);
    
    k=squeeze(sum(kdata_GE(1:16,:,:,:),1));
    %         figure; imshow(permute(sos(ifft2c(k)),[2 1]),[0 5]); colormap('gray');
    
    echotype='GE';
    %         load([directory,'data/Data_acq/',filename_calib,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib');
    load([CalibP,filename_calib,'_Slice',num2str(slice),'.mat'],'kdata_calib');
    kdata_calib=double(kdata_calib);
    
    [recon] = EPTI_recon_SMS1_GenCalib(kdata_GE,kdata_calib,echotype,param);
%     close all;
%     im_recon=ifft2c(recon);
%     im_recon_show = sos(im_recon(:,:,13:4:end-14,:),4);
%     figure; imshow3(im_recon_show,[0 8],[2,size(im_recon_show,3)/2]);
    k_recon(:,:,1,:,:)=single(recon);
    [tmp,B0_variation_esti] = B0_variationCorrection(recon,param,parameters,0);
    kdata_correct(:,:,1,:,:) = single(tmp);
    disp(datestr(now));
%     reconS(:,:,:,:,slice)=recon;
%     B0_variation_estiS(:,:,:,slice)=B0_variation_esti;
%     kdata_correctS(:,:,:,:,slice)=single(tmp);

    disp(['Post Processing Start: apodization and coil combination']);
    disp(datestr(now));
    [im_EPTI] = fMRI_PostRecon_process(k_recon,apodization_para);
    [im_EPTI_correct] = fMRI_PostRecon_process(kdata_correct,apodization_para);
    disp(datestr(now));

    save([OutP,'Recon_EPTI','_Dyn',num2str(dyn),'_Slice_',num2str(slice),'_',echotype,'.mat'],'im_EPTI','im_EPTI_correct',...
        'recon','kdata_correct','B0_variation_esti','-v7.3');
end
%%
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%%
load([OutP,'parameters.mat'],'parameters');
echotype='GE';
dyn=2;
slice=5;
load([OutP,'Recon_EPTI','_Dyn',num2str(dyn),'_Slice_',num2str(slice),'_',echotype,'.mat']);
[nRO,nPE,nEchos,nChannels]=size(recon);
ES_ms=parameters.iEffectiveEpiEchoSpacing/1e3;
FirstTE_ms=2;
% dThickness: 92
% dPhaseFOV: 228
% dReadoutFOV: 228
%%
IdxForSens=20;
im_recon_forSens=squeeze(ifft2c(recon(:,:,IdxForSens,:)));
disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
CurSens=RunESPIRiTForSensMapsMultiMap(im_recon_forSens,0,[nRO,nPE]);
CurSens=CurSens(:,:,:,1);
SensMsk=grmss(CurSens,3)>0.01;
%%
RecAllEchosChannels=ifft2c(recon);
RecAllEchos=sum(RecAllEchosChannels.*conj(perm43(CurSens)),4);
%%
WhichEchosToUse=15:55;
% WhichEchosToUse=10:14;
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(RecAllEchos,WhichEchosToUse,ES_ms,FirstTE_ms);
%%




