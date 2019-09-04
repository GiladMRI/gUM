%% demo for SMS PSF processing
% Loop by slice since the dataset could be very large
clear;
close all;
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
FNBase='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
%% Get Parameters
% filename = [directory_rawdata,'meas_MID04688_FID21650_ep2d_ge_EPTI_1p3_7shot_calib.dat'];
% save_filename = 'Human_Calib_1p3_7shot';
filename = [directory_rawdata,FNBase,'.dat'];
save_filename = 'Human_Calib_1p9';

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
    kdata  = EPTI_SMS_Preprocess_Imgscan_memorySave_GESE(filename,1,RepsToRead,SMS_data,slice,pf_echo);
    kdata=kdata(3:3+22,:,:,:,:);
    kdata_calib=single(kdata);
    
%     save([OutP,'data/Data_acq/',save_filename,'_Slice',num2str(slice),'.mat'],'kdata_calib','-v7.3');
    save([OutP,save_filename,'_Slice',num2str(slice),'.mat'],'kdata_calib','-v7.3');
    % generating calibration data
    if(generate_calib == 1)
        sData=size(kdata);
        kdata=crop(kdata,sData);
        [kdata_calib,B0fit] = generate_Calib(kdata,parameters,0.3);
        kdata_calib = single(kdata_calib);

        k=permute(kdata_calib,[2 3 1 4]);
        im=ifft2c(k);
%         figure; imshow(permute(sos(im(:,end:-1:1,14,:),4),[2 1 3]),[]);

%         save([OutP,'data/Data_acq/',save_filename,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib','B0fit','-v7.3');
        save([OutP,save_filename,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib','B0fit','-v7.3');
    end
end
%%
parameters.nechoGE=size(kdata,1);
save([OutP,'meas_prot_',save_filename,'.mat'],'parameters');
%%
k=permute(kdata_calib,[2 3 1 4]);
kExt=padarray(padarray(k,[0 40 0 0],'post'),[0 41 0 0],'pre');
idata_calib=ifft2c(k);
iExt=ifft2c(kExt);

CurSensC=RunESPIRiTForSensMapsMultiMap(squeeze(iExt(:,:,1,:)),0,[nRO,nPE]);
CurSensC=CurSensC(:,:,:,1);
SensMskC=grmss(CurSensC,3)>0.01;
%%
iExtC=sum(iExt.*conj(perm43(CurSensC)),4);
%%
WhichEchosToUse=3:23;
% WhichEchosToUse=10:14;
[PDBaseC, UpdatedB0Map_HzC, UpdatedT2SMap_msC, s_valsC, Fitted0C, PDBase0C]=...
    FitToModel_MPBD1CSf(iExtC,WhichEchosToUse,ES_ms,FirstTE_ms);

%%
idata_calib=ifft2c(kdata_calib);

CurSens=RunESPIRiTForSensMapsMultiMap(im_recon_forSens,0,[nRO,nPE]);
CurSens=CurSens(:,:,:,1);
SensMsk=grmss(CurSens,3)>0.01;
