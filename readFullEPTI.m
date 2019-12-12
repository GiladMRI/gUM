%% demo for SMS PSF processing
% Loop by slice since the dataset could be very large
clear;
close all;
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
FNBase='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';

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
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%% Get Parameters
% filename = [directory_rawdata,'meas_MID04688_FID21650_ep2d_ge_EPTI_1p3_7shot_calib.dat'];
% save_filename = 'Human_Calib_1p3_7shot';

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
%%
RepsToRead=1:Nreps;
for slice= 1:nslice_group
    delete([OutP 'kdata_Slice' num2str(slice) '.mat']);
    kdata  = EPTI_SMS_Preprocess_Imgscan_memorySave_GESE(filename,1,RepsToRead,SMS_data,slice,pf_echo);
    kdata=single(permute(kdata,[2 3 1 4]));
    save([OutP 'kdata_Slice' num2str(slice) '.mat'],'kdata','-v7.3');
    disp(['saved ' OutP 'kdata_Slice' num2str(slice) '.mat']);
end
%% Now sens
for slice= 1:nslice_group
    load([OutP 'kdata_Slice' num2str(slice) '.mat'],'kdata');
end
%%