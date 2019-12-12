setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TS')
%%
SensR=SensP(:,:,:,:,1);
% SensW=repmat(SensR,[1 1 1 1 2]);
%%
clc
Rec=bart(['pics -S -m -R T:7:0:' num2str(1e-6) ' -t'],BARTTrajAct, nukDataP, SensR);
%%
clc
SensW=repmat(SensR,      [1 1 1 1 1 3]);
DataW=repmat(nukDataP,   [1 1 1 1 1 3]);
TrajW=repmat(BARTTrajAct,[1 1 1 1 1 3]);
Rec=bart(['pics -S -m -R T:7:0:' num2str(1e-6) ' -t'],TrajW, DataW, SensW);
% fgmontage(Rec)
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TS')

nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
nukDataP=permute(nukData,[3 2 4 1]);
SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

% SensW=repmat(SensR,      [1 1 1 1 1 3]);
% DataW=repmat(nukDataP,   [1 1 1 1 1 3]);
% TrajW=repmat(BARTTrajAct,[1 1 1 1 1 3]);

%%
% 3TS is 6.564922 sec
% 10TS is 20.351443 sec
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TS')

nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
nukDataP=permute(nukData,[3 2 4 1]);
SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

% writecfl('/tmp/TSB',ones([1 5114 1 32 1 nTS])/nTS);
TSBFX=repmat(permute(TSBF,[6 2 5 4 3 1]),[1 1 1 nChannels 1]);
writecfl('/tmp/TSB',TSBFX);


clc
% SensW=repmat(SensP(:,:,:,:,1),      [1 1 1 1 1 nTS]);
SensW=SensP(:,:,:,:,1).*permute(TSC,[1 2 6 5 4 3]);
DataW=repmat(nukDataP,   [1 1 1 1 1 1]);
TrajW=repmat(BARTTrajAct,[1 1 1 1 1 nTS]);
% Rec=bart(['pics -S -m -R T:7:0:' num2str(1e-6) ' -t'],TrajW, DataW, SensW);
RecTS=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],TrajW, DataW, SensW);
%        ['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1)
%%
fgmontage(rot90(cat(3,im_res,RecS(:,:,SliI),RecTS)),'Size',[1 3])
title([PadStringWithBlanks('SparseMRI',100) PadStringWithBlanks('BART no B_0',100) 'BART TS'])
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03')

nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
nukDataP=permute(nukData,[3 2 4 1]);
SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);
    
RecA=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));
%%
% const struct linop_s* fft_op = nufft_create(DIMS, ksp_dims, coilim_dims, traj_dimsX, traj, weights, conf);
% ksp_dims    [    1  5114     1    32     1     3     1     1     1     1     1     1     1     1     1     1 ]
% coilim_dims [  128   128     1    32     1     3     1     1     1     1     1     1     1     1     1     1 ]
% traj_dimsX  [    3  5114     1     1     1     3     1     1     1     1     1     1     1     1     1     1 ]
% const struct linop_s* maps_op = mapsX2_create(coilim_dims, map_dims, img_dims, maps);
% coilim_dims [  128   128     1    32     1     3     1     1     1     1     1     1     1     1     1     1 ]
% img_dims    [  128   128     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% map_dims    [  128   128     1    32     1     3     1     1     1     1     1     1     1     1     1     1 ]
% const struct linop_s* TSB_op = mapsR2_create(TSB_dims, TSB_dims, TSB_dimsX, TSB);
% TSB_dims [    1  5114     1    32     1     3     1     1     1     1     1     1     1     1     1     1 ]
% TSB_dimsX [    1  5114     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% const struct linop_s* TSB_op = mapsY2_create(TSB_dimsX, TSB_dimsX, TSB_dims, TSB);
%%
% 2 maps.
% ESPIRiT reconstruction.
% sense_ncTS_init
% max_dimsX [  128   128     1    32     2     1     1     1     1     1     1     1     1     1     1     1 ]
% sense_nc_init
% coilim_dims [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% img_dims    [  128   128     1     1     2     1     1     1     1     1     1     1     1     1     1     1 ]
% map_dims    [  128   128     1    32     2     1     1     1     1     1     1     1     1     1     1     1 ]
% ksp_dims    [    1  5114     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% traj_dims   [    3  5114     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% max_dims    [  128   128     1    32     2     1     1     1     1     1     1     1     1     1     1     1 ]
%%
% coilim_dims [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% img_dims    [  128   128     1     1     2     1     1     1     1     1     1     1     1     1     1     1 ]
% map_dims    [  128   128     1    32     2     1     1     1     1     1     1     1     1     1     1     1 ]
% ksp_dims    [    1  5114     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% traj_dims   [    3  5114     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% max_dims    [  128   128     1    32     2     1     1     1     1     1     1     1     1     1     1     1 ]
% %%
% --maps_create_data:
% ksp_dims    [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% img_dims    [  128   128     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% mps_dims    [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
%%
% --linop_fmac_create:
% odims    [  128   128     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% idims    [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% tdims    [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
%%
% sense_ncTS_init
% max_dimsX [  128   128     1    32     1     2     1     1     1     1     1     1     1     1     1     1 ]
% sense_nc_init
% coilim_dims [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% img_dims    [  128   128     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% map_dims    [  128   128     1    32     1     1     1     1     1     1     1     1     1     1     1     1 ]
% ksp_dims    [    1  5114     1    32     1     2     1     1     1     1     1     1     1     1     1     1 ]
% traj_dims   [    3  5114     1     1     1     1     1     1     1     1     1     1     1     1     1     1 ]
% max_dims    [  128   128     1    32     1     2     1     1     1     1     1     1     1     1     1     1 ]