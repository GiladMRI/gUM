clear;
close all;
clc;
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI_3D/funcs'));
% addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI_3D/Data3DEPTI_1018/pattern_test'));
addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/pattern_test')
addpath('/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/Recon_Test_code/');
%% load data
% directory = '/autofs/cluster/kawin/FuyixueWang/EPTI_3D/Data3DEPTI_1018/';
directory = '/autofs/cluster/kawin/Zijing/Data3DEPTI_1018/';
RO_select = 110;
filename='Fully_Sample';
load([directory,'data/Data_acq/meas_prot_',filename,'.mat']);
load([directory,'Recon_data/Basis_Only_GE_PSIF','.mat']);
% load([directory,'Recon_data/Basis_Only_RealData','.mat']);

filename=['Fully_Sample_RO',num2str(RO_select)];
load([directory,'data/Data_acq/',filename,'.mat']);

%% Set Parameters
K=4;
Phi = U(:,1:K);
dt = param.dt; % echo spacing
t0 = param.t0; % % time for first echo
nechoGE=size(kdata_full,1);
%% Undersampling and Generate Calibration data
Block_size_y=12;
Block_size_z=6;
mask_sample = EPTI_sampling_mask_Regular_YZ_V83(kdata_full,Block_size_y,Block_size_z);
mask_sample2 = circshift(mask_sample,[0,-2,-1,0]);
mask_sample2 = mask_sample+mask_sample2;
mask_sample = mask_sample2;
kdata=mask_sample.*kdata_full;
[nt,nx,npe,nc]=size(kdata);

numberOFecho = 6;
calib_PEsize=[48,48];
kdata_calib=kdata_full(1:numberOFecho,:,:,:);
kdata_calib=crop(kdata_calib,[numberOFecho,calib_PEsize,size(kdata_calib,4)]);
t0_calib = t0; % time before first echo

TEs_GRE_calib=(0:numberOFecho-1)*dt+t0_calib;
TEs_GRE_calib = TEs_GRE_calib(:);
TEs_GRE=(0:nt-1)*dt+t0;
TEs_GRE=TEs_GRE(:);
%% Pre-Processing: smap, background phase, B0-inhomogeneity phase
[ sens_map,P_dB,Phase0 ] = Preprocess_SubspaceRec_Use( kdata_calib,TEs_GRE_calib,[nx,npe] );
figure; imshow(abs(permute(P_dB(:,:,1),[2 1])),[-100 100]);
figure; imshow((permute(Phase0,[2 1])),[-pi pi]);
sens_map=squeeze(sens_map);
Phase0=exp(1i*Phase0);
% figure; imshow(angle(permute(Phase0.*(-P_dB*dt),[2 1])),[-pi pi]);
% figure; imshow(angle(permute(Phase0.*(-P_dB*dt),[2 1])),[-pi pi]);

Phase_T=zeros(nx,npe,nt);
for t=1:nechoGE
    Phase_T(:,:,t)=exp(1i*2*pi*P_dB*TEs_GRE(t)).*Phase0;  
end
%% Reconstruction
a0=zeros([nx,npe,K]);
disp('Reconstruction Start');
tic;
iter_ops.max_iter = 80;
llr_ops.lambda = .01;
iter_ops.rho = 0.1;
llr_ops.block_dim = [12, 12];
lsqr_ops.max_iter = 5;
lsqr_ops.tol = 1e-4;

[im_recon,a,history]=EPTI_Image_Recon_Subspace_GESE_ADMM(kdata,mask_sample,sens_map,Phase_T,Phi,iter_ops,llr_ops,lsqr_ops);
toc;
disp('Reconstruction Done');

%%

im_recon=permute(abs(im_recon),[2 1 3]);
im_full_ref=permute(sos(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);
im_full_phase=permute(sop(ifft2c(permute(kdata_full,[2,3,1,4])),4),[2 1 3]);

figure; imshow3(permute(abs(a),[2 1 3]),[],[1 size(a,3)]);
figure; imshow3(sos(im_recon(:,:,10:end-10),3),[],[1 1]);

%% 
T_to_show=[10,20,30,40];
im_recon_show = permute((im_recon(:,:,T_to_show,:)),[1,2,3]);
im_ref_show = permute((im_full_ref(:,:,T_to_show)),[1,2,3]);
im_diff = 5*abs(im_recon_show-im_ref_show);

im_to_show=cat(3,im_recon_show,im_ref_show,im_diff);
% im_to_show=crop(im_to_show,[200,180,size(im_to_show,3)]);
figure; imshow3(im_to_show,[0 3],[3,length(T_to_show)]);


RMSE_total=rmse3d(im_recon,im_full_ref,im_full_ref>0.8*mean(im_full_ref(:)));
disp(['RMSE_total is ',num2str(RMSE_total)]);
figure; imshow(im_diff(:,:,2),[0 3]);


%% single pixel
figure; 
% x_pos=82; y_pos=95; 
x_pos=50; y_pos=91; 

signal1 = squeeze(im_recon(x_pos,y_pos,:));
signal2=squeeze(im_full_ref(x_pos,y_pos,:)); 
plot(signal1,'r','LineWidth',2); hold on; plot(signal2,'b--','LineWidth',2);
h=legend('subspace','full reference');set(h,'FontSize',18);
%%
%% T2 T2* fitting
index_t=10:40;
TEs_GRE_fit=TEs_GRE(index_t);
im_for_fit=im_recon(:,:,index_t);
% im_for_fit=im_full_ref(:,:,index_t);

im_sos=sos(im_for_fit,3);
mask=im_sos>mean(im_sos(:));
threshold=0;
[S0_EPTI, T2s_EPTI, fitting_map_EPTI] = T2s_GREEPI(im_for_fit,TEs_GRE_fit,mask,threshold);

T2s_EPTI=T2s_EPTI*1e3;
T2s_EPTI(T2s_EPTI<0)=1000;

protonDensity1 = S0_EPTI;
[protonDensity1,~] = nomalization_forImage(protonDensity1,4);
% mask_forProton = mask;
% mask_forProton(60:140,80:130) = double(protonDensity1(60:140,80:130)<2);
% figure,imagesc(mask_forProton.*protonDensity1,[0,1.7]);colormap('gray');

fig1=figure; 
imagesc(T2s_EPTI,[0 250]); colormap('hot'); colorbar;h=title('T2s');set(h,'FontSize',18);
%%
apx='_48_48_6Echo_CalibS2';
save([directory,'/Recon_data/test/Recon_LLR_VDS_RO_',num2str(RO_select),'_K',num2str(K),apx,'.mat'],'im_recon','a','T2s_EPTI');