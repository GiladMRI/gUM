load('/autofs/cluster/kawin/FuyixueWang/Share/Fully_sampled_2DGESE/Full_GESE_3mm_Slice13.mat');
parameters=load('/autofs/cluster/kawin/FuyixueWang/Share/Fully_sampled_2DGESE/meas_prot_Full_GESE_3mm.mat');
parameters=parameters.parameters;
ES_ms=parameters.iEffectiveEpiEchoSpacing/1000;
%%
size(kdata_GE)

size(kdata_SE)
%%
idata_GE=fft2cg(permute(kdata_GE,[2 3 4 1]));
idata_SE=fft2cg(permute(kdata_SE,[2 3 4 1]));

idata_GErms=grmss(idata_GE,3);
idata_SErms=grmss(idata_SE,3);
%%
WhichEchosToUseGE=3:37;
[PDBase_GE, UpdatedB0Map_GE, UpdatedT2SMap_ms_GE, s_vals_GE, Fitted0_GE, PDBase0_GE]=...
    FitToModel_MPBD1CSf(idata_GErms,WhichEchosToUseGE,ES_ms,0);

WhichEchosToUseSEa=3:37;
[PDBase_SEa, UpdatedB0Map_SEa, UpdatedT2SMap_ms_SEa, s_vals_SEa, Fitted0_SEa, PDBase0_SEa]=...
    FitToModel_MPBD1CSf(idata_SErms,WhichEchosToUseSEa,ES_ms,0);

WhichEchosToUseSEb=43:77;
[PDBase_SEb, UpdatedB0Map_SEb, UpdatedT2SMap_ms_SEb, s_vals_SEb, Fitted0_SEb, PDBase0_SEb]=...
    FitToModel_MPBD1CSf(idata_SErms,WhichEchosToUseSEb,ES_ms,0);

R2_minus_R2p=1./UpdatedT2SMap_ms_SEa;
R2_plus_R2p=1./UpdatedT2SMap_ms_SEb;
R2=(R2_plus_R2p+R2_minus_R2p)/2;
R2p=(R2_plus_R2p-R2_minus_R2p)/2;
T2Est=1./R2;
T2pEst=1./R2p;

fgmontagex(UpdatedT2SMap_ms_GE,[0 250]);colormap hot

fgmontagex(UpdatedT2SMap_ms_SEb,[0 250]);colormap hot;

fgmontagex(T2Est,[0 250]);colormap hot;
%%
SE_fitting_secondHalf = 0;

SE_fitting_secondHalf = 1; % ggg

shift_TE=0;  % off-center for TE
apodization_para=0.25; %% fliter/method: tukey window
Nedge=14; %points to crop

Nedge=0; % ggg

mask_thr=1.0;
thr=0;

nt=size(idata_GErms,3); % ggg


nGEs=nt;%parameters.nechoGE %mjf??
center_index = (nGEs+1)/2+shift_TE;
Effect_Tesp=parameters.iEffectiveEpiEchoSpacing/1000; %1.05abstract %unit (ms)
TE_GRE=parameters.alTE(1)/1000;    %25 abstract                          %unit (ms)
TEs_fast=(1:nt)*Effect_Tesp+TE_GRE-center_index*Effect_Tesp;
TEs_fast=TEs_fast(:);

TEs_GRE=crop(TEs_fast,[nt-Nedge,1]); 




nt=size(idata_SErms,3); % ggg

% recon=double(recon);
% Get TEs
% [npe,nro,nt,~]=size(recon);
nSEs=nt;%parameters.nechoSE; %mjf??
center_index = (nSEs+1)/2+shift_TE;
Effect_Tesp=parameters.iEffectiveEpiEchoSpacing/1000; %1.05 abstract %unit (ms)
TE_SE=parameters.alTE(2)/1000;       %98 abstract                       %unit (ms)
TEs_fast=(1:nt)*Effect_Tesp+TE_SE-center_index*Effect_Tesp;
TEs_fast=TEs_fast(:);


TEs_SE=crop(TEs_fast,[nt-Nedge,1]); 

%%
im_GRE=idata_GErms;
im_SE=idata_SErms;
%% Fitting
im_mean=mean(im_GRE,3);
mask=im_mean>mask_thr*mean(im_mean(:));
threshold=thr*mean(im_mean(:));

% T_180RF=TE_SE/2;
im_EPTI=cat(3,im_GRE,im_SE);
disp('start Fitting!');

[S0_EPTI,S02_EPTI,T2_EPTI,T2s_EPTI,fitting_map_GRE,fitting_map_SE] = gT2T2s_GE_2SE(im_EPTI,TEs_GRE,TEs_SE,TE_SE,mask,threshold);

% if(SE_fitting_secondHalf==0)
%     [S0_EPTI,S02_EPTI,T2_EPTI,T2s_EPTI,fitting_map_GRE,fitting_map_SE] = T2T2s_GE_SE(im_EPTI,TEs_GRE,TEs_SE,TE_SE,mask,threshold);
% else
%     [S0_EPTI,S02_EPTI,T2_EPTI,T2s_EPTI,fitting_map_GRE,fitting_map_SE] = T2T2s_GE_2SE(im_EPTI,TEs_GRE,TEs_SE,TE_SE,mask,threshold,floor(center_index)-Nedge/2);
% end

fitting_map_EPTI=cat(3,fitting_map_GRE,fitting_map_SE);
TEs_EPTI=cat(1,TEs_GRE,TEs_SE);
T2s_EPTI(T2s_EPTI<0)=1000;
T2_EPTI(T2_EPTI<0)=1000;

protonDensity1 = S0_EPTI;
[protonDensity1,~] = nomalization_forImage(protonDensity1,3 );
mask_forProton = mask;
mask_forProton(60:140,80:130) = double(protonDensity1(60:140,80:130)<6);
figure,imagesc(abs(mask_forProton.*protonDensity1),[0,1.8]);colormap('gray');

protonDensity2 = S02_EPTI;
[protonDensity2,~] = nomalization_forImage(protonDensity2,3 );
figure,imagesc(abs(mask_forProton.*protonDensity2),[0,1.8]);colormap('gray');
figure,imagesc(abs(mask_forProton.*S02_EPTI),[0,8]);colormap('gray');

ratio = S02_EPTI./S0_EPTI;
figure,imagesc(abs(ratio),[0,1]);colormap('jet');
%%
fig1=figure; 
subplot(1,2,1),imagesc(abs(T2_EPTI),[0 250]); colormap('hot'); colorbar;h=title('T2');set(h,'FontSize',18);daspect([1 1 1]);removeTicks
subplot(1,2,2),imagesc(abs(T2s_EPTI),[0 200]); colormap('hot'); colorbar;h=title('T2*');set(h,'FontSize',18);daspect([1 1 1]);removeTicks