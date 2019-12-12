BaseP='/autofs/space/daisy_001/users/data/Gilad/gep_CL/';
FN='meas_MID01928_FID43869_gre_4echo_24_26_G2';
% FN='meas_MID01929_FID43870_gre_4echo_37_26_G2';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX=sTwix{end};
% Data=sTwixX.image();

RefData=sTwixX.refscan();
%%
RefDataP=permute(RefData(:,:,:,:,:,:,:,1),[1 3 2 5 4]);
iRef=ifft2cg(RefDataP); % [RO PE Ch Sli]
fgmontage(grmss(iRef,3));daspect([1 8 1])

RefDataPME=permute(RefData(:,:,:,:,:,:,:,:),[1 3 2 5 8 4 6 7]);

%%
DATA=RefDataP(:,:,:,5);
[sx,sy,Nc] = size(DATA);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size

% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);
%% ggg quick
[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
%%
SliI=5;
nEchosToUse=4;
calibME = crop(squeeze(RefDataPME(:,:,:,SliI,:)),ncalib,ncalib,Nc,nEchosToUse);
% dat2Kernel
[sx,sy,nc,ne] = size(calibME);
imSize = [sx,sy] ;

for i=1:nEchosToUse
    tmpC{i} = im2row(calibME(:,:,:,i),ksize);
end
tmpM=cat(1,tmpC{:});
[tsx,tsy,tsz] = size(tmpM);
for i=1:nEchosToUse
    AC{i} = reshape(tmpM,tsx,tsy*tsz);
end
A=cat(1,AC{:});

[U,S,V] = svd(A,'econ');
    
kernel = reshape(V,ksize(1),ksize(2),nc,size(V,2));
S = diag(S);S = S(:);
% dat2Kernel end

idx = max(find(S >= S(1)*eigThresh_1));
[M,W] = kernelEig(kernel(:,:,:,1:idx),[sx,sy]);
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
%%
SliI=5;
nEchosToUse=4;
calibME = crop(squeeze(RefDataPME(:,:,:,SliI,:)),ncalib,ncalib,Nc,nEchosToUse);

[M,W] = ME_ESPIRIT(calibME,[6,6],eigThresh_1);
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
%%
FCombined=fft2cg(Combined); % [x y Slices Echos]
Sz=gsize(FCombined,1:2);
% CalibRegion=[30 30];
CalibRegion=Sz;
SliI=12;
HankelTemporalLen=2;
CurFCombined=squeeze(FCombined(:,:,SliI,:));
[~, ~, ~,H]=ghankel(size(CurFCombined,3),HankelTemporalLen,gsize(CurFCombined,1:2));
HCurSli=H*squeeze(Combined(:,:,SliI,:));
% HCurFCombined=H*CurFCombined;
HCurFCombined=fft2cg(HCurSli);
HCurFCombinedP=perm43(HCurFCombined);
% kSize=[6,6];
% kSize=[9,9];
% kSize=[12,12];
% kSize=[16,16];
kSize=[24,24];
calibME = crop(HCurFCombinedP,CalibRegion(1),CalibRegion(2),HankelTemporalLen,size(CurFCombined,3)-HankelTemporalLen+1);

[M,W] = ME_ESPIRIT(calibME,kSize,eigThresh_1);

maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,HankelTemporalLen]);
mapsB=imresize(maps,Sz);
tmp=mapsB(:,:,2)./mapsB(:,:,1);
% fgmontagex(angle(tmp))

M1=squeeze(M(:,:,:,2)./M(:,:,:,1));
M2=squeeze(M(:,:,2,:)./M(:,:,1,:));

M1a=squeeze(M(:,:,:,2).*M(:,:,:,1));
M2a=squeeze(M(:,:,2,:).*M(:,:,1,:));

tmp1=imresize(M2(:,:,2),Sz);
tmp=gflip(angle(conj(tmp1)),[]);
% fgmontagex(tmp)
fgmontagex(tmp,[-pi pi]);title(num2str([CalibRegion kSize],' %d'));
%%
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(perm43(Combined),WhichEchosToUse,dTEs_ms(1),TEs_ms(1));
%%
% Display coil images: 
im = ifft2c(DATA);
%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

%% 
% This shows that the calibration matrix has a null space as shown in the
% paper. 
kdisp = reshape(k,[ksize(1)*ksize(2)*Nc,ksize(1)*ksize(2)*Nc]);
figure, subplot(211), plot([1:ksize(1)*ksize(2)*Nc],S,'LineWidth',2);
hold on, 
plot([1:ksize(1)*ksize(2)*Nc],S(1)*eigThresh_1,'r-','LineWidth',2);
plot([idx,idx],[0,S(1)],'g--','LineWidth',2)
legend('signular vector value','threshold')
title('Singular Vectors')
subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
xlabel('Singular value #');
title('Singular vectors')
%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);
%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps
figure, imshow3(abs(W),[],[4,8]); 
title('Eigen Values in Image space');
colormap((gray(256))); colorbar;

figure, imshow3(abs(M),[],[32,32]); 
title('Magnitude of Eigen Vectors');
colormap(gray(256)); colorbar;

figure, imshow3(angle(M),[],[32,32]); 
title('Magnitude of Eigen Vectors');
colormap(jet(256)); colorbar;
%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 


% alternate way to compute projection is:
% ESP = ESPIRiT(M);
% P = ESP'*im;

P = sum(repmat(im,[1,1,1,Nc]).*conj(M),3);
figure, imshow3(abs(P),[],[4,8]); 
title('Magnitude of Eigen Vectors');
colormap(sqrt(gray(256))); colorbar;

figure, imshow3(angle(P),[],[4,8]); 
title('Magnitude of Eigen Vectors');
colormap((jet(256))); colorbar;
%%
% crop sensitivity maps 
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);

figure, imshow3(abs(maps),[],[4,8]); 
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps),[],[4,8]); 
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;
