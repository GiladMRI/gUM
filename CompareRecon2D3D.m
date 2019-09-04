% Load M and P to C
M=loadniidata('/autofs/cluster/kawin/Gilad/037_FirstEcho_M.nii');
P=loadniidata('/autofs/cluster/kawin/Gilad/037_FirstEcho_P.nii');
C=M.*exp(1i*P);
clear M P
Cs=imresize3(C,[64 64 64]);
%% Poisson masks
Accs=1.01:0.01:1.2;
for i=1:numel(Accs)
    Acc=Accs(i);
    M(:,:,i)=squeeze(bart(['poisson -Y 64 -Z 64 -y ' num2str(Acc) ' -z ' num2str(Acc) ' -C 4']));
end
%% Recon with W 2D,3D, different lambdas
Lambdas=10.^(-7:0.5:2);
F=fft2cg(Cs);
%%
for i=1:numel(Accs)
    Data=F.*M(:,:,i);
    DataP=permute(Data,[1 2 6 5 4 3]);
    
    for l=1:numel(Lambdas)
        Rec2DT(:,:,:,i,l)=squeeze(bart(['pics -m -R T:3:3:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec2DW(:,:,:,i,l)=squeeze(bart(['pics -m -R W:3:0:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec3DT(:,:,:,i,l)=squeeze(bart(['pics -m -R T:35:35:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec3DW(:,:,:,i,l)=squeeze(bart(['pics -m -R W:35:0:' num2str(Lambdas(l))],DataP,DataP*0+1));
    end
end
%% Normalize
Rec2DTN=Rec2DT.*grmss(Cs)./grms(Rec2DT,1:3);
Rec3DTN=Rec3DT.*grmss(Cs)./grms(Rec3DT,1:3);
Rec2DWN=Rec2DW.*grmss(Cs)./grms(Rec2DW,1:3);
Rec3DWN=Rec3DW.*grmss(Cs)./grms(Rec3DW,1:3);
%% SSIM
for i=1:numel(Accs)
    for l=1:numel(Lambdas)
        for s=1:size(Cs,3)
            S2DT(s,i,l)=ssim(abs(squeeze(Rec2DTN(:,:,s,i,l))),abs(Cs(:,:,s)));
            S3DT(s,i,l)=ssim(abs(squeeze(Rec3DTN(:,:,s,i,l))),abs(Cs(:,:,s)));
            S2DW(s,i,l)=ssim(abs(squeeze(Rec2DWN(:,:,s,i,l))),abs(Cs(:,:,s)));
            S3DW(s,i,l)=ssim(abs(squeeze(Rec3DWN(:,:,s,i,l))),abs(Cs(:,:,s)));
        end
    end
end
disp('Finished ssim')
%% Compare
S2DTm=squeeze(mean(S2DT,1));
S3DTm=squeeze(mean(S3DT,1));
S2DWm=squeeze(mean(S2DW,1));
S3DWm=squeeze(mean(S3DW,1));
S2DTmm=max(S2DTm,[],2);
S3DTmm=max(S3DTm,[],2);
S2DWmm=max(S2DWm,[],2);
S3DWmm=max(S3DWm,[],2);
MM=[S2DTmm S2DWmm S3DTmm S3DWmm];
figure;plot(MM);legend({'2DT','2DW','3DT','3DW'})















