% Load M and P to C
M=loadniidata('/autofs/cluster/kawin/Gilad/037_FirstEcho_M.nii');
P=loadniidata('/autofs/cluster/kawin/Gilad/037_FirstEcho_P.nii');
C=M.*exp(1i*P);
clear M P
Cs=imresize3(C,[64 64 64]);
nSlices=size(Cs,3);
%% Poisson masks
Accs=1.01:0.01:1.2;
for i=1:numel(Accs)
    Acc=Accs(i);
    for s=1:nSlices
        M(:,:,i,s)=squeeze(bart(['poisson -Y 64 -Z 64 -y ' num2str(Acc) ' -z ' num2str(Acc) ' -C 4 -s ' num2str(randi(1e13))]));
    end
end
%% Recon with W 2D,3D, different lambdas
Lambdas=10.^(-7:0.5:2);
F=fft2cg(Cs);
%%
for i=5:numel(Accs)
    Data=F.*squeeze(M(:,:,i,:));
    DataP=permute(Data,[1 2 6 5 4 3]);
    
    for l=1:numel(Lambdas)
        Rec2DT(:,:,:,i,l)=squeeze(bart(['pics -m -R T:3:3:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec2DW(:,:,:,i,l)=squeeze(bart(['pics -m -R W:3:0:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec3DT(:,:,:,i,l)=squeeze(bart(['pics -m -R T:35:35:' num2str(Lambdas(l))],DataP,DataP*0+1));
        Rec3DW(:,:,:,i,l)=squeeze(bart(['pics -m -R W:35:0:' num2str(Lambdas(l))],DataP,DataP*0+1));
    end
end
disp('Rec finished');
%%
% save('CompareRecon2D3Ds.mat','Rec2DT','Rec2DW','Rec3DT','Rec3DW');
%% Normalize
Rec2DTN=Rec2DT.*grmss(Cs)./grms(Rec2DT,1:3);
Rec3DTN=Rec3DT.*grmss(Cs)./grms(Rec3DT,1:3);
Rec2DWN=Rec2DW.*grmss(Cs)./grms(Rec2DW,1:3);
Rec3DWN=Rec3DW.*grmss(Cs)./grms(Rec3DW,1:3);
%% SSIM
clear S2DT S3DT S2DW S3DW
for i=1:numel(Accs)
    disp(i);
    for l=1:numel(Lambdas)
        disp([i l]);
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
% figure;plot(MM);legend({'2DT','2DW','3DT','3DW'})
figure;
plot(S2DTmm,'k--','LineWidth',2);hold on;
plot(S2DWmm,'k-','LineWidth',2);
plot(S3DTmm,'r--','LineWidth',2);
plot(S3DWmm,'r-','LineWidth',2);
legend({'2D TV','2D Wavelet','3D TV','3D Wavelet'},'FontSize',20)
xlabel('Acceleration','FontSize',20)
ylabel('SSIM','FontSize',20)
title('SSIM vs. undersampling, poisson disc','FontSize',20);
%%
W2DSciptFN='wav2D.txt';
W3DSciptFN='wav3D.txt';
fid=fopen(W2DSciptFN,'wt');fprintf(fid,"wavelet 3 0\n");fclose(fid);
fid=fopen(W3DSciptFN,'wt');fprintf(fid,"wavelet 7 0\n");fclose(fid);
W2D=bart(['linopScript ' W2DSciptFN],FillOnesTo16(size(Cs)),Cs);
W3D=bart(['linopScript ' W3DSciptFN],FillOnesTo16(size(Cs)),Cs);
%%
SW2D=sort(abs(W2D(:).^2),'descend');
cSW2D=cumsum(SW2D)./sum(SW2D);
SW3D=sort(abs(W3D(:).^2),'descend');
cSW3D=cumsum(SW3D)./sum(SW3D);
figure;plot(linspace(0,1,numel(cSW2D)),cSW2D,'k');hold on
plot(linspace(0,1,numel(cSW3D)),cSW3D,'r');