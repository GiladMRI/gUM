NNRecSx=load('S01_r1e1_MLN_NN_recon.mat');
NNRecSx=NNRecSx.NNRecSx;
NNRecSx=single(NNRecSx);
MNRecSx=mean(abs(NNRecSx),4);
Msk=MNRecSx>0.5;
Ttls={'Odd','Even'};
%%
% Odd=abs(NNRecSx(:,:,:,1:2:end));
figure;
Ttls={'Odd','Even'};
for i=1:2
%     Data=(NNRecSx(:,:,:,i:2:end));
%     mData=mean(Data,4);
%     sData=std(Data,[],4);
%     tData=abs(mData./sData);
%     subplot(4,2,i*4-3);
%     gmontage(tData);title([Ttls{i} ' without abs']);colorbar
%     [N,c]=hist(tData(Msk),100);
%     subplot(4,2,i*4-1);
%     plot(c,N);
%     
    
    Data=abs(NNRecSx(:,:,:,i:2:end));
    mData=mean(Data,4);
    sData=std(Data,[],4);
    tData=mData./sData;
    subplot(2,2,i*2-1);
    gmontage(tData.*Msk);title([Ttls{i}]);colorbar
    [N,c]=hist(tData(Msk),100);
    subplot(2,2,i*2);
    plot(c,N);
end
%%
Nj=2;
NN=120;
figure;
i=2;
for j=1:Nj
    Ids=(j-1)*NN+10+(1:NN);
    Data=abs(NNRecSx(:,:,:,Ids(i:2:end)));
    mData=mean(Data,4);
    sData=std(Data,[],4);
    tData=mData./sData;
    subplot(2,Nj,j);
    gmontage(tData.*Msk);title([Ttls{i} ', Idxs [' num2str(min(Ids)) ' ' num2str(max(Ids)) ']'])
    [N,c]=hist(tData(Msk),100);
    subplot(2,Nj,j+Nj);
    plot(c,N);
    hold on;
    
    Data=abs(NNRecSx(:,:,:,Ids((3-i):2:end)));
    mData=mean(Data,4);
    sData=std(Data,[],4);
    tData=mData./sData;
    [N,c]=hist(tData(Msk),100);
    plot(c,N,'r');
    setXaxis([0 500]);
end
%%
Data=abs(NNRecSx(:,:,:,2:2:end-1));
mData=mean(Data,4);
sData=std(Data,[],4);
tData=mData./sData;
fgmontage(tData.*Msk);title([Ttls{2}]);colorbar

figure;
[N,c]=hist(tData(Msk),100);
%     subplot(2,2,i*2);
    plot(c,N);
%%
MM=(tData<100).*Msk;
% fgmontage(MM);
TwoD=Reshape4d22d(Data,MM);
figure;plot(mean(TwoD,1));
%%
Odd=abs(NNRecSx(:,:,:,1:2:end));
Even=abs(NNRecSx(:,:,:,2:2:end));
mOdd=mean(Odd,4);
mEven=mean(Even,4);
fgmontage(mEven-mOdd,[0 1e-2])
%%
Perf=abs(NNRecSx(:,:,:,2:2:end))-abs(NNRecSx(:,:,:,1:2:end));
% Perf=abs(NNRecSx(:,:,:,2:2:end)-NNRecSx(:,:,:,1:2:end));
mPerf=mean(Perf,4);
sPerf=std(Perf,[],4);
tPerf=mPerf./sPerf;
fgmontage(tPerf,[0 2])
%%
% F=fft2cg(NNRecSx);
% [X Y]=meshgrid(linspace(-1,1,96),linspace(-1,1,96));
% R=sqrt(X.^2+Y.^2);
% M=R>0.05;
% NNRecSxF=F.*repmat(M,[1 1 24 270]);
% NNRecSxF=ifft2cg(NNRecSxF);
% NNRecSx=NNRecSxF;
