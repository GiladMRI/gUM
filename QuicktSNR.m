XX=grmss(Rec_CompgB0_MX(:,:,:,10:30),4);

sSig=std(XX,0,3);
mSig=mean(XX,3);
tSig=mSig./sSig;

figure;
subplot(1,3,1);gmontage(mSig);
subplot(1,3,2);gmontage(sSig);
subplot(1,3,3);gmontage(tSig);colorbar