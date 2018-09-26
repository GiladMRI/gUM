XFMStr='Daubechies';

filterSize=4;
wavScale=4;

XFM = Wavelet(XFMStr,filterSize,wavScale);

W=XFM*double(MLN);
S=sort(abs(W2(:)),'descend');
W2(abs(W2)<S(8000))=0;