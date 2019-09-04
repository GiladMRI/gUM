Sz=[128 128];
I1=imresize(cameraman,Sz);
I1=I1./gmax(I1);
I1(1)=I1(1)+1i*0.00001;
I2=double(234-imresize(rgb2gray(imread('hands1.jpg')),Sz));
I2=I2./gmax(I2);
I2(1)=I2(1)+1i*0.00001;
N1=128;
N2=128;
v=I1;

Sz128=Sz;
Trajm2=
%% Test with same, complex, wi
% wi=ones(nTraj,1);
% wi=sin(linspace(0,3*pi,nTraj)).';
CAIPISep_mm=36;
CAIPIDelay_us=200;
CAIPIPeriod_us=200;
CAIPIVec=CAIPIBlips([paramLongSpGradAmp, paramLongSpSlewRate,CAIPISep_mm,CAIPIDelay_us,CAIPIPeriod_us,...
    2560*paramLongROSamples/1024]);
gamma=42.5774806;
cCAIPIVec=cumsum(CAIPIVec)*gamma*10*2*pi/1e6;
cCAIPIVecX=interp1(1:numel(cCAIPIVec),cCAIPIVec,1:1e5/spBW:(numel(cCAIPIVec)-0.01))*CAIPISep_mm;
wi=exp(1i*cCAIPIVecX);

% wi=sin(linspace(0,3*pi,nTraj)).'+1i*sin(linspace(0,4*pi,nTraj)).';
wiforblock1=wi.*conj(wi);
wiforblock2=wi.*conj(wi);
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om

block1=nufft_adj(wiforblock1, nufftStruct);
block2=nufft_adj(wiforblock2, nufftStruct2);

fftkern=blocksToFFTKern(block1,block2);

fftkern=single(real(fftkern));

Sig=nufft(I1,nufftStruct);
Sig=Sig.*wi;
Sig=Sig.*conj(wi);
GGv=nufft_adj(Sig, nufftStruct);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N2);

rDiff=abs(GGv(:)-CFHFPadv(:))./abs(GGv(:));
rDiff2=abs(GGv(:)-CFHFPadv(:))/mean(abs(GGv(:)));
[max(rDiff) mean(rDiff) max(rDiff2) mean(rDiff2)]*100
%%
wi=ones(nTraj,1);

block1=nufft_adj(wi, nufftStruct);
block2=nufft_adj(wi, nufftStruct2);

% fftkern=blocksToFFTKern(block1,block2);
N1=size(block1,1);
N2=size(block1,2);
z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
    [block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
    zeros(1,2*N2);
    [flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
    ]; % [(2Nd)]
n1 = 0:(N1*2-1);
n2 = 0:(N2*2-1);
[n1, n2] = ndgrid(n1, n2);
tmp1 = kern(sub2ind([N1*2 N2*2], 1+n1, 1+n2));
m1 = mod(-n1, N1*2); % circular symmetry
m2 = mod(-n2, N2*2);
tmp2 = kern(sub2ind([N1*2 N2*2], 1+m1, 1+m2));
tmp2 = conj(tmp2);
kern = (tmp1+tmp2)/2; % force it to be Hermitian
fftkern = fftn(kern);
%% Test with 2 complex wis
% wi=ones(nTraj,1);
w1=sin(linspace(0,3*pi,nTraj)).'+1i*sin(pi/4+linspace(0,pi,nTraj)).';
w2=sin(linspace(0,pi,nTraj)).'+1i*sin(linspace(0,2*pi,nTraj)).';

wx=w1.*conj(w2);
swx=sqrt(wx);
wiforblock1=w1.*conj(w2);
wiforblock2=w1.*conj(w2);
% wiforblock1=w1.*conj(w2);
% wiforblock2=w1.*(w2);
% wiforblock1=w2.*conj(w2);
% wiforblock2=w1.*conj(w1);
%%
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
% nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
fftkernBigc=NUFFT_to_Toep(nufftStruct,wx);


% X=Bigc;figure;subplot(1,2,1);imshow(abs(X),[0 400]);subplot(1,2,2);imshow(angle(X),[-pi pi]);title(Bigc)
% % X=Bigc';figure;subplot(1,2,1);imshow(abs(X),[0 400]);subplot(1,2,2);imshow(angle(X),[-pi pi]);title(Bigc)
% X=circshift(gflip(Bigc,1:2),[1 1]);figure;subplot(1,2,1);imshow(abs(X),[0 400]);subplot(1,2,2);imshow(angle(X),[-pi pi]);title(Bigc)

%%
Sig=nufft(I1,nufftStruct);
Sig=Sig.*wx;
AHAI1=nufft_adj(Sig, nufftStruct);
X=AHAI1;figure;subplot(1,2,1);imshow(abs(X),[]);subplot(1,2,2);imshow(angle(X),[-pi pi]);title(AHAI1)
subplot(1,2,1);title('AHAI1');colorbar
%%

Padv=I1;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkernBigc.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadvBigc=FHFPadv(1:N1,1:N2);
% CFHFPadv=FHFPadv(N1+1:end,N2+1:end);

X=CFHFPadvBigc;figure;subplot(1,2,1);imshow(abs(X),[]);subplot(1,2,2);imshow(angle(X),[-pi pi]);
subplot(1,2,1);title('CFHFPadvBigc');colorbar
%%
X=AHAI1-CFHFPadvBigc;figure;subplot(1,2,1);imshow(abs(X),[]);subplot(1,2,2);imshow(angle(X),[-pi pi]);
subplot(1,2,1);title('diff');colorbar
grmss(AHAI1)/grmss(AHAI1-CFHFPadvBigc)
%% now SMS
w1=sin(linspace(0,3*pi,nTraj)).'+1i*sin(pi/4+linspace(0,pi,nTraj)).';
w2=sin(linspace(0,pi,nTraj)).'+1i*sin(linspace(0,2*pi,nTraj)).';
w1=w1*0+1;
w2=w2*0+1;

w12=w1.*conj(w2);
w21=w2.*conj(w1);

fftkern11=NUFFT_to_Toep(nufftStruct,w1.*conj(w1));
fftkern22=NUFFT_to_Toep(nufftStruct,w2.*conj(w2));
fftkern12=NUFFT_to_Toep(nufftStruct,w1.*conj(w2));
fftkern21=NUFFT_to_Toep(nufftStruct,w2.*conj(w1));
%%
Sig=nufft(I1,nufftStruct);
Sig1=Sig.*w1;
Sig=nufft(I2,nufftStruct);
Sig2=Sig.*w2;
SigS=Sig1+Sig2;
AHAI1=nufft_adj(SigS.*conj(w1), nufftStruct);
AHAI2=nufft_adj(SigS.*conj(w2), nufftStruct);
AHAI12=cat(3,AHAI1,AHAI2);
X=AHAI12;figure;subplot(1,2,1);gmontage(abs(X),[0 900000]);subplot(1,2,2);gmontage(angle(X),[-pi pi]);title(AHAI1)
subplot(1,2,1);title('AHAI12');colorbar
%%
Padv=I1;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern11.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv11=FHFPadv(1:N1,1:N2);

Padv=I2;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern22.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv22=FHFPadv(1:N1,1:N2);

Padv=I1;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern12.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv12=FHFPadv(1:N1,1:N2);

Padv=I2;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern21.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv21=FHFPadv(1:N1,1:N2);

CFHFPadv1=CFHFPadv11+CFHFPadv21;
CFHFPadv2=CFHFPadv12+CFHFPadv22;

TBoth=cat(3,CFHFPadv1,CFHFPadv2);

X=TBoth;figure;subplot(1,2,1);gmontage(abs(X),[0 900000]);subplot(1,2,2);gmontage(angle(X),[-pi pi]);
subplot(1,2,1);title('TBoth');colorbar

grmss(AHAI12)/grmss(AHAI12-TBoth)