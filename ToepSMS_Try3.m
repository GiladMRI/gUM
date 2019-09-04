A=load('/autofs/space/daisy_002/users/Gilad/All_Orientation-0x.mat');
load('/autofs/space/daisy_002/users/Gilad/CCSensMaps.mat');
%%
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
% Trajm2=BARTTrajx;
%% Test with same, complex, wi
% wi=ones(nTraj,1);
% wi=sin(linspace(0,3*pi,nTraj)).';
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;
paramLongROSamples=12288;
spBW=400000;
CAIPISep_mm=36;
CAIPIDelay_us=200;
CAIPIPeriod_us=200;
load('Trajm2.mat');
Sz128=[128 128];
CAIPIVec=CAIPIBlips([paramLongSpGradAmp, paramLongSpSlewRate,CAIPISep_mm,CAIPIDelay_us,CAIPIPeriod_us,...
    2560*paramLongROSamples/1024]);
gamma=42.5774806;
cCAIPIVec=cumsum(CAIPIVec)*gamma*10*2*pi/1e6;
cCAIPIVecX=interp1(1:numel(cCAIPIVec),cCAIPIVec,1:1e5/spBW:(numel(cCAIPIVec)-0.01))*CAIPISep_mm;
wi=exp(1i*cCAIPIVecX);

% wi=sin(linspace(0,3*pi,nTraj)).'+1i*sin(linspace(0,4*pi,nTraj)).';
wiforblock1=(wi.*conj(wi)).';
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om

fftkern=NUFFT_to_Toep(nufftStruct,wiforblock1);
fftkern=single(real(fftkern));

SN=nufftStruct.sn;
try
    P=nufftStruct.p.G;
catch
    P=nufftStruct.p;
end

Fac=prod(Sz128)*4;

% TakeTopLeftBlock=@(x,Sz) x(1:Sz(1),1:Sz(2));
TakeTopLeftBlock=@(x,Sz) x(1:Sz(1),1:Sz(2),:,:,:,:);
NUFT_forw=@(x) P*reshape(fft2(padarray(x.*SN,Sz128,'post')),[],1)*Fac;
NUFT_adj=@(x) TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128*2)),Sz128).*conj(SN);
% NUFT_adj=@(x) subsref(ifft2(reshape(P'*x,Sz128*2)),struct('type','()','subs',{{1:Sz128(1),1:Sz128(2)}})).*conj(SN);

NUFT_adj_forw=@(x,wa,wb) NUFT_adj((NUFT_forw(x).*wa).*wb);
GGv=NUFT_adj_forw(I1,wi.',conj(wi.'));

NUFT_T=@(x,fkern) TakeTopLeftBlock(ifft2(fkern.*fft2(padarray(x,Sz128,'post'))),Sz128);

Tv=NUFT_T(I1,fftkern);

rDiff=abs(GGv(:)-Tv(:))./abs(GGv(:));
rDiff2=abs(GGv(:)-Tv(:))/mean(abs(GGv(:)));
[max(rDiff) mean(rDiff) max(rDiff2) mean(rDiff2)]*100
%% now SMS
w1=exp(1i*cCAIPIVecX*0).';
w2=exp(1i*cCAIPIVecX).';

fftkern11=NUFFT_to_Toep(nufftStruct,w1.*conj(w1));
fftkern22=NUFFT_to_Toep(nufftStruct,w2.*conj(w2));
fftkern12=NUFFT_to_Toep(nufftStruct,w1.*conj(w2));
fftkern21=NUFFT_to_Toep(nufftStruct,w2.*conj(w1));
%% MB
Sig=NUFT_forw(I1);
Sig1=Sig.*w1;
Sig=NUFT_forw(I2);
Sig2=Sig.*w2;
SigS=Sig1+Sig2;save
AHAI1=NUFT_adj(SigS.*conj(w1));
AHAI2=NUFT_adj(SigS.*conj(w2));
AHAI12=cat(3,AHAI1,AHAI2);
X=AHAI12;figure;subplot(1,2,1);gmontage(abs(X),[0 900000]);subplot(1,2,2);gmontage(angle(X),[-pi pi]);title(AHAI1)
subplot(1,2,1);title('AHAI12');colorbar
%% Toep MB
Tv11=NUFT_T(I1,fftkern11);
Tv12=NUFT_T(I1,fftkern12);
Tv21=NUFT_T(I2,fftkern21);
Tv22=NUFT_T(I2,fftkern22);

Tv1=Tv11+Tv21;
Tv2=Tv12+Tv22;

TBoth=cat(3,Tv1,Tv2);save

X=TBoth;figure;subplot(1,2,1);gmontage(abs(X),[0 900000]);subplot(1,2,2);gmontage(angle(X),[-pi pi]);
subplot(1,2,1);title('TBoth');colorbar

grmss(AHAI12)/grmss(AHAI12-TBoth)
%% With sensitivity maps, no MB
% Sens [X Y channels]
wi=ones(nTraj,1);

% NUFT_forwS=@(x,Sens) P*reshape(fft2(padarray((x.*SN).*Sens,Sz128,'post')),prod(Sz128)*4,[])*Fac;
NUFT_forwS=@(x,Sens) P*reshape(fft2((x.*SN).*Sens,Sz128(1)*2,Sz128(2)*2),prod(Sz128)*4,[])*Fac;
NUFT_adjS=@(x,Sens) sum((TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128(1)*2,Sz128(2)*2,[])),Sz128).*conj(SN)).*conj(Sens),3);

NUFT_adj_forwS=@(x,Sens,wa,wb) NUFT_adjS((NUFT_forwS(x,Sens).*wa).*wb,Sens);

% NUFT_TS=@(x,Sens,fkern) sum(TakeTopLeftBlock(ifft2(fkern.*fft2(padarray(x.*Sens,Sz128,'post'))),Sz128).*conj(Sens),3);
NUFT_TS=@(x,Sens,fkern) sum(TakeTopLeftBlock(ifft2(fkern.*fft2(x.*Sens,Sz128(1)*2,Sz128(2)*2)),Sz128).*conj(Sens),3);

GGvS=NUFT_adj_forwS(I1,Sens,wi,conj(wi));
TvS=NUFT_TS(I1,Sens,fftkern);

rDiffS=abs(GGvS(:)-TvS(:))./abs(GGvS(:));
rDiff2S=abs(GGvS(:)-TvS(:))/mean(abs(GGvS(:)));
[max(rDiffS) mean(rDiffS) max(rDiff2S) mean(rDiff2S)]*100
%% With sensitivity maps, MB
SensMB=repmat(Sens,[1 1 1 2]);
SensMB=double(permute(SensCC(1:8,:,:,[10 20]),[2 3 1 4]));

IMB=double(cat(4,I1,I2));
wC=cat(3,w1,w2);

fftkerns(:,:,1,1,1)=fftkern11;
fftkerns(:,:,1,1,2)=fftkern12;
fftkerns(:,:,1,2,1)=fftkern21;
fftkerns(:,:,1,2,2)=fftkern22;

% NUFT_forwSMB=@(x,SensMB) reshape(P*reshape(fft2(padarray((x.*SN).*SensMB,Sz128,'post')),prod(Sz128)*4,[])*Fac,[],size(SensMB,3),size(SensMB,4));
NUFT_forwSMB=@(x,SensMB) reshape(P*reshape(fft2((x.*SN).*SensMB,Sz128(1)*2,Sz128(2)*2),prod(Sz128)*4,[])*Fac,[],size(SensMB,3),size(SensMB,4));
NUFT_adjSMB=@(x,SensMB) sum((TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128(1)*2,Sz128(2)*2,size(SensMB,3),size(SensMB,4))),Sz128).*conj(SN)).*conj(SensMB),3);

NUFT_forwSMBC=@(x,SensMB,wa) sum(NUFT_forwSMB(x,SensMB).*wa,3);
NUFT_adjSMBC=@(x,SensMB,wb) NUFT_adjSMB(reshape(x.*wb,[],size(SensMB,3)*size(SensMB,4)),SensMB);

NUFT_adj_forwSMBC=@(x,SensMB,wa,wb) NUFT_adjSMBC(NUFT_forwSMBC(x,SensMB,wa),SensMB,wb);

% NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(padarray(x.*SensMB,Sz128,'post'))),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);
% NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(x.*SensMB,Sz128(1)*2,Sz128(2)*2)),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);
% NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(x.*SensMB,Sz128(1)*2,Sz128(2)*2)),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);

NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(x.*SensMB,Sz128(1)*2,Sz128(2)*2)),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);

permute54=@(x) permute(x,[1 2 3 5 4]);
TakeTopLeftBlocka=@(x) x(1:Sz128(1),1:Sz128(2),:,:,:,:);
Padded2DFT=@(x) fft2(x,Sz128(1)*2,Sz128(2)*2);
Truncated2DIFT=@(x) TakeTopLeftBlocka(ifft2(x));
NUFT_TSMBC=@(x,SensMB,fkerns) sum(Truncated2DIFT(permute54(sum(fkerns.*Padded2DFT(x.*SensMB),4))).*conj(SensMB),3);


GGvSMBC=NUFT_adj_forwSMBC(IMB,SensMB,wC,conj(wC));
TvSMBC=NUFT_TSMBC(IMB,SensMB,fftkerns);

FF=fft2(IMB.*SensMB,Sz128(1)*2,Sz128(2)*2);
FFH=fftkerns.*FF;
IFFH=TakeTopLeftBlock(ifft2(FFH),Sz128);
% SnIFFH=IFFH.*conj(SensMB);
SnIFFH=IFFH.*conj(permute(SensMB,[1 2 3 5 4]));
SSnIFFH=sum(SnIFFH,3);
SSSnIFFH=permute(sum(SSnIFFH,4),[1 2 3 5 4]);

mean(abs(GGvSMBC(:)-SSSnIFFH(:))/mean(abs(GGvSMBC(:))))*100

FF=fft2(IMB.*SensMB,Sz128(1)*2,Sz128(2)*2);
FFH=fftkerns.*FF;
IFFH=TakeTopLeftBlock(ifft2(FFH),Sz128);
SIFFH=permute(sum(IFFH,4),[1 2 3 5 4]);
SnSIFFH=sum(SIFFH.*conj(SensMB),3);

mean(abs(GGvSMBC(:)-SnSIFFH(:))/mean(abs(GGvSMBC(:))))*100

FF=fft2(IMB.*SensMB,Sz128(1)*2,Sz128(2)*2);
FFH=fftkerns.*FF;
SFFH=permute(sum(FFH,4),[1 2 3 5 4]);
ISFFH=TakeTopLeftBlock(ifft2(SFFH),Sz128);
SnISFFH=sum(ISFFH.*conj(SensMB),3);

mean(abs(GGvSMBC(:)-SnISFFH(:))/mean(abs(GGvSMBC(:))))*100

fftkerns3=fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,1)=fftkerns(:,:,1,1,2)-fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,2)=fftkerns(:,:,1,2,1)-fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,3)=fftkerns(:,:,1,1,1);
tic
FF=fft2(IMB.*SensMB,Sz128(1)*2,Sz128(2)*2);
EFF=cat(4,FF,sum(FF,4));
HEFF=EFF.*fftkerns3;
% XHEFF=cat(4,HEFF(:,:,:,2)+HEFF(:,:,:,3),HEFF(:,:,:,1)+HEFF(:,:,:,3));
XHEFF=HEFF(:,:,:,[2 1])+HEFF(:,:,:,3);
IXHEFF=Truncated2DIFT(XHEFF);
SnIXHEFF=sum(IXHEFF.*conj(SensMB),3);
toc

mean(abs(GGvSMBC(:)-SnIXHEFF(:))/mean(abs(GGvSMBC(:))))*100

tic
TvSMBC=NUFT_TSMBC(IMB,SensMB,fftkerns);
toc

rDiffSMBC=abs(GGvSMBC(:)-TvSMBC(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBC(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDiffSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100
%% With trick
%% ADMM
iter_ops.max_iter = 13;
iter_ops.rho = .1;

llr_ops.lambda = 0.400;tic

lsqr_ops.max_iter = 10;
lsqr_ops.tol = 1e-4;

alpha_ref = RefValue;
alpha_ref.data = zeros(Sz128);

% AHA_I=NUFT_adj(SigS);
AHA_I=NUFT_TS(I1,Sens,fftkernClean);
% AHA=@(x) NUFT_adj_forwS(x,Sens,1,1)
AHA=@(x) NUFT_TS(x,Sens,fftkernClean);

% proxf=@(x,lambda) TV_thresh(x,lambda);
proxf=@(x,lambda) TV_thresh_rotationInvariant(x,lambda);
history = iter_admm_1prox(alpha_ref, iter_ops, llr_ops, lsqr_ops, AHA, AHA_I, proxf, @admm_callback);

admmRes=alpha_ref.data;

ShowAbsAngle(admmRes)
%% Speed test
nSpeesTest=10;
% MB=5;
for MB=1:16
    disp(MB);
    clear fftkernsX
    for i=1:MB
        for j=1:MB
            wi=randn(nTraj,1)+1i*randn(nTraj,1);
            wj=randn(nTraj,1)+1i*randn(nTraj,1);
            fftkernsX(:,:,1,i,j)=NUFFT_to_Toep(nufftStruct,wi.*conj(wj));
        end
    end
    
    SensMBX=repmat(Sens,[1 1 1 MB]);
    IMBX=repmat(I1,[1 1 1 MB]);
    wCX=repmat(w1,[1 1 MB]);
    
    % NUFT_adj_forwSMBC=@(x,SensMB,wa,wb) NUFT_adjSMBC(NUFT_forwSMBC(x,SensMB,wa),SensMB,wb);
    
    % NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(padarray(x.*SensMB,Sz128,'post'))),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);
    
    tic
    for i=1:nSpeesTest
        GGvSMBC_SpeedTest=NUFT_adj_forwSMBC(IMBX,SensMBX,wCX,conj(wCX));
    end
    tFB(MB)=toc;
    tic
    for i=1:nSpeesTest
        TvSMBC_SpeedTest=NUFT_TSMBC(IMBX,SensMBX,fftkernsX);
    end
    tToep(MB)=toc;
end
%%
figure;plot(tFB,'k','LineWidth',2);hold on;plot(tToep,'b','LineWidth',2);
legend({'Forw-adjoint','Toeplitz'});
xlabel(MB);
ylabel('Seconds');
title(['Time for ' num2str(nSpeesTest) ' runs, 8 channels, CPU']);
%%
figure;plot(tFB./tToep,'k','LineWidth',2);hold on; plot([1 16],[1 1],'r');
xlabel(MB);
ylabel('Seconds');
setXaxis([0 17]);
title(['Calculation time Forward-Adjoint/Toeplitz ratio, ' num2str(nSpeesTest) ' runs, 8 channels, CPU']);
%% ADMM MB
SensMB=permute(SensCC(:,:,:,[10 20]),[2 3 1 4]);
IMB=imresizeBySlices(permute(A.CurSetAll(1,:,:,[100 104]).*exp(1i*A.CurSetAll(2,:,:,[100 104])),[2 3 1 4]),Sz128);

iter_ops.max_iter = 13;
iter_ops.rho = .1;

llr_ops.lambda = 0.400;

lsqr_ops.max_iter = 10;
lsqr_ops.tol = 1e-4;

alpha_ref = RefValue;
alpha_ref.data = zeros([Sz128 1 2]);

% AHA_I=NUFT_adj(SigS);
AHA_I=NUFT_TSMBC(IMB,SensMB,fftkerns);
% AHA=@(x) NUFT_adj_forwS(x,Sens,1,1)
AHA=@(x) NUFT_TSMBC(x,SensMB,fftkerns);

% proxf=@(x,lambda) TV_thresh(x,lambda);
proxf=@(x,lambda) TV_thresh_rotationInvariant16D(x,lambda);
history = iter_admm_1prox(alpha_ref, iter_ops, llr_ops, lsqr_ops, AHA, AHA_I, proxf, @admm_callback);

admmRes=alpha_ref.data;

ShowAbsAngle(admmRes)
%% Now with time segments
% Dims: [H,W,Channels,MB,MBaux,TS,batch]
IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,[100 104]).*exp(1i*A.CurSetAll(2,:,:,[100 104])),[2 3 1 4]),Sz128));
T2Map_ms=imresizeBySlices(double(permute(A.CurSetAll(3,:,:,[100 104]),[2 3 1 4])),Sz128);
SensMB=permute(SensCC(:,:,:,[10 20]),[2 3 1 4]);

% SensMB=repmat(Sens,[1 1 1 2]);
wi=exp(1i*cCAIPIVecX);
wC=cat(3,wi*0+1,wi);
% TSC=
% TSB=

fftkerns(:,:,1,1,1)=fftkern11;
fftkerns(:,:,1,1,2)=fftkern12;
fftkerns(:,:,1,2,1)=fftkern21;
fftkerns(:,:,1,2,2)=fftkern22;

NUFT_forwSMB=@(x,SensMB) reshape(P*reshape(fft2(padarray((x.*SN).*SensMB,Sz128,'post')),prod(Sz128)*4,[])*Fac,[],size(SensMB,3),size(SensMB,4));
NUFT_adjSMB=@(x,SensMB) sum((TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128(1)*2,Sz128(2)*2,size(MB,3),size(SensMB,4))),Sz128).*conj(SN)).*conj(SensMB),3);

NUFT_forwSMBC=@(x,SensMB,wa) sum(NUFT_forwSMB(x,SensMB).*wa,3);
NUFT_adjSMBC=@(x,SensMB,wb) NUFT_adjSMB(reshape(x.*wb,[],size(SensMB,3)*size(SensMB,4)),SensMB);

NUFT_adj_forwSMBC=@(x,SensMB,wa,wb) NUFT_adjSMBC(NUFT_forwSMBC(x,SensMB,wa),SensMB,wb);

NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(padarray(x.*SensMB,Sz128,'post'))),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);

GGvSMBC=NUFT_adj_forwSMBC(IMB,SensMB,wC,conj(wC));
TvSMBC=NUFT_TSMBC(IMB,SensMB,fftkerns);

rDiffSMBC=abs(GGvSMBC(:)-TvSMBC(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBC(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDifsavefSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100


%% forw adj Single band
% NUFT_forw=@(x) P*reshape(fft2(padarray(x.*SN,Sz128,'post')),[],1)*Fac;
% NUFT_adj=@(x) TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128*2)),Sz128).*conj(SN);
RI1=double(IMB(:,:,:,1));
RI1SN=RI1.*SN;
Sig1=NUFT_forw(RI1);
PSig=reshape(P'*Sig1,Sz128*2);
FPSig=ifft2(PSig);
CFPSig=TakeTopLeftBlock(FPSig,Sz128);
SNCFPSig=conj(SN).*CFPSig;
fgmontage(fftshift(log(abs(PSig))));removeTicks; axis equal
fgmontage(FPSig);removeTicks; axis equal
fgmontage(circshift(FPSig,[64 64]));removeTicks;axis equal
fgmontage(CFPSig);removeTicks; axis equal
fgmontage(SNCFPSig);removeTicks; axis equal

figure;TT=linspace(0,10,1000);plot(exp(-TT)+randn(1,numel(TT))/10,'LineWidth',2);removeTicks;
%% Toep single band
RI1=double(IMB(:,:,:,1));
PRI1=padarray(RI1,Sz128,'post');
FPRI1=fft2(PRI1);
HFPRI1=fftkernClean.*FPRI1;
FHFPRI1=ifft2(HFPRI1);
CFHFPRI1=TakeTopLeftBlock(FHFPRI1,Sz128);
fgmontage(circshift(PRI1,[64 64]));removeTicks;axis equal
fgmontage(fftshift(log(abs(FPRI1))));removeTicks; axis equal
fgmontage(fftshift(log(abs(HFPRI1))));removeTicks; axis equal
fgmontage(circshift(FHFPRI1,[64 64]));removeTicks;axis equal
fgmontage(CFHFPRI1);removeTicks; axis equal

fgmontage(fftshift(log(abs(fftkernClean))));removeTicks; axis equal
%% forw adj Multi band
IMBWithS=double(IMB).*SensMB(:,:,2,:);
IMBWithSSN=IMBWithS.*SN;
PIMBWithSSN=padarray(IMBWithSSN,Sz128,'post');
FPIMBWithSSN=fft2(PIMBWithSSN);
Sig11=P*(double(reshape(FPIMBWithSSN(:,:,:,1),[],1)));
Sig12=P*(double(reshape(FPIMBWithSSN(:,:,:,2),[],1)));
Sig11C=Sig11.*wC(:,:,1);
Sig12C=Sig12.*wC(:,:,2);
SumSig=Sig11C+Sig12C;
CSig11C=SumSig.*conj(wC(:,:,1));
CSig12C=SumSig.*conj(wC(:,:,2));

PCSig11C=reshape(P'*CSig11C,Sz128*2);
PCSig12C=reshape(P'*CSig12C,Sz128*2);
PCSigC=cat(4,PCSig11C,PCSig12C);
FPCSigC=ifft2(PCSigC);
CFPCSigC=TakeTopLeftBlock(FPCSigC,Sz128);
SNCFPCSigC=SN.*CFPCSigC;
SSNCFPCSigC=double(SNCFPCSigC).*conj(SensMB(:,:,2,:));
fgmontage(IMB);removeTicks;axis equal
fgmontage(IMBWithS);removeTicks;axis equal
fgmontage(SensMB(:,:,2,:));removeTicks;axis equal
fgmontage(IMBWithSSN);removeTicks;axis equal
fgmontage(circshift(PIMBWithSSN,[64 64]));removeTicks;axis equal
fgmontage(fftshift(log(abs(FPIMBWithSSN))));removeTicks; axis equal
fgmontage(fftshift(log(abs(PCSigC))));removeTicks; axis equal
fgmontage(circshift(FPCSigC,[64 64]));removeTicks;axis equal
fgmontage(CFPCSigC);removeTicks;axis equal
fgmontage(SNCFPCSigC);removeTicks;axis equal
fgmontage(SSNCFPCSigC);removeTicks;axis equal


figure;TT=linspace(0,7,1000);plot(exp(-TT)+randn(1,numel(TT))/10,'b','LineWidth',2);removeTicks;
figure;TT=linspace(0,5,1000);plot(exp(-TT)+randn(1,numel(TT))/10,'r','LineWidth',2);removeTicks;
figure;plot(cCAIPIVecX(1:800),'Color',[0.5 0 0],'LineWidth',4);setYaxis([-0.3 4])
figure;plot(cCAIPIVecX(1:800)*0,'Color',[0 0 0.5],'LineWidth',6);setYaxis([-0.3 4])
figure;TT=linspace(0,5,1000);plot(exp(-TT)+randn(1,numel(TT))/10,'r','LineWidth',2);removeTicks;

figure;plot(abs(SumSig(1:1000)),'k','LineWidth',4);setYaxis([-150000 2000000])
%% Toep Multi band
% NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(padarray(x.*SensMB,Sz128,'post'))),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);

PIMBWithS=padarray(IMBWithS,Sz128,'post');
FPIMBWithS=fft2(PIMBWithS);
HFPIMBWithS=FPIMBWithS.*fftkerns;
FHFPIMBWithS=ifft2(HFPIMBWithS);
CFHFPIMBWithS=TakeTopLeftBlock(FHFPIMBWithS,Sz128);
SCFHFPIMBWithS=CFHFPIMBWithS.*conj(SensMB(:,:,2,:));
SmSCFHFPIMBWithS=sum(SCFHFPIMBWithS,5);
fgmontage(circshift(PIMBWithS,[64 64]));removeTicks;axis equal
fgmontage(fftshift(log(abs(FPIMBWithS))));removeTicks; axis equal
fgmontage(fftshift(log(abs(HFPIMBWithS))));removeTicks; axis equal
fgmontage(circshift(FHFPIMBWithS,[64 64]));removeTicks;axis equal
fgmontage(CFHFPIMBWithS);removeTicks;axis equal
fgmontage(SCFHFPIMBWithS);removeTicks;axis equal
fgmontage(SmSCFHFPIMBWithS);removeTicks;axis equal


fgmontage(fftshift(log(abs(fftkerns))));removeTicks; axis equal
%%
QQ1=fftshift(log(abs(ifft2(fftkern11))));
QQ2=circshift(gflip(fftshift(log(abs(ifft2(fftkern11)))),1:2),[1 1]);
%%
QQ1b=fftshift(log(abs(ifft2(fftkern12))));
QQ2b=circshift(gflip(fftshift(log(abs(ifft2(fftkern12)))),1:2),[1 1]);

%%
I2=IMB(:,:,:,2);
fgmontage(I2);removeTicks; axis equal
I2S=I2.*SensMB(:,:,2,2);
I2SSN=I2S.*SN;
I2SSN=I2S;
PadI2SSN=padarray(I2SSN,Sz128,'post');
FPadI2SSN=fft2(PadI2SSN);
HFPadI2SSN=fftkerns(:,:,1,2,:).*FPadI2SSN;
FHFPadI2SSN=ifft2(HFPadI2SSN);
CFHFPadI2SSN=permute(TakeTopLeftBlock(FHFPadI2SSN,Sz128),[1 2 3 5 4]);
SCFHFPadI2SSN=CFHFPadI2SSN.*conj(SensMB(:,:,2,:));
fgmontage(I2S);removeTicks; axis equal
fgmontage(I2SSN);removeTicks; axis equal
fgmontage(PadI2SSN);removeTicks; axis equal
fgmontage(fftshift(log(abs(FPadI2SSN))));removeTicks; axis equal
fgmontage(fftshift(log(abs(HFPadI2SSN))));removeTicks; axis equal
fgmontage(circshift(FHFPadI2SSN,[64 64]));removeTicks;axis equal
fgmontage(CFHFPadI2SSN);removeTicks; axis equal
fgmontage(SCFHFPadI2SSN);removeTicks; axis equal



% RI1=
%% images
fgmontage(RI1);axis equal
whos SN
fgmontage(RI1.*SN);axis equal
fgmontage(SN);axis equal
fgmontage(SN);removeTicks; axis equal
fgmontage(log(abs(fft2cg(padarray(RI1.*SN,Sz128,'both')))));axis equal
fgmontage((((padarray(RI1.*SN,Sz128,'both')))));axis equal
fgmontage((((padarray(RI1.*SN,Sz128,'post')))));axis equal
fgmontage(log(abs(fft2cg(padarray(RI1.*SN,Sz128,'post')))));axis equal
fgmontage(fftkernClean);axis equal
The table is broken
10:24
I don’t want to let them go before declaring
fgmontage(abs(fftkernClean));axis equal
fgmontage(log(abs(fftkernClean)));axis equal
fgmontage(fftshift(log(abs(fftkernClean))));axis equal
fgmontage(ifft2cg(fftkernClean));axis equal
fgmontage(abs(ifft2cg(fftkernClean)));axis equal
fgmontage(log(abs(ifft2cg(fftkernClean))));axis equal
fgmontage(log(abs(ifft2cg(fftkernClean))));removeTicks;axis equal
% fgmontage(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'));axis equal
fgmontage(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post')));axis equal
fgmontage(abs(log(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post')))));removeTicks;axis equal
fgmontage(fftshift(abs(log(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'))))));removeTicks;axis equal
close all
fgmontage(((ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'))))));removeTicks;axis equal
fgmontage(fftshift((ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'))))));removeTicks;axis equal
fgmontage(((ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'))))));removeTicks;axis equal
fgmontage(circshift((ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post')))),[64 64]));removeTicks;axis equal
fgmontage(TakeTopLeftBlock(ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post')))));removeTicks;axis equal
fgmontage(TakeTopLeftBlock(ifft2(fftkernClean.*fft2(padarray(RI1.*SN,Sz128,'post'))),Sz128));removeTicks;axis equal
%% Padded fft2 time test
nFFTSpeedTest=10000;
% SzStart=Sz128;SzEnd=Sz128*2;
A=rand(SzStart)+1i*rand(SzStart);
A = gpuArray(rand(SzStart)+1i*rand(SzStart));
tic
for i=1:nFFTSpeedTest
    B=fft2(A,SzEnd(1),SzEnd(2));
end
toc
tic
for i=1:nFFTSpeedTest
    BN=fftn(A,SzEnd);
end
toc
tic
for i=1:nFFTSpeedTest
%     BS1=fft(A,Sz128(1)*2,1);
%     BS2=fft(BS1,Sz128(2)*2,2);
    BS2=fft(fft(A,SzEnd(1),1),SzEnd(2),2);
end
toc
tic
for i=1:nFFTSpeedTest
    BP=fftn(padarray(A,Sz128,'post'));
end
toc
tic
for i=1:nFFTSpeedTest
    BP2=fft2(padarray(A,Sz128,'post'));
end
toc
%% Truncation
nFFTSpeedTest=10000;
SzStart=Sz128*2;SzEnd=Sz128;
% A=rand(SzStart)+1i*rand(SzStart);
A = gpuArray(rand(SzStart)+1i*rand(SzStart));
tic
for i=1:nFFTSpeedTest
    BN=TakeTopLeftBlock(fftn(A),SzEnd);
end
toc
tic
for i=1:nFFTSpeedTest
    B2=TakeTopLeftBlock(fft2(A),SzEnd);
end
toc
tic
for i=1:nFFTSpeedTest
%     BS1=fft(A,[],1);
%     BS2=TakeTopLeftBlock(fft(BS1(1:SzEnd(1),:),[],2),SzEnd);
    BS=TakeTopLeftBlock(fft(TakeTopLeftBlock(fft(A,[],1),[SzEnd(1),SzStart(2)]),[],2),SzEnd);
end
toc
%% Padded fft2 time test
nFFTSpeedTest=2000;
SzStart=Sz128;SzEnd=Sz128*2;
% A=rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]);
A = gpuArray(rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]));
tic
B=fft2(A,SzEnd(1),SzEnd(2));
toc
tic
%     BS1=fft(A,Sz128(1)*2,1);
%     BS2=fft(BS1,Sz128(2)*2,2);
BS2=fft(fft(A,SzEnd(1),1),SzEnd(2),2);
toc
%% Truncated fft2 time test
nFFTSpeedTest=1000;
SzStart=Sz128*2;SzEnd=Sz128;
% A=rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]);
A = gpuArray(rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]));
tic
B=TakeTopLeftBlock(fft2(A),SzEnd);
toc
tic
BS=TakeTopLeftBlock(fft(TakeTopLeftBlock(fft(A,[],1),[SzEnd(1),SzStart(2)]),[],2),SzEnd);
toc
%%
nFFTSpeedTest=1000;
SzStart=Sz128*2;SzEnd=Sz128;
% A=rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]);
A = gpuArray(rand([SzStart nFFTSpeedTest])+1i*rand([SzStart nFFTSpeedTest]));
tic
B=fft(A,[],2);
toc
tic
B1=fft(A(1:40,:,:),[],2);
toc