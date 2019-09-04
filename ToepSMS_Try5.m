A=load('/autofs/space/daisy_002/users/Gilad/All_Orientation-0x.mat');
load('/autofs/space/daisy_002/users/Gilad/CCSensMaps.mat');

permute54=@(x) permute(x,[1 2 3 5 4 6:12]);
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
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2, Sz128),Sz128, [6 6], Sz128*2,Sz128/2); % , [0 0] st.om

SN=nufftStruct.sn;
try
    P=nufftStruct.p.G;
catch
    P=nufftStruct.p;
end

fftkern=NUFFT_to_Toep(nufftStruct,wiforblock1);
fftkern=single(real(fftkern));


Fac=prod(Sz128)*4;

% TakeTopLeftBlock=@(x,Sz) x(1:Sz(1),1:Sz(2));
TakeTopLeftBlock=@(x,Sz) x(1:Sz(1),1:Sz(2),:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
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
% Dims: [H,W,Channels,MB,MBaux]
IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,[100 104]).*exp(1i*A.CurSetAll(2,:,:,[100 104])),[2 3 1 4]),Sz128));
SensMB=double(permute(SensCC(1:8,:,:,[10 20]),[2 3 1 4]));

w1=exp(1i*cCAIPIVecX*0).';
w2=exp(1i*cCAIPIVecX).';

wC=cat(3,w1,w2);

clear fftkerns
for i=1:size(IMB,4)
    for j=1:size(IMB,4)
        fftkerns(:,:,1,i,j)=NUFFT_to_Toep(nufftStruct,wC(:,1,i).*conj(wC(:,1,j)));
    end
end

MultSpTensor=@(Sp,x,Sz) reshape(Sp*reshape(x,size(Sp,2),[]),Sz);
fft2os=@(x) fft2(x,size(x,1)*2,size(x,2)*2);
ifft2os=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2]);
fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));
% ifft2osN(fft2osN(magic(5)))-magic(5)
%
NUFT_forwSMB=@(x,SensMB) MultSpTensor(P,fft2osN((x.*SN).*SensMB),[size(P,1),gsize(SensMB,3:4)]);
NUFT_adjSMB=@(x,SensMB) sum((ifft2osN(MultSpTensor(P',x,[gsize(SensMB,1:2)*2,gsize(SensMB,3:4)])).*conj(SN)).*conj(SensMB),3);

NUFT_forwSMBC=@(x,SensMB,wa) sum(NUFT_forwSMB(x,SensMB).*wa,3);
NUFT_adjSMBC=@(x,SensMB,wb) NUFT_adjSMB(x.*wb,SensMB);

NUFT_adj_forwSMBC=@(x,SensMB,wa,wb) NUFT_adjSMBC(NUFT_forwSMBC(x,SensMB,wa),SensMB,wb);

% That's L^2
% NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(ifft2osN(fkerns.*fft2osN(x.*SensMB)).*conj(SensMB),3),4),[1 2 3 5 4]);
NUFT_TSMBC_L2=@(x,SensMB,fkerns) sum(permute54(sum(ifft2osN(fkerns.*fft2osN(x.*SensMB)),4)).*conj(SensMB),3);

% That's O(L) + L^2 Hadamards
NUFT_TSMBC_L=@(x,SensMB,fkerns) sum(ifft2osN(permute54(sum(fkerns.*fft2osN(x.*SensMB),4))).*conj(SensMB),3);

% Verify operator normalization
ImSize=[gsize(SensMB,1:2) 1 size(SensMB,4)];
DataSize=[size(P,1) gsize(SensMB,3:4)];
Op=OpFromFunc(@(x) NUFT_forwSMB(x,SensMB),@(x) NUFT_adjSMB(x,SensMB));
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
disp(['Operator conj test: ' num2str(Out)]);
%

GGvSMBC=NUFT_adj_forwSMBC(IMB,SensMB,wC,conj(wC));
TvSMBCL2=NUFT_TSMBC_L2(IMB,SensMB,fftkerns);
TvSMBCL=NUFT_TSMBC_L(IMB,SensMB,fftkerns);

rDiffSMBC=abs(GGvSMBC(:)-TvSMBCL2(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBCL2(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDiffSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100

rDiffSMBC=abs(GGvSMBC(:)-TvSMBCL(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBCL(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDiffSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100
%%
FF=fft2osN(IMB.*SensMB);
FFH=fftkerns.*FF;
IFFH=ifft2osN(FFH);
SnIFFH=IFFH.*conj(permute54(SensMB));
SSnIFFH=sum(SnIFFH,3);
SSSnIFFH=permute54(sum(SSnIFFH,4));

mean(abs(GGvSMBC(:)-SSSnIFFH(:))/mean(abs(GGvSMBC(:))))*100

SIFFH=permute54(sum(IFFH,4));
SnSIFFH=sum(SIFFH.*conj(SensMB),3);

mean(abs(GGvSMBC(:)-SnSIFFH(:))/mean(abs(GGvSMBC(:))))*100

SFFH=permute(sum(FFH,4),[1 2 3 5 4]);
ISFFH=ifft2osN(SFFH);
SnISFFH=sum(ISFFH.*conj(SensMB),3);

mean(abs(GGvSMBC(:)-SnISFFH(:))/mean(abs(GGvSMBC(:))))*100

fftkerns3=fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,1)=fftkerns(:,:,1,1,2)-fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,2)=fftkerns(:,:,1,2,1)-fftkerns(:,:,1,1,1);
fftkerns3(:,:,1,3)=fftkerns(:,:,1,1,1);
% tic
EFF=cat(4,FF,sum(FF,4));
HEFF=EFF.*fftkerns3;
XHEFF=HEFF(:,:,:,[2 1])+HEFF(:,:,:,3);
% IXHEFF=Truncated2DIFT(XHEFF);
IXHEFF=ifft2osN(XHEFF);
SnIXHEFF=sum(IXHEFF.*conj(SensMB),3);
% toc

mean(abs(GGvSMBC(:)-SnIXHEFF(:))/mean(abs(GGvSMBC(:))))*100

% tic
TvSMBC=NUFT_TSMBC(IMB,SensMB,fftkerns);
% toc

rDiffSMBC=abs(GGvSMBC(:)-TvSMBC(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBC(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDiffSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100

QQ=(GGvSMBC(:)./TvSMBC(:));
mean(QQ(isfinite(QQ)))
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
for MB=1:4
    disp(MB);
    clear fftkernsX
    for i=1:MB
        for j=1:MB
            wi=randn(nTraj,1)+1i*randn(nTraj,1);
            wj=randn(nTraj,1)+1i*randn(nTraj,1);
            fftkernsX(:,:,1,i,j)=NUFFT_to_Toep(nufftStruct,wi.*conj(wj));
        end
    end
    
%     SensMBX=repmat(Sens,[1 1 1 MB]);
    SensMBX=double(permute(SensCC(:,:,:,10+(1:MB)),[2 3 1 4]));
    
%     IMBX=repmat(I1,[1 1 1 MB]);
    IIdxs=100+(1:MB);
    IMBX=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,IIdxs).*exp(1i*A.CurSetAll(2,:,:,IIdxs)),[2 3 1 4]),Sz128));
%     wCX=repmat(w1,[1 1 MB]);
    wCX=rand(nTraj,1,MB)+1i*rand(nTraj,1,MB);
    
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
% Dims: [H,W,Channels,MB,MBaux,TS]
MB=1;
nTimeSegmetns=94;
nChannels=3;
SliIdxs=[100 104 105:120];
SensIdx=[10 20 30:10:80];
CurSliIdx=SliIdxs(1:MB);
CurSensIdx=SensIdx(1:MB);
IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,CurSliIdx).*exp(1i*A.CurSetAll(2,:,:,CurSliIdx)),[2 3 1 4]),Sz128));
SensMB=double(permute(SensCC(1:nChannels,:,:,CurSensIdx),[2 3 1 4]));

T2Map_ms=imresizeBySlices(double(permute(A.CurSetAll(3,:,:,CurSliIdx),[2 3 1 4])),Sz128);
T2Map_ms=max(T2Map_ms,20);
B0Map_Hz=imresizeBySlices(double(permute(A.CurSetAll(4,:,:,CurSliIdx),[2 3 1 4])),Sz128);
B0Map_Hz(~isfinite(B0Map_Hz))=0;
B0Map_Hz(abs(IMB)<10)=0;
deltaT_ms=2.5/1000;
TimePoints_ms=(0:(nTraj-1))*deltaT_ms;
TimeSegmentsPoints_ms=permute(linspace(0,TimePoints_ms(end),nTimeSegmetns),[1 6 5 4 3 2]);
TSCM=exp(-TimeSegmentsPoints_ms./T2Map_ms);
TSCP=exp(1i*2*pi*B0Map_Hz.*TimeSegmentsPoints_ms/1000);
TSC=TSCM.*TSCP;
TSB=permute(GetTSCoeffsByLinear(nTraj,nTimeSegmetns),[1 6 5 4 3 2]);

% TSC=ones(size(TSC));
% TSC=TSCM;

SensMBTSC=SensMB.*TSC;

w1=exp(1i*cCAIPIVecX*0).';
w2=exp(1i*cCAIPIVecX).';

wC=cat(3,w1,w2);
wC=wC(:,:,1:MB);
wCTS=wC.*TSB;
fftkernsTS=NUFFT_to_Toep(nufftStruct,permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB);


MultSpTensor=@(Sp,x,Sz) reshape(Sp*reshape(x,size(Sp,2),[]),Sz);
fft2os=@(x) fft2(x,size(x,1)*2,size(x,2)*2);
ifft2os=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2]);
fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));

NUFT_forwSMBTS=@(x,SensMBTSC) MultSpTensor(P,fft2osN((x.*SN).*SensMBTSC),[size(P,1),gsize(SensMBTSC,3:4),1,1,size(SensMBTSC,6)]);
NUFT_adjSMBTS=@(x,SensMBTSC) sum(sum(ifft2osN(MultSpTensor(P',x,[gsize(SensMBTSC,1:2)*2,gsize(SensMBTSC,3:6)])).*conj(SN).*conj(SensMBTSC),3),6);

NUFT_forwSMBCTS=@(x,SensMBTSC,wa) sum(sum(NUFT_forwSMBTS(x,SensMBTSC).*wa,3),6);
NUFT_adjSMBCTS=@(x,SensMBTSC,wb) NUFT_adjSMBTS(x.*wb,SensMBTSC);

NUFT_adj_forwSMBCTS=@(x,SensMBTSC,wa,wb) NUFT_adjSMBCTS(NUFT_forwSMBCTS(x,SensMBTSC,wa),SensMBTSC,wb);

% That's L^2
NUFT_TSMBCTS_L2=@(x,SensMBTSC,fkerns) sum(sum(permute54(sum(ifft2osN(fkerns.*fft2osN(x.*SensMBTSC)),4)).*conj(SensMBTSC),3),6);

% That's O(L) + L^2 Hadamards
NUFT_TSMBCTS_L=@(x,SensMBTSC,fkerns) sum(sum(ifft2osN(permute54(sum(fkerns.*fft2osN(x.*SensMBTSC),4))).*conj(SensMBTSC),3),6);

% Verify operator normalization
ImSize=[gsize(SensMB,1:2) 1 size(SensMB,4)];
DataSize=[size(P,1) gsize(SensMB,3)];
% Op=OpFromFunc(@(x) NUFT_forwSMB(x,SensMB),@(x) NUFT_adjSMB(x,SensMB));
Op=OpFromFunc(@(x) NUFT_forwSMBCTS(x,SensMBTSC,wCTS),@(x) NUFT_adjSMBCTS(x,SensMBTSC,conj(wCTS)));
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
disp(['Operator conj test: ' num2str(Out)]);

GGvSMBCTS=NUFT_adj_forwSMBCTS(IMB,SensMBTSC,wCTS,conj(wCTS));
TvSMBCTSL2=NUFT_TSMBCTS_L2(IMB,SensMBTSC,fftkernsTS);
TvSMBCTSL=NUFT_TSMBCTS_L(IMB,SensMBTSC,fftkernsTS);

rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC) 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100

rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC)] 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100
%% Speed test
nSpeesTest=10;
nTimeSegmetns=10;
nCh=8;
% MB=5;
for MB=9:16
    disp(MB);
    SliIdxs=[100 104 105:120];
    SensIdx=[10 20 30:10:80 81:87];
    CurSliIdx=SliIdxs(1:MB);
    CurSensIdx=SensIdx(1:MB);
    
    IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,CurSliIdx).*exp(1i*A.CurSetAll(2,:,:,CurSliIdx)),[2 3 1 4]),Sz128));
    SensMB=double(permute(SensCC(1:nCh,:,:,CurSensIdx),[2 3 1 4]));
    
    T2Map_ms=imresizeBySlices(double(permute(A.CurSetAll(3,:,:,CurSliIdx),[2 3 1 4])),Sz128);
    T2Map_ms=max(T2Map_ms,20);
    B0Map_Hz=imresizeBySlices(double(permute(A.CurSetAll(4,:,:,CurSliIdx),[2 3 1 4])),Sz128);
    B0Map_Hz(~isfinite(B0Map_Hz))=0;
    B0Map_Hz(abs(IMB)<10)=0;
    deltaT_ms=2.5/1000;
    TimePoints_ms=(0:(nTraj-1))*deltaT_ms;
    TimeSegmentsPoints_ms=permute(linspace(0,TimePoints_ms(end),nTimeSegmetns),[1 6 5 4 3 2]);
    TSCM=exp(-TimeSegmentsPoints_ms./T2Map_ms);
    TSCP=exp(1i*2*pi*B0Map_Hz.*TimeSegmentsPoints_ms/1000);
    TSC=TSCM.*TSCP;
    TSB=permute(GetTSCoeffsByLinear(nTraj,nTimeSegmetns),[1 6 5 4 3 2]);

    SensMBTSC=SensMB.*TSC;

    w1=exp(1i*cCAIPIVecX*0).';
    w2=exp(1i*cCAIPIVecX).';

%     wC=cat(3,w1,w2);
%     wC=wC(:,:,1:MB);
    wC=rand(nTraj,1,MB)+1i*rand(nTraj,1,MB);
    wCTS=wC.*TSB;
    fftkernsTS=NUFFT_to_Toep(nufftStruct,permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB);
    
    tic
    for i=1:nSpeesTest
        GGvSMBCTS=NUFT_adj_forwSMBCTS(IMB,SensMBTSC,wCTS,conj(wCTS));
    end
    tFB(MB)=toc;
    tic
    for i=1:nSpeesTest
        TvSMBCTSL2=NUFT_TSMBCTS_L2(IMB,SensMBTSC,fftkernsTS);
    end
    tToepL2(MB)=toc;
    tic
    for i=1:nSpeesTest
        TvSMBCTSL2=NUFT_TSMBCTS_L(IMB,SensMBTSC,fftkernsTS);
    end
    tToepL(MB)=toc;
end
disp('Finished speed test');
%%
figure;plot(tFB,'k','LineWidth',3);hold on;plot(tToepL2,'b','LineWidth',3);plot(tToepL,'g','LineWidth',3);
legend({'Forw-adjoint','Toeplitz','Toeplitz acc'},'FontSize',22);
xlabel(MB);
ylabel('Seconds');
title(['Time for ' num2str(nSpeesTest) ' runs, ' num2str(nCh) ' channels, ', num2str(nTimeSegmetns) ' time segmetns, CPU'],'FontSize',22);
%% Now with time segments NOT TO USE
% Dims: [H,W,Channels,MB,MBaux,TS,batch]
MB=2;
nTimeSegmetns=4;
nChannels=3;
SliIdxs=[100 104 105:120];
SensIdx=[10 20 30:10:80];
CurSliIdx=SliIdxs(1:MB);
CurSensIdx=SensIdx(1:MB);
IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,CurSliIdx).*exp(1i*A.CurSetAll(2,:,:,CurSliIdx)),[2 3 1 4]),Sz128));
SensMB=double(permute(SensCC(1:nChannels,:,:,CurSensIdx),[2 3 1 4]));

% IMB=repmat(IMB,[1 1 1 1 1 1 2 3]);
% SensMB=repmat(SensMB,[1 1 1 1 1 1 2 3]);

T2Map_ms=imresizeBySlices(double(permute(A.CurSetAll(3,:,:,CurSliIdx),[2 3 1 4])),Sz128);
T2Map_ms=max(T2Map_ms,20);

B0Map_Hz=imresizeBySlices(double(permute(A.CurSetAll(4,:,:,CurSliIdx),[2 3 1 4])),Sz128);
B0Map_Hz(~isfinite(B0Map_Hz))=0;
% B0Map_Hz(abs(IMB)<10)=0;
deltaT_ms=2.5/1000;
TimePoints_ms=(0:(nTraj-1))*deltaT_ms;
TimeSegmentsPoints_ms=permute(linspace(0,TimePoints_ms(end),nTimeSegmetns),[1 6 5 4 3 2]);
TSCM=exp(-TimeSegmentsPoints_ms./T2Map_ms);
TSCP=exp(1i*2*pi*B0Map_Hz.*TimeSegmentsPoints_ms/1000);
TSC=TSCM.*TSCP;
TSB=permute(GetTSCoeffsByLinear(nTraj,nTimeSegmetns),[1 6 5 4 3 2]);


SensMBTSC=SensMB.*TSC;
SensMBTSC=repmat(SensMBTSC,[1 1 1 1 1 1 3]);
IMB=repmat(IMB,[1 1 1 1 1 1 3]);

w1=exp(1i*cCAIPIVecX*0).';
w2=exp(1i*cCAIPIVecX).';

wC=cat(3,w1,w2);
wC=wC(:,:,1:MB);
wCTS=wC.*TSB;
% fftkernsTS=NUFFT_to_Toep(nufftStruct,permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB);
wCTSE0=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*TSB;
wCTSE1=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*circshift(TSB,1,6);
wCTSEm1=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*circshift(TSB,-1,6);
wCTSE=cat(7,wCTSE0,wCTSE1,wCTSEm1);
fftkernsTS=NUFFT_to_Toep(nufftStruct,permute(wCTSE,[1 2 3 4 5 8 6 7]));

MultSpTensor=@(Sp,x,Sz) reshape(Sp*reshape(x,size(Sp,2),[]),Sz);
fft2os=@(x) fft2(x,size(x,1)*2,size(x,2)*2);
ifft2os=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2]);
fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));

NUFT_forwSMBTS=@(x,SensMBTSC) MultSpTensor(P,fft2osN((x.*SN).*SensMBTSC),[size(P,1),gsize(SensMBTSC,3:4),1,1,gsize(SensMBTSC,6:9)]);
NUFT_adjSMBTS=@(x,SensMBTSC) sum(sum(ifft2osN(MultSpTensor(P',x,[gsize(SensMBTSC,1:2)*2,gsize(SensMBTSC,3:9)])).*conj(SN).*conj(SensMBTSC),3),6);

NUFT_forwSMBCTS=@(x,SensMBTSC,wa) sum(sum(NUFT_forwSMBTS(x,SensMBTSC).*wa,3),6);
NUFT_adjSMBCTS=@(x,SensMBTSC,wb) NUFT_adjSMBTS(x.*wb,SensMBTSC);

NUFT_adj_forwSMBCTS=@(x,SensMBTSC,wa,wb) NUFT_adjSMBCTS(NUFT_forwSMBCTS(x,SensMBTSC,wa),SensMBTSC,wb);

% That's L^2
NUFT_TSMBCTS_L2=@(x,SensMBTSC,fkerns) sum(sum(permute54(sum(ifft2osN(fkerns.*fft2osN(x.*SensMBTSC)),4)).*conj(SensMBTSC),3),6);

% That's O(L) + L^2 Hadamards
NUFT_TSMBCTS_L=@(x,SensMBTSC,fkerns) sum(sum(ifft2osN(permute54(sum(fkerns.*fft2osN(x.*SensMBTSC),4))).*conj(SensMBTSC),3),6);

% Verify operator normalization
ImSize=[gsize(SensMB,1:2) 1 size(SensMB,4) 1 1 gsize(SensMBTSC,7:9)];
DataSize=[size(P,1) gsize(SensMB,3) 1 1 1 1 gsize(SensMBTSC,7:9)];
% Op=OpFromFunc(@(x) NUFT_forwSMB(x,SensMB),@(x) NUFT_adjSMB(x,SensMB));
Op=OpFromFunc(@(x) NUFT_forwSMBCTS(x,SensMBTSC,wCTS),@(x) NUFT_adjSMBCTS(x,SensMBTSC,conj(wCTS)));
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
disp(['Operator conj test: ' num2str(Out)]);

GGvSMBCTS=NUFT_adj_forwSMBCTS(IMB,SensMBTSC,wCTS,conj(wCTS));
SensMBTSC=permute(SensMBTSC,[1 2 3 4 5 8 6 7]);
IMB=permute(IMB,[1 2 3 4 5 6 8 7]);

TvSMBCTSL2=NUFT_TSMBCTS_L2(IMB,SensMBTSC,fftkernsTS);
TvSMBCTSL=NUFT_TSMBCTS_L(IMB,SensMBTSC,fftkernsTS);

TSCx=permute(TSC,[1 2 3 4 5 7 6]);
TSCx1=circshift(TSCx,1,7);
TSCxm1=circshift(TSCx,-1,7);
TSCxE=cat(8,TSCx,TSCx1,TSCxm1);
TvSMBCTSL2=TvSMBCTSL2.*conj(TSCxE)./conj(TSCx) ;
TvSMBCTSL2=gsum(TvSMBCTSL2,7:8);
GGvSMBCTS=GGvSMBCTS(:,:,:,:,:,:,1);
%%
rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC) 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100

rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC)] 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100
%% Now with TSaux
% Dims: [H,W,Channels,MB,MBaux,TS,TSaux,batch]
MB=2;
nTimeSegmetns=7;
nChannels=8;
SliIdxs=[100 104 105:120];
SensIdx=[10 20 30:10:80];
CurSliIdx=SliIdxs(1:MB);
CurSensIdx=SensIdx(1:MB);
IMB=double(imresizeBySlices(permute(A.CurSetAll(1,:,:,CurSliIdx).*exp(1i*A.CurSetAll(2,:,:,CurSliIdx)),[2 3 1 4]),Sz128));
SensMB=double(permute(SensCC(1:nChannels,:,:,CurSensIdx),[2 3 1 4]));

T2Map_ms=imresizeBySlices(double(permute(A.CurSetAll(3,:,:,CurSliIdx),[2 3 1 4])),Sz128);
T2Map_ms=max(T2Map_ms,20);

B0Map_Hz=imresizeBySlices(double(permute(A.CurSetAll(4,:,:,CurSliIdx),[2 3 1 4])),Sz128);
B0Map_Hz(~isfinite(B0Map_Hz))=0;
% B0Map_Hz(abs(IMB)<10)=0;
deltaT_ms=2.5/1000;
TimePoints_ms=(0:(nTraj-1))*deltaT_ms;
TimeSegmentsPoints_ms=permute(linspace(0,TimePoints_ms(end),nTimeSegmetns),[1 6 5 4 3 2]);
TSCM=exp(-TimeSegmentsPoints_ms./T2Map_ms);
TSCP=exp(1i*2*pi*B0Map_Hz.*TimeSegmentsPoints_ms/1000);
TSC=TSCM.*TSCP;
TSB=permute(GetTSCoeffsByLinear(nTraj,nTimeSegmetns),[1 6 5 4 3 2]);

SensMBTSC=SensMB.*TSC;

w1=exp(1i*cCAIPIVecX*0).';
w2=exp(1i*cCAIPIVecX).';

wC=cat(3,w1,w2);
wC=wC(:,:,1:MB);
wCTS=wC.*TSB;
% fftkernsTS=NUFFT_to_Toep(nufftStruct,permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB);
wCTSE0=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*TSB;
wCTSE1=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*circshift(TSB,1,6);
wCTSEm1=permute(wC,[1 2 4 3]).*permute(conj(wC),[1 2 5 4 3]).*TSB.*circshift(TSB,-1,6);
wCTSE=cat(7,wCTSE0,wCTSE1,wCTSEm1);
fftkernsTS=NUFFT_to_Toep(nufftStruct,wCTSE);

TSCx1=circshift(TSC,1,6);
TSCxm1=circshift(TSC,-1,6);
TSCxE=cat(7,TSC,TSCx1,TSCxm1);

cSensMBTSC=conj(SensMB.*TSCxE);

MultSpTensor=@(Sp,x,Sz) reshape(Sp*reshape(x,size(Sp,2),[]),Sz);
fft2os=@(x) fft2(x,size(x,1)*2,size(x,2)*2);
ifft2os=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2]);
fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));

NUFT_forwSMBTS=@(x,SensMBTSC) MultSpTensor(P,fft2osN((x.*SN).*SensMBTSC),[size(P,1),gsize(SensMBTSC,3:4),1,1,gsize(SensMBTSC,6:9)]);
NUFT_adjSMBTS=@(x,SensMBTSC) sum(sum(ifft2osN(MultSpTensor(P',x,[gsize(SensMBTSC,1:2)*2,gsize(SensMBTSC,3:9)])).*conj(SN).*conj(SensMBTSC),3),6);

NUFT_forwSMBCTS=@(x,SensMBTSC,wa) sum(sum(NUFT_forwSMBTS(x,SensMBTSC).*wa,3),6);
NUFT_adjSMBCTS=@(x,SensMBTSC,wb) NUFT_adjSMBTS(x.*wb,SensMBTSC);

NUFT_adj_forwSMBCTS=@(x,SensMBTSC,wa,wb) NUFT_adjSMBCTS(NUFT_forwSMBCTS(x,SensMBTSC,wa),SensMBTSC,wb);

% That's L^2
NUFT_TSMBCTS_L2=@(x,SensMBTSC,fkerns,cSensMBTSC) sum(sum(sum(permute54(sum(ifft2osN(fkerns.*fft2osN(x.*SensMBTSC)),4)).*cSensMBTSC,3),6),7);

% That's O(L) + L^2 Hadamards
NUFT_TSMBCTS_L=@(x,SensMBTSC,fkerns,cSensMBTSC) sum(sum(sum(ifft2osN(permute54(sum(fkerns.*fft2osN(x.*SensMBTSC),4))).*cSensMBTSC,3),6),7);

% Verify operator normalization
ImSize=[gsize(SensMB,1:2) 1 size(SensMB,4) 1 1 gsize(SensMBTSC,7:9)];
DataSize=[size(P,1) gsize(SensMB,3) 1 1 1 1 gsize(SensMBTSC,7:9)];
% Op=OpFromFunc(@(x) NUFT_forwSMB(x,SensMB),@(x) NUFT_adjSMB(x,SensMB));
Op=OpFromFunc(@(x) NUFT_forwSMBCTS(x,SensMBTSC,wCTS),@(x) NUFT_adjSMBCTS(x,SensMBTSC,conj(wCTS)));
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
disp(['Operator conj test: ' num2str(Out)]);

GGvSMBCTS=NUFT_adj_forwSMBCTS(IMB,SensMBTSC,wCTS,conj(wCTS));
%
TvSMBCTSL2=NUFT_TSMBCTS_L2(IMB,SensMBTSC,fftkernsTS,cSensMBTSC);
TvSMBCTSL=NUFT_TSMBCTS_L(IMB,SensMBTSC,fftkernsTS,cSensMBTSC);

rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL2(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC) 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100

rDiffSMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))./abs(GGvSMBCTS(:));
rDiff2SMBC=abs(GGvSMBCTS(:)-TvSMBCTSL(:))/mean(abs(GGvSMBCTS(:)));
% [max(rDiffSMBC) mean(rDiffSMBC)] 
[max(rDiff2SMBC) mean(rDiff2SMBC)]*100

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

% Shortcut
SHFPIMBWithS=sum(HFPIMBWithS,5);
FSHFPIMBWithS=ifft2(SHFPIMBWithS);
CFSHFPIMBWithS=TakeTopLeftBlock(FSHFPIMBWithS,Sz128);

fgmontage(fftshift(log(abs(SHFPIMBWithS))));removeTicks; axis equal
fgmontage(circshift(FSHFPIMBWithS,[64 64]));removeTicks;axis equal
fgmontage(CFSHFPIMBWithS);removeTicks;axis equal

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