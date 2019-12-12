345

I=double(rgb2gray(imread('beach.jpg')))/256;

S=[50 60];

I=imresize(I,S);

[X, Y]=ndgrid(-S(1)/2:S(1)/2 -1, -S(2)/2:S(2)/2 -1);
XY=cat(3,X,Y);
Idxs2=reshape(XY,[prod(S) 2]);

Coeffs=gFTCoeffs2D(S,Idxs2);
FTVec=squeeze(gsum(RepDotMult(I,Coeffs),1:2));
FTM=reshape(FTVec,S);

B=-mod(X+Y,2)*2+1;
FTMX=FTM.*B;
FFTRes=fft2cg(I);
gsum(abs(FTMX-FFTRes))
%%
% SR is image->Data

FT2DIxs=randperm(size(Coeffs,3));

CoeffsX=Coeffs(:,:,FT2DIxs(1:1500));
%%
N=prod(S);
R=2;
SpiralPow=1.3;
SpiralLoops=21;
C=gSpiral(N/R,SpiralLoops,SpiralPow,0);
% figure;plot(real(C),imag(C),'*-')
X=real(C)*0.5*S(1);
Y=imag(C)*0.5*S(2);
figure;plot(X,Y,'*-')

Idxs2=[X.' Y.'];
CoeffsX=gFTCoeffs2D(S,Idxs2);
%%
SR=RepDotMultOp(CoeffsX,1:2);

Data=SR*I;
XX=SR'*Data;

TVW=.001;
SensMFull=ones(S);
FigH=100;
nfnlCgIters=10;

ImSize=S;
DataSize=size(Data);


OperatorTest(SR,ImSize,DataSize,DataSize(3),1);
param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',SR,'Verbose',false,'XFM',1,'TV',TVOP_MSlice,'xfmWeight',0),init);
param.data =     Data;
res=zeros(ImSize);
if(isfield(param,'WarmStart'))
    res=param.WarmStart;
end
RunFnlViewAmp=1;
RunfnlCgIterations;

% im_res=RunfnlCgIterationsf(FigH,SR,nfnlCgIters,SensMFull,TVW,Data,2,struct('ShowFig',true));

figure; 
subplot(1,3,1);
imagesc(I);colormap gray;
title('Original');
subplot(1,3,2);
imagesc(abs(im_res));colormap gray;
title(['Reconstructed, R=' num2str(R)]);
xlabel(['TVW=' num2str(TVW)]);
subplot(1,3,3);
plot(X,Y,'*-')
title([num2str(SpiralLoops) ' loops, Pow=' num2str(SpiralPow)]);
axis square
%%
setenv('TOOLBOX_PATH','/home/a/bart-0.4.01X/')

DataP=gpermute(Data,[6 3]);
SensP=ones(1,50,60);
SRP=permute(CoeffsX,[5 1 2 4 6 3]);


BartRecon=bart('pics -l1 -r 1e-5 -u 1e-4 -m',DataP*1e0,SensP*1e0,SRP*1e0);
BartRecon=squeeze(BartRecon);
%

figure; 
subplot(2,2,1);
imagesc(I);colormap gray;
title('Original');
subplot(2,2,2);
plot(X,Y,'*-')
title([num2str(SpiralLoops) ' loops, Pow=' num2str(SpiralPow)]);
axis square
subplot(2,2,3);
imagesc(abs(im_res));colormap gray;
title(['Reconstructed, R=' num2str(R)]);
xlabel(['TVW=' num2str(TVW)]);

subplot(2,2,4);
imagesc(abs(BartRecon));colormap gray;
title(['BART, R=' num2str(R)]);
% xlabel(['TVW=' num2str(TVW)]);


%% Multi channel
load('DemoBrainAndMps.mat','DATA','maps');
I=fft2cg(DATA);
S=[80 80];

SensF=RunESPIRiTForSensMaps(I,12);

I=imresizeBySlices(I,S);

I=grmss(I,3);
Sens=imresizeBySlices(SensF,S);
Sens=I*0+1;

% Sens=RunESPIRiTForSensMaps(I,17);

%%
N=prod(S);
R=1;
SpiralPow=1.3;
SpiralLoops=17;
C=gSpiral(N/R,SpiralLoops,SpiralPow,0);
% figure;plot(real(C),imag(C),'*-')
X=real(C)*0.5*S(1);
Y=imag(C)*0.5*S(2);
figure;plot(X,Y,'*-')

Idxs2=[X.' Y.'];
CoeffsX=gFTCoeffs2D(S,Idxs2);
%% use Native bart
setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensorClean')
SensP=permute(Sens,[1 2 4 3]);

nChannels=size(I,3);

% Data=gsum(RepDotMult(I,CoeffsX3),1:2);
Data=gsum(RepDotMult(I,gpermute(CoeffsX,[4 3])),1:2);
DataP=permute(Data,[1 4 2 3]);
Trajectory=[X;Y;X*0];
TrajectoryP=Trajectory;
reco = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',TrajectoryP, DataP, SensP);
fgmontage(reco)
% reco = bart('pics -S -r:0.01 -RT:1024:0:0.0001 -p Weights  -t ',Trajectory, Data, SensP);
%%
AA=readcfl('/home/a/bart-gpu_tensorClean/AA');
%%
NN=74;

% Idxs=Idxs+0.5;
% NN=113;
Idxs=-NN/2:NN/2-1;

Idxs=-NN/2:0.27:NN/2;

Idxs=rand(1,101)*NN-NN/2;
% Idxs=sort(Idxs);

% Idxs=Idxs+0.2;
% Idxs=Idxs(2:2:end);

A=rand(NN,1)+ 1i*(rand(NN,1)-0.5);

IdxsB=Idxs;
if(mod(NN,2)==1)
    IdxsB=Idxs-0.5;
end
NUbyBA=bart('nufft ',[IdxsB;Idxs*0;Idxs*0],A);
if(mod(NN,2)==0)
    NUbyBA=NUbyBA.*exp(-1i* (Idxs*pi/NN));% linspace(pi/2,-pi/2,numel(Idxs)));
end

CoeffsY=gFTCoeffs2D([NN 1],[Idxs.' Idxs.'*0]);

NUbyM=squeeze(gsum(RepDotMult(A,CoeffsY),[1]));
% NUbyM=squeeze(gsum(RepDotMult(A,gpermute(CoeffsY,[4 3])),[1 3]));

% NUbyFT=fft1cg(A,1);

figure; 
subplot(2,2,1);
plot(abs(NUbyBA),'k.');hold on;
% plot(abs(NUbyFT),'b*');
plot(abs(NUbyM),'ro');
subplot(2,2,2);
plot(angle(NUbyBA),'k.');hold on;
% plot(angle(NUbyFT),'b*')
plot(angle(NUbyM),'ro')
subplot(2,2,3);
plot(real(NUbyBA),'k.');hold on;
% plot(real(NUbyFT),'b*')
plot(real(NUbyM),'ro')
subplot(2,2,4);
plot(imag(NUbyBA),'k.');hold on;
% plot(imag(NUbyFT),'b*')
plot(imag(NUbyM),'ro')
%%
% setenv('TOOLBOX_PATH','/home/a/bart-0.4.01X/')
% setenv('TOOLBOX_PATH','/home/a/bart-0.4.02MC/')
setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensor/')

% I=I(:,:,1:2);
% Sens=Sens(:,:,1:2);
% Sens=Sens*0+1;
nChannels=size(I,3);

CoeffsX3=RepDotMult(Sens,gpermute(CoeffsX,[4 3]));

% Data=gsum(RepDotMult(I,CoeffsX3),1:2);
Data=gsum(RepDotMult(I,gpermute(CoeffsX,[4 3])),1:2);

DataP=permute(Data,[1 2 5 3 6 7 8 4]);
% SensP=ones(1,S(1),S(2),nChannels);
SRP=permute(CoeffsX3,[1 2 5 3 6 7 8 4]);

% BartRecon=bart('pics -g -l1 -r 1e-7',DataP*1e5,SRP*1e0);

% BartRecon=bart('pics -l1 -r 1 -u 1e-3 -m',DataP*1e1,SRP*1e0);
BartRecon=bart('pics -l1 -r 1e-3 -u 1e-8 -m',DataP*1e1,SRP*1e0);
% BartRecon=bart('pics -g -l1 -r 1e-7 -u 1e-8 -m',DataP*1e1,SRP*1e0);

% BartRecon=bart('pics -l1 -r 1e-7 -u 1e-8 -m',DataP*1e1,SensP*1e0,SRP*1e0);
% BartRecon=bart('pics -l2 -r 1e2',DataP*1e0,SensP*1e0,SRP*1e2);
BartRecon=squeeze(BartRecon);
%

figure; 
subplot(2,2,1);
gmontage(I,'Size',[2 4]);
title('Original');
subplot(2,2,2);
gmontage(angle(I),'Size',[2 4]);
title('Phase');
subplot(2,2,3);
plot(X,Y,'*-')
title([num2str(SpiralLoops) ' loops, Pow=' num2str(SpiralPow)]);
axis square

subplot(2,2,4);
imagesc(abs(BartRecon));colormap gray;
title(['BART, R=' num2str(R)]);
% xlabel(['TVW=' num2str(TVW)]);


%% BART NUFT check
NN=62;
A=zeros(NN,1);
A=rand(NN,1)-0.1+ 1i*(rand(NN,1)-0.5);
%
Idxs=-NN/2:NN/2-1;

Idxs=Idxs+0.2;

Idxs=rand(1,50)*NN-NN/2;
% Idxs=-30:1:30;

FesNUFTOp = nuFTOperator([(  (Idxs+0).'*2*pi/NN) , Idxs.'*0],[NN 5]); % Orange


% A(31)=1;
% NUbyBA1=bart('nufft ',[Idxs;Idxs*0;Idxs*0],A);
% A=zeros(NN,1);
% A(32)=1;
% NUbyBA2=bart('nufft ',[Idxs;Idxs*0;Idxs*0],A);
% A=zeros(NN,1);
% A(33)=1;
if(mod(NN,2)==1)
    Idxs=Idxs-0.5;
end
NUbyBA3=bart('nufft ',[Idxs;Idxs*0;Idxs*0],A);
if(mod(NN,2)==1)
    NUbyBA3=NUbyBA3.*exp(1i* 2*((Idxs+0.5)*pi/NN));
end

NUbyFS3=(FesNUFTOp.'*repmat(A/sqrt(5),[1 5])).';

% NUbyFS3=NUbyFS3.*exp(-1i*linspace(-pi,pi,numel(NUbyFS3)));
QQ=[NUbyBA3; NUbyFS3];

figure;subplot(1,2,1);
plot(abs(NUbyBA3),'.');
hold on
plot(abs(NUbyFS3),'o');

subplot(1,2,2);
plot(angle(NUbyBA3),'.');
hold on
plot(angle(NUbyFS3),'o');

%% BART NUFT 2D check
Sz=[68 60];

A=rand(Sz)-0.5+ 1i*(rand(Sz)-0.5);
%
Idxs1=-Sz(1)/2:Sz(1)/2-1;
Idxs2=-Sz(2)/2:Sz(2)/2-1;

P1=getKrandomSamples(Sz(1),30);
P2=getKrandomSamples(Sz(2),30);
% P1=1:67;
% P2=1:67;

IdxsX=[Idxs1(P1); Idxs2(P2)*0;Idxs1(P1)*0];

% IdxsX=[-20:2:20; zeros(2,21)];
% IdxsX=[zeros(1,21);-20:2:20; zeros(1,21)];

IdxsX=[rand(1,101)*Sz(1)-Sz(1)/2;rand(1,101)*Sz(2)-Sz(2)/2;zeros(1,101)];

[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(IdxsX,Sz),Sz); % Orange
NUbyFS3=(FesNUFTOp.'*A).';

x = A .* st.nufftStruct.sn;		% apply scaling factors
Xk = col(fftn(x, st.nufftStruct.Kd));	% [*Kd] oversampled FFT, padded at end
X = st.nufftStruct.p * Xk;					% [M,*L]
NUbyFS3=X.'/sqrt(prod(Sz));

% Step1 = A .* st.nufftStruct.sn;		% apply scaling factors
% Step2=fftn(Step1, st.nufftStruct.Kd);
% Xk = col();	% [*Kd] oversampled FFT, padded at end
% X = st.nufftStruct.p * Xk;					% [M,*L]
% NUbyFS3=X.'/sqrt(prod(Sz));

Op=gNUFT(st.nufftStruct.sn,st.nufftStruct.Kd,st.nufftStruct.p/sqrt(prod(Sz)));

NUbyFS3=Op*A;

Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz));
save('ForTFNUFT.mat','SN','Kd','P','A','NUbyFS3');
CC=load('FromTFNUFT.mat');
C=CC.C;
figure;plot(abs(C),'o');hold on;plot(abs(NUbyFS3),'r.')

QQ=Op'*NUbyFS3;

NUbyBA3=gBARTnuFT(IdxsX,A,Sz);

CoeffsY=gFTCoeffs2D(Sz,[IdxsX(1,:).',IdxsX(2,:).']);

NUbyM=squeeze(gsum(RepDotMult(A,CoeffsY),[1:2]));
NUbyM=NUbyM.*exp(1i*pi*(IdxsX(1,:).')/Sz(1));
NUbyM=NUbyM.*exp(1i*pi*(IdxsX(2,:).')/Sz(2));

QQ=[NUbyBA3; NUbyFS3];

figure;subplot(1,2,1);
plot(abs(NUbyFS3),'ko');
hold on
plot(abs(NUbyM),'g*');
plot(abs(NUbyBA3),'r.');
title('Magnitude')
legend({'Fes','M','BART'});

subplot(1,2,2);
plot(angle(NUbyFS3),'ko');
hold on
plot(angle(NUbyM),'g*');
plot(angle(NUbyBA3),'r.');
title('Phase')
%%
% function  [A, s] = nuFTOperator(trajectory, imageDim, sensmaps, os, Neighborhood)

%%
ImSize=size(A);
x = randn(ImSize) + 1j*randn(ImSize);
DataSize=size(NUbyFS3);
y = randn([DataSize(1:2)]) + 1j*randn([DataSize(1:2)]);
A=Op;

Ax = A*x;
Aty = A'*y;
aa=x(:)'*Aty(:);
bb=conj(y(:)'*Ax(:));
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
if(nargin<5)
    disp(['Operator conj test: ' num2str(Out)]);
else
    Err=Out;
    Out=Out<1e-10;
    if(~Out)
%         error('Failed OperatorTest');
        disp(['Operator conj test: ' num2str(Err)]);
    end
end
%%
SzZ=[50 50];
I=rand(SzZ);
K=rand(SzZ);
aM=rand(SzZ)>0.8;
a=rand(SzZ).*aM;
C=conv2(I,K,'same');

b=conv2(a,K,'same');
aa1=gsum(C.*a)

aa2=gsum(I.*b)
%% From EPFL

if numel(FOV)==1
    FOV = FOV*[1 1];
end
param.k = w*diag(FOV)/2/pi;
param.res = size(x);
param.method = 'nuft gg';

m = E0(x,param);
        % MEX implementation of Greengard's fast gaussian gridding.
        % Complexity N^2log(N).
        if ~isfield(param,'st')
            param.st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], res, 12, 4*res);
        end
        m = nuft_gg_forw(x, param.st);
