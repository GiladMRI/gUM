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
Trajm2=BARTTrajx;
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

NUFT_forwS=@(x,Sens) P*reshape(fft2(padarray((x.*SN).*Sens,Sz128,'post')),prod(Sz128)*4,[])*Fac;
NUFT_adjS=@(x,Sens) sum((TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128(1)*2,Sz128(2)*2,[])),Sz128).*conj(SN)).*conj(Sens),3);

NUFT_adj_forwS=@(x,Sens,wa,wb) NUFT_adjS((NUFT_forwS(x,Sens).*wa).*wb,Sens);

NUFT_TS=@(x,Sens,fkern) sum(TakeTopLeftBlock(ifft2(fkern.*fft2(padarray(x.*Sens,Sz128,'post'))),Sz128).*conj(Sens),3);
save
GGvS=NUFT_adj_forwS(I1,Sens,wi.',conj(wi.'));
TvS=NUFT_TS(I1,Sens,fftkern);

rDiffS=abs(GGvS(:)-TvS(:))./abs(GGvS(:));
rDiff2S=abs(GGvS(:)-TvS(:))/mean(abs(GGvS(:)));
[max(rDiffS) mean(rDiffS) max(rDiff2S) mean(rDiff2S)]*100
%% With sensitivity maps, MB
SensMB=repmat(Sens,[1 1 1 2]);
IMB=cat(4,I1,I2);
wC=cat(3,w1,w2);

fftkerns(:,:,1,1,1)=fftkern11;
fftkerns(:,:,1,1,2)=fftkern12;
fftkerns(:,:,1,2,1)=fftkern21;
fftkerns(:,:,1,2,2)=fftkern22;

NUFT_forwSMB=@(x,SensMB) reshape(P*reshape(fft2(padarray((x.*SN).*SensMB,Sz128,'post')),prod(Sz128)*4,[])*Fac,[],size(SensMB,3),size(SensMB,4));
NUFT_adjSMB=@(x,SensMB) sum((TakeTopLeftBlock(ifft2(reshape(P'*x,Sz128(1)*2,Sz128(2)*2,size(SensMBsave,3),size(SensMB,4))),Sz128).*conj(SN)).*conj(SensMB),3);

NUFT_forwSMBC=@(x,SensMB,wa) sum(NUFT_forwSMB(x,SensMB).*wa,3);
NUFT_adjSMBC=@(x,SensMB,wb) NUFT_adjSMB(reshape(x.*wb,[],size(SensMB,3)*size(SensMB,4)),SensMB);

NUFT_adj_forwSMBC=@(x,SensMB,wa,wb) NUFT_adjSMBC(NUFT_forwSMBC(x,SensMB,wa),SensMB,wb);

NUFT_TSMBC=@(x,SensMB,fkerns) permute(sum(sum(TakeTopLeftBlock(ifft2(fkerns.*fft2(padarray(x.*SensMB,Sz128,'post'))),Sz128).*conj(SensMB),3),4),[1 2 3 5 4]);

GGvSMBC=NUFT_adj_forwSMBC(IMB,SensMB,wC,conj(wC));
TvSMBC=NUFT_TSMBC(IMB,SensMB,fftkerns);

rDiffSMBC=abs(GGvSMBC(:)-TvSMBC(:))./abs(GGvSMBC(:));
rDiff2SMBC=abs(GGvSMBC(:)-TvSMBC(:))/mean(abs(GGvSMBC(:)));
[max(rDiffSMBC) mean(rDifsavefSMBC) max(rDiff2SMBC) mean(rDiff2SMBC)]*100
%% CGP no MB, L2 Finite Difference regularization
Sig=NUFT_forw(I1);
AHA_I=NUFT_adj(SigS);


NUFT_adj_forw_clean=@(x) NUFT_adj(NUFT_forw(x));

fftkernClean=NUFFT_to_Toep(nufftStruct,ones(nTraj,1));
fftkernClean=single(real(fftkernClean));

Lambda=1e1;

Dx = @(x) diff(x,1,2);
Dy = @(x) diff(x,1,1);
Dxc=@(x) [Dx(x) zeros(128,1)]-[zeros(128,1) Dx(x)];save
Dyc=@(x) [Dy(x); zeros(1,128)]-[zeros(1,128); Dy(x)];

AHAReg=@(x) NUFT_adj_forw_clean(x)-Lambda*(Dxc(x)+Dyc(x));
% AHAReg=@(x) NUFT_T(x,fftkernClean)-Lambda*(Dxc(x)+Dyc(x));
tic;
x0=zeros(size(I1));
[x, ItersUsed] = cgp2D(x0, 0, 0, AHA_I, MaxIter, 10^-8, @(Z, x) AHAReg(x), @(Z, o) o);
toc
ShowAbsAngle(x)
%% Single band, with sens, cgp
SigS=NUFT_forwS(I1,Sens);
AHA_I=NUFT_adjS(SigS,Sens);

Lambda=1e-13;
save
Dx = @(x) diff(x,1,2);
Dy = @(x) diff(x,1,1);
Dxc=@(x) [Dx(x) zeros(128,1)]-[zeros(128,1) Dx(x)];
Dyc=@(x) [Dy(x); zeros(1,128)]-[zeros(1,128); Dy(x)];

% AHAReg=@(x) NUFT_adj_forwS(x,Sens,1,1)-Lambda*(Dxc(x)+Dyc(x));
AHAReg=@(x) NUFT_TS(x,Sens,fftkernClean)-Lambda*(Dxc(x)+Dyc(x));
tic;
x0=zeros(size(I1));
[x, ItersUsed] = cgp2D(x0, 0, 0, AHA_I, MaxIter, 10^-8, @(Z, x) AHAReg(x), @(Z, o) o);
toc
ShowAbsAngle(x)
%% Single-band, Sens, fnlCg with TV/Wavelet
DataP=SigS;save

AOdd = OpFromFunc(@(x) NUFT_forwS(x,Sens),@(x) NUFT_adjS(x,Sens));

TVW=1e-7;

% XFMStr='Daubechies';
XFMStr='None';

filterSize=4;
wavScale=4;

if(strcmp(XFMStr,'None'))save
    XFM=1;save
    xfmWeight=0;
else
    XFM = Wavelet(XFMStr,filterSize,wavScale);
    xfmWeight = 1e-7;	% Weight for Transform L1 penalty
end

param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight),init);

param.data = DataP;

param.ShowFig=false;

nfnlCgIters=40;
RunFnlViewAmp=1;save
res=zeros([Sz128]);

FigH=4000;
if(param.ShowFig)save
    figure(FigH);close(FigH);
end

if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
StartTime_fnl=now;
param.Verbose=false;save
clear ObjConv Score
for n=1:nfnlCgIters
    [res, CurObj] = fnlCg(res,param);
    ObjConv(n)=CurObj;
    im_res = param.XFM'*res;
    if(param.ShowFig)
        figure(FigH); subplot(1,3,1);
        gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,3,2);
        gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
        subplot(1,3,3);
        plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
        save
    end
    if(n>1)
        dObjP=(ObjConv(n-1)-CurObj)/ObjConv(n-1);
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
        if(dObjP<2e-3)
            disp('Not advancing. Stopping.');
            break;
        end
    else
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
    end
   
end
if(param.ShowFig)
    close(FigH);save
end
disp('ok im_res');
%% ADMM
iter_ops.max_iter = 13;
iter_ops.rho = 111.001;
% iter_ops.objfun = @(a, sv, lam) 0.5*norm_mat(ksp - A_for(a))^2 + lam*sum(sv(:));

llr_ops.lambda = 0.400;
% llr_ops.block_dim = [10, 10];

lsqr_ops.max_iter = 10;
lsqr_ops.tol = 1e-4;

alpha_ref = RefValue;
alpha_ref.data = zeros(Sz128);

AHA=@(x) NUFT_adj_forwS(x,Sens,1,1)
% AHA=@(x) NUFT_TS(x,Sens,fftkernClean);

proxf=@(x,lambda) TV_thresh(x,lambda);
% proxf=@(x,lambda) TV_thresh_rotationInvariant(x,lambda);
history = iter_admm_1prox(alpha_ref, iter_ops, llr_ops, lsqr_ops, AHA, AHA_I, proxf, @admm_callback);

admmRes=alpha_ref.data;

ShowAbsAngle(admmRes)
%% Speed test
nSpeesTest=100;
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
%%
% CGP MB