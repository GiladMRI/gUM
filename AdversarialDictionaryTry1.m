%% get data
SzI=[128 128];
Ni=500;
A=load('/media/a/H1/First3kIm128x128MagSinglex.mat');
A=A.First3kIm128x128MagSingle;
A=permute(A,[2 3 1]);
A=A(:,:,getKrandomSamples(3000,Ni));
% AR=imresizeBySlices(A,SzI);
disp('ok');
%% add phase
for i=1:Ni
    if(mod(i,100)), disp(i), end
    Out=GenerateRandomSinPhase(SzI(1),5,.1);
    A(:,:,i)=A(:,:,i).*single(Out);
end
disp('Added phase');
%%
load brain_8ch
nChannels=8;
DATA=crop(DATA,[SzI nChannels]);
I=ifft2c(DATA);
Sens=RunESPIRiTForSensMaps(I);
SMsk=grmss(Sens,3)>0.01;
Im=CalcSENSE1f(I,Sens);
ImWithSens=Im.*Sens;
FData=fft2cg(ImWithSens);
imSize=gsize(Sens,1:2);
DataSize=size(Sens);
%%
Acc=3;
% Msk=bart(['poisson -Y ' num2str(SzI(1)) ' -Z ' num2str(SzI(1)) ' -y ' num2str(Acc) ' -z ' num2str(Acc) ' -v -e -C 8']);
% Msk=squeeze(Msk);

Msk=zeros(SzI);
Msk(1:6:end,:)=1;
Msk(64+(-2:2),:)=1;

Op = p2DFTC(Msk,DataSize,Sens,2);

x = randn(imSize) + 1j*randn(imSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))

TVW=1e-33;
XFM=1;
xfmWeight=0;
Sz2=imSize;
%%
AllRes=NaN([SzI Ni]);
ASMsk=A.*SMsk;
%%
ShowFig=false;
for i=1:Ni
    disp(['Image #' num2str(i)]);
    Im=A(:,:,i).*SMsk;
    ImWithSens=Im.*Sens;
    FData=fft2cg(ImWithSens);
    MData=Msk.*FData;
    %
    DataP=MData;
    AOdd = Op;
    RunfnlCgIterationsx;
    AllRes(:,:,i)=im_res;
end
%%
fgmontage(im_res);
fgmontage(Im);
%%
nImToUse=30;
ResToUse=AllRes(:,:,1:nImToUse);
RefImToUse=ASMsk(:,:,1:nImToUse);
AllResN=ResToUse.*grms(RefImToUse,1:2)./grms(ResToUse,1:2);
blockSizeR = 5; % Rows in block.
blockSizeC = blockSizeR; % Columns in block.
blockSize=blockSizeR*blockSizeC;
clear RefPatches RecPatches
for x=1:(SzI(1)+1-blockSizeR)
    for y=1:(SzI(2)+1-blockSizeC)
        RefPatches(:,:,:,x,y)=RefImToUse(x:(x+blockSizeR-1),y:(y+blockSizeC-1),:);
        RecPatches(:,:,:,x,y)=AllResN(x:(x+blockSizeR-1),y:(y+blockSizeC-1),:);
    end
end
RefPatchesM=CombineDims(CombineDims(RefPatches,[5,3]),[4,3]);
RecPatchesM=CombineDims(CombineDims(RecPatches,[5,3]),[4,3]);
disp('ok');
%%
nPatchToUse=5000;
RefPatchesMV=reshape(RefPatchesM(:,:,getKrandomSamples(size(RefPatchesM,3),nPatchToUse)),blockSize,[]);
RecPatchesMV=reshape(RecPatchesM(:,:,getKrandomSamples(size(RefPatchesM,3),nPatchToUse)),blockSize,[]);
%%
[COEFFRef, SCORERef, LATENTRef] = pca(RefPatchesMV.');
[COEFFRec, SCORERec, LATENTRec] = pca(RecPatchesMV.');
% Each column of COEFF contains coefficients for one principal component
COEFFRefM=reshape(COEFFRef,blockSizeR,blockSizeC,[]);
COEFFRecM=reshape(COEFFRec,blockSizeR,blockSizeC,[]);
% ShowAbsAngle(COEFFRefM)
% ShowAbsAngle(COEFFRecM)
figure;plot(LATENTRef,'k');hold on;plot(LATENTRec,'b')
%%
% x = quadprog(H,f) returns a vector x that minimizes 1/2*x'*H*x + f'*x. H must be positive definite for the problem to have a finite minimum.
% Multiply the columns of coeffre by the sigular value
% we want v such that
% (v*COEFFRec)^2 is big
% (v*COEFFRef)^2 is small
% v is unit norm
% so we want to maximize
% COEFFRec*COEFFRec'-COEFFRef*COEFFRef' - lambda(eye)
% so H is
clc
% lambda=10;
WRef=.1;
LATENTRecX=LATENTRec;
LATENTRefX=LATENTRef;
LATENTRecX(23:end)=0;
LATENTRefX(23:end)=0;
COEFFRecX=COEFFRec.*(LATENTRecX.');
COEFFRefX=COEFFRef.*(LATENTRefX.');
H= (COEFFRecX*COEFFRecX') - WRef*(COEFFRefX*COEFFRefX');

% x0=rand(1,blockSize*2);
try
    x0=[real(COEFFRefX(:,40)); -imag(COEFFRefX(:,40))].';
catch
    x0=[real(COEFFRefX(:,end)); -imag(COEFFRefX(:,end))].';
end
x0=double(x0);
RIToC=@(x) x(1:blockSize)+1i*x(blockSize+(1:blockSize));
NormX=@(x) x/(grmss(x)+0.00001);
ByH=@(x) -real(x*H*x');
% ByH2=@(x) x*H*x';
VarW=.1;
CostFunc=@(x) ByH(NormX(RIToC(x)))+VarW*var(NormX(RIToC(x)));
% options = optimset('MaxFunEvals',1000);
options = optimset();
[BestX, BestVal]=fminsearch(CostFunc,x0);
CostFuncC=@(x) ByH(NormX(x))+VarW*var(NormX(x));
BestVal
BestXc=RIToC(BestX);
BestXc=BestXc./grmss(BestXc);
BestXcM=reshape(BestXc,[blockSizeR,blockSizeC]);
BestXc1=BestXc;
ShowAbsAngle(BestXcM)
%%
H2= H-1*(BestXc'*BestXc);
ByH2=@(x) -real(x*H2*x');
% ByH2=@(x) x*H*x';
VarW=.1;
CostFunc2=@(x) ByH2(NormX(RIToC(x)))+VarW*var(NormX(RIToC(x)));
% options = optimset('MaxFunEvals',1000);
options = optimset();
x0=BestX;
try
    x0=[real(COEFFRefX(:,40)); -imag(COEFFRefX(:,40))].';
catch
    x0=[real(COEFFRefX(:,end)); -imag(COEFFRefX(:,end))].';
end
x0=double(x0);
[BestX2, BestVal2]=fminsearch(CostFunc2,x0);
CostFuncC2=@(x) ByH2(NormX(x))+VarW*var(NormX(x));
BestVal2
BestXc2=RIToC(BestX2);
BestXc2=BestXc2./grmss(BestXc2);
BestXcM2=reshape(BestXc2,[blockSizeR,blockSizeC]);
BestXc2=BestXc;
ShowAbsAngle(BestXcM2)
%%
% K=zeros(11,1);
% K(6,1)=1;
% XOP=TVXOP(K);
% XOP=TVXOP({BestXcM});
% Fltrs={[1 -1 0],[1 -1 0].',100*BestXcM,1*BestXcM2};
% Fltrs={[1 -1 0],[1 -1 0].',0.1*gflip((BestXcM),1:2)}; % 0.80286
Fltrs={[1 -1 0],[1 -1 0].',1*gflip((BestXcM),1:2)}; % 0.80286
Fltrs={COEFFRefM(:,:,2) COEFFRefM(:,:,3) COEFFRefM(:,:,4) COEFFRefM(:,:,5)};
% Fltrs={[1 -1 0],[1 -1 0].',0.1*gflip((BestXcM),1:2),0.1*gflip((BestXcM2),1:2)};
% Fltrs={[1 -1 0],[1 -1 0].',1*gflip((BestXcM),1:2),1*gflip((BestXcM2),1:2)};
% Fltrs={[1 -1 0],[1 -1 0].',.1*gflip((BestXcM),1:2),.1*gflip((BestXcM2),1:2)};
% Fltrs={[1 -1 0],[1 -1 0].',1*gflip((BestXcM3),1:2),1*gflip((BestXcM),1:2)};
% Fltrs={[1 -1 0],[1 -1 0].'}; % 0.80186
nFltrs=numel(Fltrs);
XOP=TVXOP(Fltrs);
Sz2=size(Im);
x = randn(Sz2) + 1j*randn(Sz2);
y = randn([Sz2 nFltrs]) + 1j*randn([Sz2 nFltrs]);
Ax = XOP*x;
Aty = XOP'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
%%
TVW=1e-7;
XFM=1;
xfmWeight=0;
Sz2=imSize;

    disp(['Image #' num2str(i)]);
    Im=A(:,:,i).*SMsk;
    ImWithSens=Im.*Sens;
    FData=fft2cg(ImWithSens);
    MData=Msk.*FData;
    %
    DataP=MData;
    AOdd = Op;
%     RunfnlCgIterationsx;
%
param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,...
    'TV',XOP,'xfmWeight',xfmWeight,'ShowFig',true),init);
param.data =     DataP;
nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
FigH=4000;
if(param.ShowFig)
    figure(FigH);close(FigH);
end
StartTime_fnl=now;
param.Verbose=false;
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
    end
%     t=toc;
    if(n>1)
        dObjNom=ObjConv(n-1)-CurObj;
        dObjP=dObjNom/ObjConv(n-1);
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
        if(dObjP<2e-3 || dObjNom<1e-16)
            disp('Not advancing. Stopping.');
            break;
        end
    else
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
    end
end
if(param.ShowFig)
    close(FigH)
end
disp('ok im_res');
%     AllRes(:,:,i)=im_res;
fgmontage(cat(3,im_res,AllRes(:,:,i),ASMsk(:,:,i)),'Size',[1 3])
title(num2str([ssim(double(abs(im_res)),double(abs(ASMsk(:,:,i)))) ssim(abs(AllRes(:,:,i)),double(abs(ASMsk(:,:,i))))]));
%%
% x0=rand(1,blockSize*2);
x0=[real(COEFFRefX(:,40)); -imag(COEFFRefX(:,40))].';
x0=double(x0);
RIToC=@(x) x(1:blockSize)+1i*x(blockSize+(1:blockSize));
NormX=@(x) x/(grmss(x)+0.00001);
ByH=@(x) -real(x*H*x');
% ByH2=@(x) x*H*x';
VarW=20;
CostFunc=@(x) ByH(NormX(RIToC(x)))+VarW*var(NormX(RIToC(x)));
% options = optimset('MaxFunEvals',1000);
options = optimset();
[BestX, BestVal]=fminsearch(CostFunc,x0);
CostFuncC=@(x) ByH(NormX(x))+VarW*var(NormX(x));
BestVal
BestXc=RIToC(BestX);
BestXc=BestXc./grmss(BestXc);
BestXcM=reshape(BestXc,[blockSizeR,blockSizeC]);
BestXc1=BestXc;

%%
% Figure out the size of each block in rows. 
% Most will be blockSizeR but there may be a remainder amount of less than that.
wholeBlockRows = floor(rows / blockSizeR);
blockVectorR = [blockSizeR * ones(1, wholeBlockRows), rem(rows, blockSizeR)];
% Figure out the size of each block in columns. 
wholeBlockCols = floor(columns / blockSizeC);
blockVectorC = [blockSizeC * ones(1, wholeBlockCols), rem(columns, blockSizeC)];

% Create the cell array, ca.  
% Each cell (except for the remainder cells at the end of the image)
% in the array contains a blockSizeR by blockSizeC by 3 color array.
% This line is where the image is actually divided up into blocks.
if numberOfColorBands > 1
	% It's a color image.
	ca = mat2cell(rgbImage, blockVectorR, blockVectorC, numberOfColorBands);
else
	ca = mat2cell(rgbImage, blockVectorR, blockVectorC);
end

%%
% 
% 
% 
% Lam=1e-7;
% Cmd=['pics -m -R W:7:0:' num2str(Lam) ' -p'];
% 
% FDataP=permute(Msk.*FData,[1 2 4 3]);
% Rec=bart(Cmd,Msk,FDataP,SensP);
% RecN=Rec*grmss(Im)/grmss(Rec);
% %%
% Lams=10.^([-8,-7,-6,-5,-4,-3.5,-3,-2.5,-2,-1.5,-1]);
% %%
% RecNM=zeros(200,200,numel(Accs),numel(Lams));