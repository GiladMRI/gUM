QQ=load('/media/a/DATA/2018_01_25/meas_MID135_BP_gre_Sri_2Dsingleslice_FatSat_FID279_Both.mat');
I=imresizeBySlices(QQ.Combined(:,:,1),Sz1);
B0Real=imresizeBySlices(QQ.B0_Hz(:,:,1),Sz1);

TmpMaskA=imfillholesBySlices(I>5e-5);
TmpMaskE=TmpMaskA;
% se = strel('disk',1);
% TmpMaskE=imerode(TmpMaskA,se);

Z=zeros(Sz1);
MidSz=floor(Sz1/2);
Area=[5 5];
Reg=false(Sz1);
Reg((MidSz(1)-Area(1)):(MidSz(1)+Area(1)),(MidSz(2)-Area(2)):(MidSz(2)+Area(2)))=true;
F=find(Reg);
n=numel(F);
ParamToField=@(x) abs(fft2c(PutValsInIdxs(Z,x(1:n)+1i*x(n+(1:n)),F)));
GetValsInLoc=@(A,B) A(B);
FMsk=find(TmpMaskE);



% CostFunc=@(x) grmss(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+200,FMsk));

W=GetValsInLoc(I,FMsk).';
CostFunc=@(x) W*abs(GetValsInLoc(ParamToField(x),FMsk) - GetValsInLoc(B0Real+200,FMsk));

options = optimset('Display','iter');
B0M=(B0Real+200).*TmpMaskE;
B0M(~isfinite(B0M))=0;
C=ifft2c(B0M);
X0=[real(C(Reg)).' imag(C(Reg)).'];
[BestX, BestCost]=fminsearch(CostFunc,X0,options);

fgmontage(ParamToField(BestX)-200,[-200 300]);colorbar
fgmontage((B0Real).*TmpMaskE,[-200 300]);colorbar

fgmontage(ParamToField(BestX).*TmpMaskA,[0 500]);colorbar

fgmontage(ParamToField(BestX),[0 500]);colorbar

B0RealEx=ParamToField(BestX)-200;
%%
fgmontage(B0Real,[-400 300])
fgmontage(B0RealEx,[-400 300])

fgmontage(B0Real-B0RealEx,[-100 100])