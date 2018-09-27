FN='/media/a/DATA/FC/meas_MID187_mp2rage_iso0_7_iPAT3_FID24573.dat';

AData = mapVBVD(FN);

BData=AData.image();
BDataRef=AData.refscan();
WhichIn=1;
Data=AData.image(:,:,:,:,1,1,1,1,1,WhichIn,1);
Data=squeeze(Data(:,:,:,:,1,1,1,1,1,1));
Data=permute(Data,[1 3 4 2]); % 1 and 3 are full
% So 2 is the subsampled dim

RefData=AData.refscan(:,:,:,:,1,1,1,1,1,1);
RefData=permute(RefData,[1 3 4 2]);
%%
ZData=ifft1cg(Data,1);

ZData1=ZData(50:50:450,:,:,:); 
%%
ZData=ifft1cg(Data,3);

ZData3=ZData(:,:,30:30:220,:);

clear ZData
%%
save('XCCTry_Data.mat');
%%
ZData3P=permute(ZData3,[3 2 1 4]);

ZRefData=ifft1cg(RefData,1);
ZRefData1=ZRefData(50:50:450,:,:,:); 

ZRefData=ifft1cg(RefData,3);
ZRefData3=ZRefData(:,:,30:30:220,:);
ZRefData3P=permute(ZRefData3,[3 2 1 4]);

D1s=exp(1i*(-(1:240).'/AData.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV)*2*pi*AData.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor);
D1sx=repmat(D1s,[1 240]);

D3s=exp(1i*(-(1:240).'/AData.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV)*2*pi*AData.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor);
D3s2=exp(1i*((1:504)/AData.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness    )*2*pi*AData.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra);
D3sx=D3s.*D3s2;

DataSets={ZData1(4:7,:,:,:) ZData3P(2:5,:,:,:)};
RefDataSets={ZRefData1(4:7,:,:,:) ZRefData3P(2:5,:,:,:)};

DShiftC={D1sx D3sx};
disp('ok');
%%
% 1: 4:7, 2: 2:5
WhichDataset=2;
WhichI=3;

for WhichDataset=1:2
    for WhichI=1:size(DataSets{WhichDataset},1)
        disp([WhichDataset WhichI]);
        disp([WhichDataset WhichI]);
        disp([WhichDataset WhichI]);
        disp([WhichDataset WhichI]);
        disp([WhichDataset WhichI]);
% PData=squeeze(ZData(250,:,:,:));

PData=squeeze(DataSets{WhichDataset}(WhichI,:,:,:));
PDataOrig=PData;

% ZRefData=ifft1cg(RefData,1);
% PRefData=squeeze(ZRefData(250,:,:,:,:));
PRefData=squeeze(RefDataSets{WhichDataset}(WhichI,:,:,:,:));

DShift=DShiftC{WhichDataset};

% fgmontage(fft2cg(PRefData))
%%
coils = size(PRefData,3);
kCalib = crop(PRefData,[24,24,coils]);
kSize = [5,5];

[res] = GRAPPA(PData,kCalib, kSize, 0.01);
Rec=ifft2c(res);
% fgmontage(grmss(Rec,3),[0 5e-5]);
%%
Nfull=size(PData);
Ncalib=24;
PadCalib=padBoth(padBoth(kCalib,(Nfull(1)-Ncalib)/2,1),(Nfull(2)-Ncalib)/2,2);
QuickSensM=ifft2c(PadCalib);

calib = bart('ecalib -r 24 -k 6 ',permute(PadCalib,[4 1 2 3]));
ESens = permute(bart('slice 4 0', calib),[2 3 4 1]);

RecE=sum(Rec.*conj(ESens),3);
% ShowAbsAngle(RecE.*DShift);
%% XCC to 3
AccK=3;
ToShift=size(ESens,1)/AccK;
clear ShiftedESens mShiftedESens ShiftedESensX mShiftedESensX ESensX mESensX ShRecE3 Step2 ShXCCFPData
for i=1:AccK
    ShiftedESens(:,:,:,i)=circshift(ESens,ToShift*(i-1),1);
    mShiftedESens(:,:,:,i)=circshift(ESens,-ToShift*(i-1),1);
    ShiftedESensX(:,:,i)=sum(ShiftedESens(:,:,:,i).*conj(ESens),3);
    mShiftedESensX(:,:,i)=sum(mShiftedESens(:,:,:,i).*conj(ESens),3);
    ESensX(:,:,i)=circshift(ShiftedESensX(:,:,i),-ToShift*(i-1),1);
    mESensX(:,:,i)=circshift(mShiftedESensX(:,:,i),ToShift*(i-1),1);
end
% ShowAbsAngle(ESensX,[],'Size',[1 AccK])
CShiftedESens=conj(ShiftedESens);
mCShiftedESens=conj(mShiftedESens);
CESensX=conj(ESensX);

RecE3=sum(Rec.*CShiftedESens,3);
SRecE3=sum(RecE3,3);

for i=1:AccK
    ShRecE3(:,:,i,:)=circshift(RecE3,-ToShift*(i-1),1);
end

RecI=RecE3(:,:,1);
Step1=RecI.*mESensX; % RecE3
Shifts=-ToShift*((1:3)-1);
for i=1:AccK
    Step2(:,:,i)=circshift(Step1(:,:,i),Shifts(i),1);
end
XCCSp=permute(exp(1i*2*pi*(0:2)/3) ,[1 3 2]);
Step3=sum(Step2.*XCCSp,3);
XCCOp= XCCSenseOp(mESensX,Shifts,XCCSp);

FOP=XCCOp*RecI;
AOP=XCCOp'*FOP;

FPadCalib=ifft2c(PadCalib);
XCCFPadCalib=squeeze(sum(FPadCalib.*CShiftedESens,3));
FXCCFPadCalib=fft2c(XCCFPadCalib);
kFXCCFPadCalib = crop(FXCCFPadCalib,[24,24,AccK]);

% ShowAbsAngle(RecE3,[],'Size',[1 AccK]);ShowAbsAngle(XCCFPadCalib,[],'Size',[1 AccK])

FPData=ifft2c(PData);
% FPData=ifft2c(circshift(PData,-1,1));
XCCFPData=squeeze(sum(FPData.*CShiftedESens,3));

for i=1:AccK
    ShXCCFPData(:,:,i)=circshift(XCCFPData(:,:,i),-ToShift*(i-1),1);
end

A1=ShRecE3(:,:,1,1)+ShRecE3(:,:,2,2)*exp(1i*2*pi/3)+ShRecE3(:,:,3,3)*exp(-1i*2*pi/3);
B1=XCCFPData(:,:,1);
A2=ShRecE3(:,:,1,2)+ShRecE3(:,:,2,3)*exp(1i*2*pi/3)+ShRecE3(:,:,3,1)*exp(-1i*2*pi/3);
B2=XCCFPData(:,:,2);
A3=ShRecE3(:,:,1,3)+ShRecE3(:,:,2,1)*exp(1i*2*pi/3)+ShRecE3(:,:,3,2)*exp(-1i*2*pi/3);
B3=XCCFPData(:,:,3);

% fgmontage(ShRecE3(:,:,1,1)+ShRecE3(:,:,2,2)*exp(1i*2*pi/3)+ShRecE3(:,:,3,3)*exp(-1i*2*pi/3),[0 25e-5]);title('From Recx')
% fgmontage(XCCFPData(:,:,1),[0 9e-5]);title('XCCFPData(:,:,1)')
% 
% fgmontage(ShRecE3(:,:,1,2)+ShRecE3(:,:,2,3)*exp(1i*2*pi/3)+ShRecE3(:,:,3,1)*exp(-1i*2*pi/3),[0 25e-5]);title('From Recx2')
% fgmontage(XCCFPData(:,:,2),[0 9e-5]);title('XCCFPData(:,:,2)')
% 
% fgmontage(ShRecE3(:,:,1,3)+ShRecE3(:,:,2,1)*exp(1i*2*pi/3)+ShRecE3(:,:,3,2)*exp(-1i*2*pi/3),[0 25e-5]);title('From Recx3')
% fgmontage(XCCFPData(:,:,3),[0 9e-5]);title('XCCFPData(:,:,3)')
disp('ok');
%%
Sz2=size(RecI);
x = randn(Sz2) + 1j*randn(Sz2);
y = randn(Sz2) + 1j*randn(Sz2);
Ax = XCCOp*x;
Aty = XCCOp'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)))
%%
TVW=1e-5;
XFM=1;
xfmWeight=0;

% XFMStr='Daubechies';
% filterSize=4;
% wavScale=4;
% XFM = Wavelet(XFMStr,filterSize,wavScale);
% xfmWeight = 1e-6;	% Weight for Transform L1 penalty
%%
PData=PDataOrig;
NoiseVal=grmss(PData)*0;
PDataMsk=abs(PData)>0;
PDataVals=PData(PDataMsk);
PData(PDataMsk)=PDataVals+randn(size(PDataVals))*NoiseVal+1j*randn(size(PDataVals))*NoiseVal;
disp('Added noise');
%% XCC
ShowFig=false;
FPData=ifft2c(PData)*3;
DataP=squeeze(sum(FPData.*conj(ESens),3));
AOdd = XCCOp;
RunfnlCgIterationsx;
im_resXCC=im_res;
% ShowAbsAngle(im_resXCC.*DShift);
%% Standard sens, all maps
DataP=PData;
Msk=zeros(size(PData));
Msk(3:3:end,:,:)=1;
AOdd = p2DFTC(Msk,size(PData),ESens,2);
RunfnlCgIterationsx;
im_resAllMaps=im_res;
% fgmontage(im_resAllMaps.*DShift);

[res_AllMaps] = GRAPPA(PData,kCalib, kSize, 0.01);
Rec_AllMaps=ifft2c(res_AllMaps);
RecGRAPPA_AllMaps=sum(Rec_AllMaps.*conj(ESens),3);
disp('OK all maps');
%% SCC
nScc=11;
sccmtx = calcSCCMtx(kCalib);
SCCDATA = CC(PData,sccmtx(:,1:nScc));
ESensSCC=CC(ESens,sccmtx(:,1:nScc));
DataP=SCCDATA;
Msk=zeros(size(DataP));
Msk(3:3:end,:,:)=1;
AOdd = p2DFTC(Msk,size(DataP),ESensSCC,2);
RunfnlCgIterationsx;
im_resSCC=im_res;
% ShowAbsAngle(im_resSCC.*DShift);

kCalib_SCC=CC(kCalib,sccmtx(:,1:nScc));
[res_SCC] = GRAPPA(SCCDATA,kCalib_SCC, kSize, 0.01);
Rec_SCC=ifft2c(res_SCC);
RecGRAPPA_SCC=sum(Rec_SCC.*conj(ESensSCC),3);
disp('OK SCC');
%% GCC
nGcc=5;
dim=2;
slwin = 5; % sliding window length for GCC
% Pad2Calib=padBoth(kCalib,(Nfull-Ncalib)/2,2);
gccmtx = calcGCCMtx(PRefData,dim,slwin);
% gccmtx = calcGCCMtx(Pad2Calib,dim,slwin);
gccmtx_aligned = alignCCMtx(gccmtx(:,1:nGcc,:));
GCCDATA = CC(PData,gccmtx_aligned, dim);
ESensGCC=ifft1cg(CC(fft1cg(ESens,dim),gccmtx_aligned, dim),dim);
% ESensGCC=ESensGCC.*grmss(ESens,3)*sqrt(size(ESens,3));
DataP=GCCDATA;
Msk=zeros(size(DataP));
Msk(3:3:end,:,:)=1;
AOdd = p2DFTC(Msk,size(DataP),ESensGCC,2);
RunfnlCgIterationsx;
im_resGCC=im_res;
% ShowAbsAngle(im_resGCC.*DShift);

kCalib_GCC=CC(PRefData,gccmtx_aligned,dim);
[res_GCC] = GRAPPA(GCCDATA,kCalib_GCC, kSize, 0.01);
Rec_GCC=ifft2c(res_GCC);
RecGRAPPA_GCC=sum(Rec_GCC.*conj(ESensGCC),3);

disp('OK GCC');
%%
AllRes=cat(3,im_resXCC,im_resAllMaps,im_resSCC,im_resGCC,RecGRAPPA_AllMaps,RecGRAPPA_SCC,RecGRAPPA_GCC).*DShift;

AllResC{WhichDataset}(:,:,:,WhichI)=AllRes;

mESensXC{WhichDataset}(:,:,:,WhichI)=mESensX;
ShiftedESensXC{WhichDataset}(:,:,:,WhichI)=ShiftedESensX;
ESensXC{WhichDataset}(:,:,:,WhichI)=ESensX;
Step2C{WhichDataset}(:,:,:,WhichI)=Step2;
Step3C{WhichDataset}(:,:,:,WhichI)=Step3;

end
end
%%
ShowAbsAngle(AllResC{1})
gprint(get(gcf,'Number'),'ACC1');
close(gcf)
ShowAbsAngle(AllResC{2})
gprint(get(gcf,'Number'),'ACC2');
close(gcf)
%%
WhichI=1;
ShowAbsAngle(Step3C{1}(:,:,WhichI).*DShiftC{1});
gprint(get(gcf,'Number'),['Step3_1_' num2str(WhichI)]);close(gcf)
ShowAbsAngle(Step2C{1}(:,:,:,WhichI).*DShiftC{1},1,[],'Size',[1 3]);
gprint(get(gcf,'Number'),['Step2_1_' num2str(WhichI)]);close(gcf)
ShowAbsAngle(ShiftedESensXC{1}(:,:,:,WhichI),1,[],'Size',[1 3]);
gprint(get(gcf,'Number'),['ShiftedESensXC_1_' num2str(WhichI)]);close(gcf)
%%
im_resXCC0=im_resXCC.*D1s;
im_resAllMaps0=im_resAllMaps.*D1s;
im_resSCC0=im_resSCC.*D1s;
im_resGCC0=im_resGCC.*D1s;
RecGRAPPA_AllMaps0=RecGRAPPA_AllMaps.*D1s;
RecGRAPPA_SCC0=RecGRAPPA_SCC.*D1s;
RecGRAPPA_GCC0=RecGRAPPA_GCC.*D1s;
%%
D(1)=grmss(im_resXCC-im_resXCC0);
D(2)=grmss(im_resAllMaps-im_resAllMaps0);
D(3)=grmss(im_resSCC-im_resSCC0);
D(4)=grmss(im_resGCC-im_resGCC0);
D(5)=grmss(RecGRAPPA_AllMaps-RecGRAPPA_AllMaps0);
D(6)=grmss(RecGRAPPA_SCC-RecGRAPPA_SCC0);
D(7)=grmss(RecGRAPPA_GCC-RecGRAPPA_GCC0);
D
%% GRAPPA direction - work on non-ZF data
WhichDataset=2;
WhichI=3;


PData=squeeze(DataSets{WhichDataset}(WhichI,:,:,:));
PDataOrig=PData;

PRefData=squeeze(RefDataSets{WhichDataset}(WhichI,:,:,:,:));

DShift=DShiftC{WhichDataset};


% ZData1x=squeeze(ZData1(4,:,:,:));
PDataU=PData(AccK:AccK:end,:,:);

% ZRefData=ifft1cg(RefData,1);
% ZRefData1=ZRefData(50:50:450,:,:,:); 

% ZRefData1x=squeeze(ZRefData1(4,:,:,:));

% PData=ZData1x;
% PDataOrig=PData;

% PRefData=ZRefData1x;

coils = size(PRefData,3);
kCalib = crop(PRefData,[24,24,coils]);
kSize = [5,5];

[res] = GRAPPA(PData,kCalib, kSize, 0.01);
Rec=ifft2c(res);
% fgmontage(grmss(Rec,3),[0 5e-5]);
%
Nfull=size(PData);
Ncalib=24;
PadCalib=padBoth(padBoth(kCalib,(Nfull(1)-Ncalib)/2,1),(Nfull(2)-Ncalib)/2,2);
QuickSensM=ifft2c(PadCalib);

calib = bart('ecalib -r 24 -k 6 ',permute(PadCalib,[4 1 2 3]));
ESens = permute(bart('slice 4 0', calib),[2 3 4 1]);

% calibx = bart('ecalib -r 24 -k 6 ',permute(kCalib,[4 1 2 3]));
% ESensx = permute(bart('slice 4 0', calibx),[2 3 4 1]);


RecE=sum(Rec.*conj(ESens),3);

AccK=3;
ToShift=size(ESens,1)/AccK;

ESensK=PartitionDim(ESens,1,AccK);

ESensKK=repmat(ESensK,[AccK 1 1 1]);
ESensKBig=imresizeBySlices(ESensK,Nfull);

clear ShiftedESens mShiftedESens ShiftedESensX mShiftedESensX ESensX mESensX ShRecE3 Step2 ShXCCFPData
for i=1:AccK
    ShiftedESens(:,:,:,i)=circshift(ESens,ToShift*(i-1),1);
    mShiftedESens(:,:,:,i)=circshift(ESens,-ToShift*(i-1),1);
    ShiftedESensX(:,:,i)=sum(ShiftedESens(:,:,:,i).*conj(ESens),3);
    mShiftedESensX(:,:,i)=sum(mShiftedESens(:,:,:,i).*conj(ESens),3);
    ESensX(:,:,i)=circshift(ShiftedESensX(:,:,i),-ToShift*(i-1),1);
    mESensX(:,:,i)=circshift(mShiftedESensX(:,:,i),ToShift*(i-1),1);
end
% ShowAbsAngle(ESensX,[],'Size',[1 AccK])
CShiftedESens=conj(ShiftedESens);
mCShiftedESens=conj(mShiftedESens);
CESensX=conj(ESensX);

disp('ok 1');
%%
% IkCalib=ifft2c(kCalib);
RIdxs=linspaceWithHalfStep(0,1,ToShift);
for i=1:ToShift
    RMat(:,:,i)=(grotx(-RIdxs(i)*pi/2)*grotz(-RIdxs(i)*pi/2));
end
ESensKA=permute(sum(ESensK.*permute(RMat,[3 4 5 1 2]),4),[1 2 3 5 4]);
fgmontage(AddSepLines(grmss(ESensKA,3)))
% ESensKA=sum(ESensK.*permute(RMat,[3 4 5 1 2]),5);
ESensKKA=repmat(ESensKA,[AccK 1 1 1]);
%%
clear ESensKM ESensKM
for i=1:AccK
    for j=1:AccK
        ESensKM(:,:,i,j)=sum(ESensK(:,:,:,i).*conj(ESensK(:,:,:,j)),3);
        ESensKMA(:,:,i,j)=sum(ESensKA(:,:,:,i).*conj(ESensKA(:,:,:,j)),3);
    end
end
ESensKMC=CombineDims(ESensKM,[3 1]);
ESensKMCA=CombineDims(ESensKMA,[3 1]);

FPDataU=ifft2cg(PDataU);

FACC_Step1=squeeze(sum(FPDataU.*conj(ESensK),3));
FACC_Step2=fft2cg(FACC_Step1);

FACC_Step1A=squeeze(sum(FPDataU.*conj(ESensKA),3));
FACC_Step2A=fft2cg(FACC_Step1A);
%%
FACC_ZF=zeros(Nfull(1),Nfull(2),AccK);
FACC_ZF(AccK:AccK:end,:,:)=FACC_Step2A;

QuickSensM_ACC=squeeze(sum(QuickSensM.*conj(ESensKKA),3));

% PadCalibS=padBoth(padBoth(kCalib,(ToShift-Ncalib)/2,1),(Nfull(2)-Ncalib)/2,2);
% QuickSensMS=ifft2c(PadCalibS);

% QuickSensMS_ACC=squeeze(sum(QuickSensMS.*conj(ESensK),3));

kCalibACCBig=fft2cg(QuickSensM_ACC);
kCalibACC=crop(kCalibACCBig,[24,Nfull(2),AccK]);

[resACC] = GRAPPA(FACC_ZF,kCalibACC, kSize, 0.01);
RecACC=ifft2c(resACC);

RecACCK=PartitionDim(RecACC,1,3);

RecACCR=cat(1,RecACCK(:,:,1,1),RecACCK(:,:,2,2),RecACCK(:,:,3,3));
% ShowAbsAngle(AddSepLines(RecACC))

% ShowAbsAngle(RecACC,1,[],'Size',[1 3])
% ShowAbsAngle(AddSepLines(ESensKM),1,[0 1.2])

% ShowAbsAngle(AddSepLines(RecACC.*DShift),1,[0 2e-4],'Size',[1 3]);subplot(121);title('Reconstructed channels');

RecACCR2=sum(RecACC.*conj(ESensKMCA),3);
% RecACCR2=sum(RecACC./conj(ESensKMCA),3);

% ShowAbsAngle(RecACCR2.*DShift)

fgmontage(grmss(RecACC,3))
fgmontage(RecE.*DShift);title('GRAPPA+GCC');
%% Try on RO
% ESensGCC_RO1=PartitionDim(ESensGCC,1,AccK);
% ESensGCC_RO=CombineDims(ESensGCC_RO1,[4 2]);
% FESensGCC_RO=fft2cg(ESensGCC_RO);
% kCalibGCC_RO=crop(FESensGCC_RO,[8,72,nGcc]);
% GCCDATA_RO=zeros(ToShift,Nfull(2)*AccK,nGcc);
% GCCDATA_RO(:,AccK:AccK:end,:)=GCCDATA(AccK:AccK:end,:,:);
% 
% kSize_RO=[3 15];
% [res_GCC_RO] = GRAPPA(GCCDATA_RO,kCalibGCC_RO, kSize_RO, 0.01);
% Rec_GCC_RO=ifft2c(res_GCC_RO);
% 
% Rec_GCC_RO1=PartitionDim(Rec_GCC_RO,2,AccK);
% Rec_GCC_ROC=CombineDims(Rec_GCC_RO1,[4 1]);
% 
% Rec_GCC_ROCX=sum(Rec_GCC_ROC.*conj(ESensGCC),3);

% ShowAbsAngle(Rec_GCC_ROCX)
% ShowAbsAngle(RecGRAPPA_GCC)
%% Now with ACC on RO:
% ESensKM_RO=CombineDims(ESensKM,[3 2]);
% FESensKM_RO=fft2cg(ESensKM_RO);
% kCalibACC_RO=crop(FESensKM_RO,[8,72,AccK]);

QuickSensM_RO1=PartitionDim(QuickSensM,1,AccK);
QuickSensM_RO=CombineDims(QuickSensM_RO1,[4 2]);
% ESensK=PartitionDim(ESens,1,AccK);
% ESensK_RO=repmat(ESensK,[1 AccK 1]);

ShiftedESens_RO1=PartitionDim(mShiftedESens,1,AccK);
ShiftedESens_RO=CombineDims(ShiftedESens_RO1,[4 2]);

% QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO.*conj(ESensK_RO),3));
QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO.*conj(ShiftedESens_RO),3));
FQuickSensM_RO_ACC=fft2cg(QuickSensM_RO_ACC);
kCalibACC_RO=crop(FQuickSensM_RO_ACC,[8,72,AccK]);

% QuickSensMP_ACC=squeeze(sum(QuickSensMS.*conj(ESensK),3));



% QuickSensM_RO1=PartitionDim(QuickSensM,1,AccK);
% QuickSensM_RO=CombineDims(QuickSensM_RO1,[4 2]);
% 
% QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO.*conj(ESens_RO),3));
% 
% FQuickSensM_RO_ACC=fft2cg(QuickSensM_RO_ACC);


FACC_ZF_RO=zeros(ToShift,Nfull(2)*AccK,AccK);
FACC_ZF_RO(:,AccK:AccK:end,:)=FACC_Step2;

[res_ACC_RO] = GRAPPA(FACC_ZF_RO,kCalibACC_RO, kSize_RO, 0.01);
Rec_ACC_RO=ifft2c(res_ACC_RO);

% Rec_ACC_ROC=sum(Rec_ACC_RO.*(ESensKM_RO),3);
Rec_ACC_ROC=sum(Rec_ACC_RO.*conj(ESensKM_RO),3);
Rec_ACC_ROC1=PartitionDim(Rec_ACC_ROC,2,AccK);
Rec_ACC_ROCx=CombineDims(Rec_ACC_ROC1,[3 1]);
ShowAbsAngle(Rec_ACC_ROCx.*DShift);title('GRAPPA+ACC');
ShowAbsAngle(RecE.*DShift);title('GRAPPA+GCC');
%%
%% GCC - to work faster
nGcc=5;
dim=2;
slwin = 5; % sliding window length for GCC
gccmtx = calcGCCMtx(PRefData,dim,slwin);
gccmtx_aligned = alignCCMtx(gccmtx(:,1:nGcc,:));
GCCDATA = CC(PData,gccmtx_aligned, dim);
ESensGCC=ifft1cg(CC(fft1cg(ESens,dim),gccmtx_aligned, dim),dim);

kCalib_GCC=CC(PRefData,gccmtx_aligned,dim);
[res_GCC] = GRAPPA(GCCDATA,kCalib_GCC, kSize, 0.01);
Rec_GCC=ifft2c(res_GCC);
RecGRAPPA_GCC=sum(Rec_GCC.*conj(ESensGCC),3).*DShift;

disp('OK GCC');
%% Now with ACC on RO:
PadCalib_GCC=padBoth(padBoth(kCalib_GCC,(Nfull(1)-Ncalib)/2,1),0,2);
QuickSensM_GCC=ifft2c(PadCalib_GCC);

QuickSensM_RO1_GCC=PartitionDim(QuickSensM_GCC,1,AccK);
QuickSensM_RO_GCC=CombineDims(QuickSensM_RO1_GCC,[4 2]);

FQuickSensM_RO_GCC=fft2cg(QuickSensM_RO_GCC);
% kCalibGCC_RO=crop(FQuickSensM_RO_GCC,[8,72,nGcc]);
kCalibGCC_RO=crop(FQuickSensM_RO_GCC,[8,Nfull(2)*AccK,nGcc]);
kCalibGCC_RO=crop(FQuickSensM_RO_GCC,[24,Nfull(2)*AccK,nGcc]);


% QuickSensT_GCC=ifft2c(kCalib_GCC);
% QuickSensT_RO1_GCC=PartitionDim(QuickSensT_GCC,1,AccK);
% QuickSensT_RO_GCC=CombineDims(QuickSensT_RO1_GCC,[4 2]);
% 
% FQuickSensT_RO_GCC=fft2cg(QuickSensT_RO_GCC);
% kCalibGCC_RO=crop(FQuickSensT_RO_GCC,[8,Nfull(2)*AccK,nGcc]);


ESensK_GCC=PartitionDim(ESensGCC,1,AccK);
ESensK_GCC_RO=repmat(ESensK_GCC,[1 AccK 1]);

% ShiftedESens_RO1=PartitionDim(mShiftedESens,1,AccK);
% ShiftedESens_RO=CombineDims(ShiftedESens_RO1,[4 2]);

QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO_GCC.*conj(ESensK_GCC_RO),3));
% QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO.*conj(ShiftedESens_RO),3));
FQuickSensM_RO_ACC=fft2cg(QuickSensM_RO_ACC);
kCalibACC_RO=crop(FQuickSensM_RO_ACC,[8,Nfull(2)*AccK,AccK]);
kCalibACC_RO=crop(FQuickSensM_RO_ACC,[24,Nfull(2)*AccK,AccK]);

% QuickSensMP_ACC=squeeze(sum(QuickSensMS.*conj(ESensK),3));



% QuickSensM_RO1=PartitionDim(QuickSensM,1,AccK);
% QuickSensM_RO=CombineDims(QuickSensM_RO1,[4 2]);
% 
% QuickSensM_RO_ACC=squeeze(sum(QuickSensM_RO.*conj(ESens_RO),3));
% 
% FQuickSensM_RO_ACC=fft2cg(QuickSensM_RO_ACC);


% FACC_ZF_RO=zeros(ToShift,Nfull(2)*AccK,AccK);
% FACC_ZF_RO(:,AccK:AccK:end,:)=FACC_Step2;
% 
% [res_ACC_RO] = GRAPPA(FACC_ZF_RO,kCalibACC_RO, kSize_RO, 0.01);

FGCC_DataU=GCCDATA(AccK:AccK:end,:,:);
FGCC_ZF_RO=zeros(ToShift,Nfull(2)*AccK,nGcc);
FGCC_ZF_RO(:,(AccK-1):AccK:end,:)=FGCC_DataU;

[res_GCC_RO] = GRAPPA(FGCC_ZF_RO,kCalibGCC_RO, [3 15], 0.01);

% res_GCC_RO_F1=ifft1cg(res_GCC_RO,2);
% res_GCC_RO_F1P=PartitionDim(res_GCC_RO_F1,2,AccK);
% res_GCC_RO_F1PC=CombineDims(res_GCC_RO_F1P,[1 4]);
% res_GCC_RO_F2=ifft1cg(res_GCC_RO_F1PC,1);

Rec_GCC_RO=ifft2c(res_GCC_RO);

Rec_GCC_ROC1=PartitionDim(Rec_GCC_RO,2,AccK);
Rec_GCC_ROCx=CombineDims(Rec_GCC_ROC1,[4 1]);

ShowAbsAngle(Rec_GCC)
ShowAbsAngle(Rec_GCC_ROCx)
%%
% Rec_ACC_ROC=sum(Rec_ACC_RO.*(ESensKM_RO),3);
Rec_ACC_ROC=sum(Rec_ACC_RO.*conj(ESensKM_RO),3);
Rec_ACC_ROC1=PartitionDim(Rec_ACC_ROC,2,AccK);
Rec_ACC_ROCx=CombineDims(Rec_ACC_ROC1,[3 1]);
ShowAbsAngle(Rec_ACC_ROCx.*DShift);title('GRAPPA+ACC');
ShowAbsAngle(RecE.*DShift);title('GRAPPA+GCC');
%%
% XCCFPDataR=XCCFPData*0;
% XCCFPDataC=XCCFPDataR;
% CentralIdxs=(1:ToShift)+(ToShift/2)*(AccK-1);
% XCCFPDataCn=XCCFPData(CentralIdxs,:,:);
% FXCCFPDataCn=fft2c(XCCFPDataCn);
% FXCCFPDataCnF=XCCFPDataR;
% FXCCFPDataCnF(3:3:end,:,:)=FXCCFPDataCn;
% XCCFPDataC(CentralIdxs,:,:)=XCCFPData(CentralIdxs,:,:);
% FXCCFPData=fft2c(XCCFPData);
% for i=1:AccK
%     XCCFPDataR=XCCFPDataR+circshift(XCCFPDataC,ToShift*(i-1),1);
%     ShiftedXCCFPData(:,:,i)=circshift(XCCFPData(:,:,i),-ToShift*(i-1),1);
%     ShXCCFPData(:,:,:,i)=circshift(XCCFPData,-ToShift*(i-1),1);
% end
% FXCCFPDataR=fft2c(XCCFPDataR);
% % FXCCFPDataR=circshift(FXCCFPDataR,2,1);
% % [resX] = GRAPPA(FXCCFPDataR,kFXCCFPadCalib, kSize, 0.01);
% [resX] = GRAPPA(FXCCFPDataCnF,kFXCCFPadCalib, kSize, 0.01);
% RecQQ=ifft2c(resX);
% 
% ShowAbsAngle(RecQQ,[],'Size',[1 AccK])
% % ShowAbsAngle(XCCFPData,[],'Size',[1 AccK])
% % ShowAbsAngle(ShiftedXCCFPData,[],'Size',[1 AccK])
% %%
% 
% kCalibXCC = crop(FXCCFPadCalib,[24,24,AccK]);
% 
% [resXCC] = GRAPPA(FXCCFPData,kCalibXCC, kSize, 0.01);
% RecXCC=ifft2c(resXCC);
% fgmontage(grmss(RecXCC,3),[0 5e-5]);
% %%
% % clear Data ZData
% save('XCCTry1.mat')
% %%
% % ADataI=AData.image();
% ADataIx=AData.image(:,:,:,:,:,3,:,:,:,:,:,:,:);
% 
% %%
% [x,y] = meshgrid(linspace(0,1,128));
% % Generate fake Sensitivity maps
% sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
% % generate 4 coil phantom
% BaseI=phantom(128);
% imgs = repmat(BaseI,[1,1,4]).*sMaps;
% DATA = fft2c(imgs);
% % crop 20x20 window from the center of k-space for calibration
% kCalib = crop(DATA,[20,20,4]);
% 
% %calibrate a kernel
% kSize = [5,5];
% coils = 4;
% 
% % undersample by a factor of 2
% % DATA(1:2:end,2:2:end,:) = 0;
% % DATA(2:2:end,1:2:end,:) = 0;
% DATA(2:2:end,:,:) = 0; % g
% 
% %reconstruct:
% [res] = GRAPPA(DATA,kCalib, kSize, 0.01);
% figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[0 1]);
% title('full,  zero-fill,   result')
%%
N1=11;
a=rand(1,N1*3);
af=fft1cg(a,2);
bf=af*0;
bf(2:3:end)=af(2:3:end);
ibf=ifft1cg(bf,2);

aa=a+circshift(a,N1,2)+circshift(a,N1*2,2);

abs([aa; ibf*3]);
%%
I=phantom(max(Nfull));
I=imresize(I,Nfull(1:2));
IWithS=I.*ESens;
FIWithS=fft2cg(IWithS);

