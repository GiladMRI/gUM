%% demo for SMS PSF processing
% Loop by slice since the dataset could be very large
clear;
close all;
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
FNBase='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
FNBase='meas_MID04684_FID21646_ep2d_ge_EPTI_1p3_5shot_calib';

directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
% FNBase='meas_MID00870_FID32105_ep2d_ge_EPTI_1p9_1shot_calib';
% FNBase='meas_MID00874_FID32109_ep2d_ge_EPTI_1p9_3and5shot_calib';
FNBase='meas_MID00894_FID32129_ep2d_ge_EPTI_1p9_1shot_calib_Cor';
% FNBase='meas_MID00905_FID32140_ep2d_ge_EPTI_1p9_3and5shot_calib_Cor';

directory_rawdata = '/autofs/space/daisy_001/users/data/Gilad/gep_CL/';
FNBase='meas_MID01944_FID43885_ep2d_ge_EPTI_1p9_3and5shot_calib';

filename = [directory_rawdata,FNBase,'.dat'];
save_filename = 'Calib';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%% Get Parameters
% filename = [directory_rawdata,'meas_MID04688_FID21650_ep2d_ge_EPTI_1p3_7shot_calib.dat'];
% save_filename = 'Human_Calib_1p3_7shot';

% read parameters
pf_echo=0;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo);
Nreps = meas.evp.RawRep;
Necho=meas.evp.NEcoMeas;
parameters=meas.prot;
clear meas.data;
nslice_group=meas.evp.NSlcMeas;
generate_calib = 1;


for slice= 1:nslice_group
% for slice= 6
    %% Process of Refscan
    RepsToRead=1:Nreps;
    kdata  = EPTI_SMS_Preprocess_Imgscan_memorySave_GESE(filename,1,RepsToRead,SMS_data,slice,pf_echo);
    kdata=kdata(3:3+22,:,:,:,:);
    kdata_calib=single(kdata);
    
%     save([OutP,'data/Data_acq/',save_filename,'_Slice',num2str(slice),'.mat'],'kdata_calib','-v7.3');
    save([OutP,save_filename,'_Slice',num2str(slice),'.mat'],'kdata_calib','-v7.3');
    % generating calibration data
    if(generate_calib == 1)
        sData=size(kdata);
        kdata=crop(kdata,sData);
        [kdata_calib,B0fit] = generate_Calib(kdata,parameters,0.3);
        kdata_calib = single(kdata_calib);

        k=permute(kdata_calib,[2 3 1 4]);
        im=ifft2c(k);
%         figure; imshow(permute(sos(im(:,end:-1:1,14,:),4),[2 1 3]),[]);

%         save([OutP,'data/Data_acq/',save_filename,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib','B0fit','-v7.3');
        save([OutP,save_filename,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib','B0fit','-v7.3');
    end
end
%%
parameters.nechoGE=size(kdata,1);
save([OutP,'meas_prot_',save_filename,'.mat'],'parameters');
disp(['Saved ' OutP,'meas_prot_',save_filename,'.mat']);
%%
load([OutP,'meas_prot_',save_filename,'.mat'],'parameters');
%%
slice=9;
load([OutP,save_filename,'_Slice',num2str(slice),'_generated.mat']);
%%
SlbLoc=[[parameters.sSliceArray.sPosition_dSag];[parameters.sSliceArray.sPosition_dCor];[parameters.sSliceArray.sPosition_dTra]];
%%
nRO=size(B0fit,1);
nPE=size(B0fit,1);
% nRO=size(k,1);
% nPE=size(k,1);
ES_ms=parameters.iEffectiveEpiEchoSpacing/1e3;
%%
k=permute(kdata_calib,[2 3 1 4]);
PadUp=ceil((nRO-size(k,2))/2);
PadDown=floor((nRO-size(k,2))/2);
kExt=padarray(padarray(k,[0 PadDown 0 0],'post'),[0 PadUp 0 0],'pre');
idata_calib=ifft2c(k);
iExt=ifft2c(kExt);

try
    CurSensC=RunESPIRiTForSensMapsMultiMap(squeeze(iExt(:,:,1,:)),0,[nRO,nPE]);
catch
    CurSensC=RunESPIRiTForSensMapsMultiMap(squeeze(iExt(:,:,1,:)),14,[nRO,nPE]);
end

CurSensC=CurSensC(:,:,:,1);
SensMskC=grmss(CurSensC,3)>0.01;

iExtC=sum(iExt.*conj(perm43(CurSensC)),4);
disp('ok 1');
%%
nSlices=parameters.sSliceArray_lSize;
FirstTE_ms=9;
WhichEchosToUse=3:12;

for slice=1:nSlices
    disp(slice);
    load([OutP,save_filename,'_Slice',num2str(slice),'_generated.mat']);
    
    k=permute(kdata_calib,[2 3 1 4]);
    kExt=padarray(padarray(k,[0 PadDown 0 0],'post'),[0 PadUp 0 0],'pre');
    idata_calib=ifft2c(k);
    iExt=ifft2c(kExt);
    
    try
        CurSensC=RunESPIRiTForSensMapsMultiMap(squeeze(iExt(:,:,1,:)),0,[nRO,nPE]);
    catch
        CurSensC=RunESPIRiTForSensMapsMultiMap(squeeze(iExt(:,:,1,:)),14,[nRO,nPE]);
    end
    CurSensC=CurSensC(:,:,:,1);
    CurSensCS(:,:,:,slice)=CurSensC;
    SensMskCS(:,:,slice)=grmss(CurSensCS(:,:,:,slice),3)>0.01;
    
    idata_calibe=imresizeBySlices(idata_calib,[nRO,nPE]);
    idata_calibeC=sum(idata_calibe.*conj(perm43(CurSensC)),4);
    
    [PDBaseCseS(:,:,slice), UpdatedB0Map_HzCseS(:,:,slice), UpdatedT2SMap_msCseS(:,:,slice), s_valsCseS(:,:,:,slice),...
        Fitted0CseS(:,:,:,slice), PDBase0CseS(:,:,slice)]=...
        FitToModel_MPBD1CSf(idata_calibeC,WhichEchosToUse,ES_ms,FirstTE_ms);
end
CurSensCS=single(CurSensCS);
%%
save([OutP 'Sens.mat'],'CurSensCS','SensMskCS');
save([OutP 'Fit.mat'],'PDBaseCseS','UpdatedB0Map_HzCseS','UpdatedT2SMap_msCseS','s_valsCseS','PDBase0CseS');
%% That's it for calib, Sens, B0fit, Fit





%%
CurSensCs=RunESPIRiTForSensMapsMultiMap(squeeze(idata_calib(:,:,1,:)),0,gsize(idata_calib,1:2));
CurSensCs=CurSensCs(:,:,:,1);
SensMskCs=grmss(CurSensCs,3)>0.01;
idata_calibC=sum(idata_calib.*conj(perm43(CurSensCs)),4);

idata_calibCext=imresizeBySlices(idata_calibC,[nRO,nPE]);
CurSensCse=imresizeBySlices(CurSensCs,[nRO,nPE]);
idata_calibe=imresizeBySlices(idata_calib,[nRO,nPE]);
idata_calibeC=sum(idata_calibe.*conj(perm43(CurSensC)),4);

idata_calibeCS=SmoothBySlices(idata_calibeC,[20 20],5);
disp('ok 2');
%%
FirstTE_ms=9;

WhichEchosToUse=3:12;
% WhichEchosToUse=10:14;
[PDBaseCse, UpdatedB0Map_HzCse, UpdatedT2SMap_msCse, s_valsCse, Fitted0Cse, PDBase0Cse]=...
    FitToModel_MPBD1CSf(idata_calibeC,WhichEchosToUse,ES_ms,FirstTE_ms);
%     FitToModel_MPBD1CSf(idata_calibCext,WhichEchosToUse,ES_ms,FirstTE_ms);
ShowAbsAngle(PDBase0Cse,1,[0 20]);title(FirstTE_ms);
fgmontage(UpdatedB0Map_HzCse,[-400 400]);

Msk=(mPDBase0Cse>4);

mPDBase0Cse=abs(PDBase0Cse);
mPDBase0Cse=min(mPDBase0Cse,median(mPDBase0Cse(:))*20);
mPDBase0Cse=SmoothBySlices(mPDBase0Cse,[20 20],5);
PDBase0Csex=mPDBase0Cse.*exp(1i*angle(PDBase0Cse));
% PDBase0CsexS=SmoothBySlices(PDBase0Csex,[20 20],5);
ShowAbsAngle(PDBase0Csex,1,[0 20]);title(FirstTE_ms);

UpdatedB0Map_HzCsex=UpdatedB0Map_HzCse;
UpdatedB0Map_HzCsex(UpdatedB0Map_HzCsex>150)=0;
fgmontage(UpdatedB0Map_HzCsex,[-400 400]);

UpdatedT2SMap_msCsex=min(UpdatedT2SMap_msCse,200);
UpdatedT2SMap_msCsex=SmoothBySlices(UpdatedT2SMap_msCsex.*Msk,[20 20],5);
fgmontage(UpdatedT2SMap_msCsex,[0 100]);

%%
[PDBaseC, UpdatedB0Map_HzC, UpdatedT2SMap_msC, s_valsC, Fitted0C, PDBase0C]=...
    FitToModel_MPBD1CSf(iExtC,WhichEchosToUse,ES_ms,FirstTE_ms);
ShowAbsAngle(PDBase0C,1,[0 20]);title(FirstTE_ms);
%%
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
Ch2D=reshape(k,prod(gsize(k,1:3)),nChannels);
[~,S,sccmtx] = svd(Ch2D(1:1:end,:),'econ');
clear Ch2D 

ncc=11;
kCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(k,[4 1 2 3])),[4 1 2 3]);
kCCExt=padarray(padarray(kCC,[0 40 0 0],'post'),[0 41 0 0],'pre');

CurSensCCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSensC,[3 1 2])),[3 1 2]);
disp('ok cc');
%%
sig=permute(kCCExt,[1 2 5 4 6 7 3]);
SensCSMap=perm43(CurSensCCC);
mask_sample=abs(sig(:,:,:,1,:,:,:))>0;
Sz=[nRO nPE];
nEchos=size(sig,TS_Dim);
nCS=1;
%% Create linear operators
F_for=@(x) fft2cg(x);
F_adj=@(x) ifft2cg(x);
M_for=@(x) x.*mask_sample;
M_adj=@(x) x.*mask_sample;
S_for=@(x) sum(x.*SensCSMap,CS_Dim);
S_adj=@(x) sum(x.*conj(SensCSMap),Ch_Dim);

A_for=@(x) M_for(F_for(S_for(x)));
A_adj=@(x) S_adj(F_adj(M_adj(x)));
Aop=OpFromFunc(A_for,A_adj);

DataSize=FillOnesTo16(Sz);
DataSize(Ch_Dim)=ncc;
DataSize(TS_Dim)=nEchos;
ImSize=FillOnesTo16(Sz);
ImSize(TS_Dim)=nEchos;
ImSize(CS_Dim)=nCS;

OperatorTest2(Aop,ImSize,DataSize);
%%
x = randn(ImSize) + 1j*randn(ImSize);
yy = randn(DataSize) + 1j*randn(DataSize);
Ax = Aop*x;
Aty = Aop'*yy;
Out=abs(x(:)'*Aty(:) - conj(yy(:)'*Ax(:)));

LS_ScriptFN=[pwd filesep 'Cart2CS_ITS.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% Copy with sens/CS mask, sum over CS, FFT and sample mask
ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);

ImSz16=FillOnesTo16(ImSize);
Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample);
grmss(Ax-Ax_LS)
%%
TEs_ms=(FirstTE_ms+(0:(nEchos-1))*ES_ms).';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms=TEs_ms.';
EchoTimes_ms3=permute32(EchoTimes_ms);
%%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');
%%






%%
BaseSP='/autofs/space/daisy_002/users/Gilad/';

ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
WriteLinopToFile(ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;
BARTS=struct();
BARTS.cmd=['picsS -w 1 -R W:3:0:0.001 -m ' ScriptFN];
BARTS.ImSz16=MPImSize16;
BARTS.sig=sig;
BARTS.Others={SensCSMap mask_sample};

BaseFP='/tmp/';

Flds=setdiff(fieldnames(BARTS),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTS.(Flds{i})))
        for j=1:numel(BARTS.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTS.(Flds{i}){j});
            BARTS.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTS.(Flds{i}));
        BARTS.(Flds{i})=tmpFN;
    end
end

BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap mask_sample};

% save('ForCall.mat','ToBARTP','BARTS_Aopx');

BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

disp('Prepared for BART_Aop');
%%
c0=PDBase0Csex;
b0=UpdatedB0Map_HzCsex;
t0=UpdatedT2SMap_msCsex;
t0=min(max(t0,5),200);
b0=min(max(b0,-400),400);
m0=abs(c0);
m0=min(m0,median(m0(:))*20);
p0=angle(c0);
c0=m0.*exp(1i*p0);

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Est0=Mm0.*expPp0.*expBb0.*expTt0;

% Est0_nor=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Est0,BARTS_Aop.Others{:});
Est0_sig=bart(['linopScript ' LS_ScriptFN],BARTS_Aop.ImSz16,Est0,BARTS_Aop.Others{:});

BaseFac=grmss(sig)./grmss(Est0_sig);
c0=c0*BaseFac;
m0=m0*BaseFac;
Mm0=M*m0;
Est0=Mm0.*expPp0.*expBb0.*expTt0;

Est0_nor=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Est0,BARTS_Aop.Others{:});
Est0_sig=bart(['linopScript ' LS_ScriptFN],BARTS_Aop.ImSz16,Est0,BARTS_Aop.Others{:});

grmss(Est0_sig)/grmss(sig)
grmss(Est0_nor)/grmss(ksp_adj)
%%
ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
Elem0={m0,p0,b0,t0};
ElemL={M,P,B,TT};

fIdent=@(x) x;
fOne=@(x) 1;
fExp=@(x) exp(x);
fInv=@(x) 1./x;
fMinusInvSqr=@(x) -1./(x.^2);

Elem_u={fIdent,fIdent,fIdent,fInv};
% Elem_du={fOne,fOne,fOne,fMinusInvSqr};
Elem_s={fIdent,fExp,fExp,fExp};
% Elem_ds={fOne,fExp,fExp,fExp};

for i=1:nElements
    Elem_sLu{i}=Elem_s{i}(ElemL{i}*(Elem_u{i}(Elem0{i})));
end

Est0a=1;
for i=1:nElements
    Est0a=Est0a.*Elem_sLu{i};
end
Est0Xa=squeeze(sum(Est0a,CS_Dim));

MPBT=RefValueMPBT;
for i=1:nElements
    MPBT.(ElemNames{i})=Elem0{i};
end

disp('ok initials');
%%
ElemsAlphas=[0.0323, 2e-3, 1e-1, 0.5];
ElemsLambda=[0.02,0.05,0.05,0.0003];

ElemsAlphas=[0.0323, 1e-3, 1e-1, 1];
ElemsLambda=[0.05,0.05,0.1,0.003];

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
end
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
for i=1:numel(ElemsL)
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
% Types: M=1 P=2 D/B=3 T=4 C=5
ElementTypes=[1 2 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[1 1 1 1]*2;
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%%
ContinueRun=true;
if(ContinueRun)
    % Continue?
    for i=1:numel(Elem0)
        writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end
QQ=bart(['splitProx -i 150 -s 60 -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

% XDimsOut=readcfl([ToBARTP 'XDimsOut']);
% mDimsOut=readcfl([ToBARTP 'mDimsOut']);
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<1,1))-1);
% figure;plot(ErrVec)
%
Maps=cell(1,4);
for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
ElemRanges={[0 10],[-pi pi],[-400 400],[0 100]};
%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i}.*(m0>mThresh),ElemRanges{i});removeTicks;
%     gmontage(Elem0{i}.*(m0>mThresh),ElemRanges{i});removeTicks;
end
%%
CurElem=3;
ElemRanges={[0 700],[-pi pi],[-10 10],[0 100]};
ElemRecs=cat(3,ElemGT{CurElem},Maps{CurElem},Elem0{CurElem});
% ElemRecs=cat(3,ElemGT{CurElem},Maps{CurElem},Elem0{2});
% ElemRecs=angle(exp(1i*ElemRecs));
ElemRecs=ElemRecs-ElemGT{CurElem};
% ElemRecs=ElemRecs.*(m0>50);
fgmontagex(ElemRecs,ElemRanges{CurElem},'Size',[1 3]);title([PadStringWithBlanks('GT',70) PadStringWithBlanks('splitProx',70) 'Start point'],'FontSize',16);


%% Try CBT
% ElemsAlphas=[0.0323 1.3036e-5 0.0015];
ElemsAlphas=[0.000323 0.3036e-6 0.00003];
ElemsAlphas=ElemsAlphas*50;
ElemsLambda=[25e-2,0.05,0.0003];

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';

ElemsL={T,-1i*2*pi*TDx/1000,-TDx};

Elem0{1}=Elem0{1}.*exp(1i*Elem0{2});
Elem0{2}=Elem0{3};
Elem0{3}=Elem0{4};
Elem0=Elem0(1:3);
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
% Types: M=1 P=2 D/B=3 T=4 C=5
ElementTypes=[5 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

CurTestElemIdx=[1:3];
% CurTestElemIdx=2;
ninneriterBART=[0 0 0 0];
for i=1:numel(CurTestElemIdx)
    ninneriterBART(CurTestElemIdx(i))=2;
end
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:3
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%
writecfl([ToBARTP 'GT'],GT);
%

% Continue?
% for i=1:numel(Maps)
%     writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
% end
QQ=bart(['splitProx -i 2000 -s 60 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

% XDimsOut=readcfl([ToBARTP 'XDimsOut']);
% mDimsOut=readcfl([ToBARTP 'mDimsOut']);
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<1,1))-1);
% figure;plot(ErrVec)
%
Maps=cell(1,3);
for i=1:3
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
Maps{4}=angle(Maps{1});
Maps{1}=abs(Maps{1});
Maps=Maps([1 4 2 3]);

CurElem=CurTestElemIdx;
ElemRanges={[0 700],[-pi pi],[-400 400],[0 100]};
%
figure;
for i=1:4
    subplot(2,2,i);
%     gmontage(Maps{i}.*(m0>50),ElemRanges{i});removeTicks;
%     gmontage(Elem0{i}.*(m0>50),ElemRanges{i});title('WS');removeTicks;
    gmontage(ElemGT{i},ElemRanges{i});title('GT');removeTicks;
end