clear;
close all;
%%
addpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions'));
rmpath(genpath('/autofs/cluster/kawin/FuyixueWang/EPTI/Functions/bart-0.2.06/'));
%%
% directory = '/autofs/cluster/kawin/FuyixueWang/EPTI/EPTI_Rec_pack/';
% directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
directory_rawdata = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
filename_calib = 'Calib';
%% Get Parameters
FNBase='meas_MID04680_FID21642_ep2d_ge_EPTI_1p9_3shot_4dyns';
FNBaseCalib='meas_MID04678_FID21640_ep2d_ge_EPTI_1p9_calib';
% filename_calib = 'Human_Calib_1p3_5shot';
% filename_calib = 'Human_Calib_1p9';

FNBase='meas_MID00872_FID32107_ep2d_ge_EPTI_1p9_1shot_4dyns';
FNBaseCalib='meas_MID00870_FID32105_ep2d_ge_EPTI_1p9_1shot_calib';

FNBase='meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns';
% % FNBase='meas_MID00878_FID32113_ep2d_ge_EPTI_1p9_5shot_4dyns';
FNBaseCalib='meas_MID00874_FID32109_ep2d_ge_EPTI_1p9_3and5shot_calib';
% % 
% FNBase='meas_MID00903_FID32138_ep2d_ge_EPTI_1p9_3shot_4dyns_Cor';
% % FNBase='meas_MID00908_FID32143_ep2d_ge_EPTI_1p9_5shot_4dyns_Cor';
% FNBaseCalib='meas_MID00905_FID32140_ep2d_ge_EPTI_1p9_3and5shot_calib_Cor';

OutP = [directory_rawdata FNBase filesep];
mkdir(OutP);
system(['chmod +777 -R ' OutP]);
disp([OutP ' Created']);

CalibP=[directory_rawdata FNBaseCalib filesep];
filename = [directory_rawdata,FNBase,'.dat'];
%%
apodization_para=0;

% read parameters
pf_echo=0;
MB_factor=1;
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0;
%%
[meas] = read_meas_dat_memmap_EPTI_GESE_SMSnoRef(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo,MB_factor);

Nseg = meas.prot.sWipMemBlock_alFree(34);
Nreps = Nseg;
Ndyn = meas.evp.RawRep/Nreps;
Necho=meas.evp.NEcoMeas;
parameters=meas.prot;
clear meas.data;
% Rseg=meas.prot.sWipMemBlock_alFree(33);
% Rpe=meas.prot.sWipMemBlock_alFree(35);
Rseg=parameters.sWipMemBlock_alFree(33);
Rpe=parameters.sWipMemBlock_alFree(35);
Ngroup=floor(parameters.lPhaseEncodingLines/(Rseg/Rpe));
MB=1;
param.Nseg=Nseg;
param.Rseg=Rseg;
param.Rpe=Rpe;
param.Ngroup=Ngroup;
parameters.nechoGE=meas.evp.NLinMeas;

Nseg = parameters.sWipMemBlock_alFree(34);
Nreps = Nseg;
%%
% parameters.nechoGE=size(kdata,1);
save([OutP,'parameters.mat'],'parameters');
%%
parametersCalib=load([CalibP,'meas_prot_',filename_calib,'.mat'],'parameters');
parametersCalib=parametersCalib.parameters;
SlbLocCalib=[[parametersCalib.sSliceArray.sPosition_dSag];[parametersCalib.sSliceArray.sPosition_dCor];[parametersCalib.sSliceArray.sPosition_dTra]];

%%
load([OutP,'parameters.mat'],'parameters');
SlbLoc=[[parameters.sSliceArray.sPosition_dSag];[parameters.sSliceArray.sPosition_dCor];[parameters.sSliceArray.sPosition_dTra]];

Rseg=parameters.sWipMemBlock_alFree(33);
Rpe=parameters.sWipMemBlock_alFree(35);

Ngroup=floor(parameters.lPhaseEncodingLines/(Rseg/Rpe));
MB=1;
Rpe=parameters.sWipMemBlock_alFree(35);
nSlices=numel(parameters.sSliceArray);
%% Slice Loop
nslice_group=meas.evp.NSlcMeas/MB;
slice_sep=nslice_group;

nCh=numel(parameters.sCoilElementID_tElement);
nPE=parameters.lPhaseEncodingLines-1;
nRO=parameters.lBaseResolution;
for dyn=1:4
    kdata_GES=single(zeros([nPE,nRO,nRO,nCh,nslice_group]));
    RepsToRead=1:Nreps;
    RepsToRead=RepsToRead+(dyn-1)*Nreps;
    [kdata_GE_all] = EPTI_SMS_Preprocess_Imgscan_SMS1_GESE(filename,1,RepsToRead,SMS_data,1:nslice_group,pf_echo);
%% reconstruction loop of all slices
%     k_recon=zeros(meas.evp.NImageLins,Nseg*Rseg,nslice_group*MB,Ngroup*(Rseg/Rpe),meas.evp.NChaMeas);
%     kdata_correct=zeros(size(k_recon));
    for slice=1:nslice_group
        disp(['--- dyn ' num2str(dyn) ', slice ---- ',num2str(slice)]);
        tic;
        kdata_GE=kdata_GE_all(:,:,:,:,slice);
        [kdata_GE] = putRawtoPattern(kdata_GE,Nseg,Rseg,Rpe,Ngroup,'GE',MB);
        kdata_GES(:,:,:,:,slice)=kdata_GE;
    end
    kdata_GES=single(kdata_GES);
    save([OutP 'kdata_GES_dyn' num2str(dyn) ' .mat'],'kdata_GES','-v7.3');
end
%%
% load([OutP 'kdata_GES.mat'],'kdata_GES');
%%
kdata_GE=kdata_GES(:,:,:,:,1);

kMsk2D=abs(squeeze(kdata_GE(:,10,:,1)).')~=0;


nRO=size(kdata_GE,2);
nPE=size(kdata_GE,3);
nEchos=size(kdata_GE,1);
nChannels=size(kdata_GE,4);
% nSlices=nslice_group;
nSlices=parameters.sSliceArray_lSize;
% nShots=Nseg;
nShots=nPE/parameters.sWipMemBlock_alFree(33);

% Nseg = parameters.sWipMemBlock_alFree(34);
% Nreps = Nseg;
% % Ndyn = parameters.RawRep/Nreps;
% % Necho=meas.evp.NEcoMeas;
% parameters=meas.prot;
% clear meas.data;
% Rseg=parameters.sWipMemBlock_alFree(33);
% Rpe=parameters.sWipMemBlock_alFree(35);
% Ngroup=floor(parameters.lPhaseEncodingLines/(Rseg/Rpe));
% MB=1;
% param.Nseg=Nseg;
% param.Rseg=Rseg;
% param.Rpe=Rpe;
% param.Ngroup=Ngroup;
% parameters.nechoGE=meas.evp.NLinMeas;
%% Fuyixue GRAPPA "Script1"
for dyn=1:4
    clear kdata_GES
    load([OutP 'kdata_GES_dyn' num2str(dyn) ' .mat']);
    for slice=1:nSlices
        disp([dyn slice]);
        disp(datestr(now));
        kdata_GE=kdata_GES(:,:,:,:,slice);
        
        k=squeeze(sum(kdata_GE(1:16,:,:,:),1));
        %         figure; imshow(permute(sos(ifft2c(k)),[2 1]),[0 5]); colormap('gray');
        
        echotype='GE';
        %         load([directory,'data/Data_acq/',filename_calib,'_Slice',num2str(slice),'_generated.mat'],'kdata_calib');
        load([CalibP,filename_calib,'_Slice',num2str(slice),'.mat'],'kdata_calib');
        kdata_calib=double(kdata_calib);
        
        [recon] = EPTI_recon_SMS1_GenCalib(kdata_GE,kdata_calib,echotype,param);
        %     close all;
        %     im_recon=ifft2c(recon);
        %     im_recon_show = sos(im_recon(:,:,13:4:end-14,:),4);
        %     figure; imshow3(im_recon_show,[0 8],[2,size(im_recon_show,3)/2]);
        k_recon(:,:,1,:,:)=single(recon);
        [tmp,B0_variation_esti] = B0_variationCorrection(recon,param,parameters,0);
        kdata_correct(:,:,1,:,:) = single(tmp);
        disp(datestr(now));
        %     reconS(:,:,:,:,slice)=recon;
        %     B0_variation_estiS(:,:,:,slice)=B0_variation_esti;
        %     kdata_correctS(:,:,:,:,slice)=single(tmp);
        
        disp(['Post Processing Start: apodization and coil combination']);
        disp(datestr(now));
        [im_EPTI] = fMRI_PostRecon_process(k_recon,apodization_para);
        [im_EPTI_correct] = fMRI_PostRecon_process(kdata_correct,apodization_para);
        disp(datestr(now));
        
        save([OutP,'Recon_EPTI','_Dyn',num2str(dyn),'_Slice_',num2str(slice),'_',echotype,'.mat'],'im_EPTI','im_EPTI_correct',...
            'recon','kdata_correct','B0_variation_esti','-v7.3');
    end
end
disp('Finished GRAPPA Recon');
%%
load([OutP,'parameters.mat'],'parameters');
echotype='GE';
dyn=2;
slice=6;
slice=9;
load([OutP,'Recon_EPTI','_Dyn',num2str(dyn),'_Slice_',num2str(slice),'_',echotype,'.mat']);
[nRO,nPE,nEchos,nChannels]=size(recon);
ES_ms=parameters.iEffectiveEpiEchoSpacing/1e3;
FirstTE_ms=9;
% dThickness: 92
% dPhaseFOV: 228
% dReadoutFOV: 228
%%
IdxForSens=20;
im_recon_forSens=squeeze(ifft2c(recon(:,:,IdxForSens,:)));
% disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
CurSens=RunESPIRiTForSensMapsMultiMap(im_recon_forSens,0,[nRO,nPE]);
CurSens=CurSens(:,:,:,1);
SensMsk=grmss(CurSens,3)>0.01;
%%
RecAllEchosChannels=ifft2c(recon);
RecAllEchos=sum(RecAllEchosChannels.*conj(perm43(CurSens)),4);
%%
WhichEchosToUse=15:55;
% WhichEchosToUse=10:14;
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(RecAllEchos,WhichEchosToUse,ES_ms,FirstTE_ms);
%%
ShowAbsAngle(PDBase0,1,[0 20]);title(FirstTE_ms);

ShowAbsAngle(gflip(PDBase0,2),1,[0 20]);title(FirstTE_ms);

fgmontage(gflip(Fitted0(:,:,11:11:end),2));title('GRAPPA, 1CS Fitted');
%%
load([CalibP 'Sens.mat'],'CurSensCS','SensMskCS');
load([CalibP 'Fit.mat'],'PDBaseCseS','UpdatedB0Map_HzCseS','UpdatedT2SMap_msCseS','s_valsCseS','PDBase0CseS');
%% Run on all slices
slice=9;
for slice=[1:nSlices]
    disp(slice);
    load([OutP,'Recon_EPTI','_Dyn',num2str(dyn),'_Slice_',num2str(slice),'_',echotype,'.mat']);
    IdxForSens=20;
    im_recon_forSens=squeeze(ifft2c(recon(:,:,IdxForSens,:)));
    % disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
    CurSens=RunESPIRiTForSensMapsMultiMap(im_recon_forSens,0,[nRO,nPE]);
    CurSens=CurSens(:,:,:,1);
    SensMsk=grmss(CurSens,3)>0.01;
    CurSensSG(:,:,:,slice)=CurSens;
    SensMskSG(:,:,slice)=SensMsk;
    CurSens=CurSensSG(:,:,:,slice);
    %
    RecAllEchosChannels=ifft2c(recon);
    RecAllEchos=sum(RecAllEchosChannels.*conj(perm43(CurSens)),4);
    RecAllEchosS(:,:,:,slice)=RecAllEchos;
    %
    WhichEchosToUse=15:55;
    % WhichEchosToUse=10:14;
    [PDBaseSG(:,:,slice), UpdatedB0Map_HzSG(:,:,slice), UpdatedT2SMap_msSG(:,:,slice), s_valsSG(:,:,:,slice),...
        FittedSG(:,:,:,slice), PDBase0SG(:,:,slice)]=...
        FitToModel_MPBD1CSf(RecAllEchosS(:,:,:,slice),WhichEchosToUse,ES_ms,FirstTE_ms);
end
%%
save([OutP 'Sens.mat'],'CurSensSG','SensMskSG');
save([OutP 'GFit.mat'],'PDBaseSG','UpdatedB0Map_HzSG','UpdatedT2SMap_msSG','s_valsSG','PDBase0SG');
save([OutP 'RecAllEchosS.mat'],'RecAllEchosS');
%%
load([OutP 'RecAllEchosS.mat'],'RecAllEchosS');
load([OutP 'Sens.mat'],'CurSensSG','SensMskSG');
load([OutP 'GFit.mat'],'PDBaseSG','UpdatedB0Map_HzSG','UpdatedT2SMap_msSG','s_valsSG','PDBase0SG');
%% Recon basic
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
%%
for dyn=1:4
    clear kdata_GES
    load([OutP 'kdata_GES_dyn' num2str(dyn) ' .mat']);
    for SliI=1:nSlices
        disp([dyn SliI]);
        kdata_GE=kdata_GES(:,:,:,:,SliI);
        kdata_GE=permute(kdata_GE,[2 3 1 4]);
        
        CurSens=CurSensSG(:,:,:,SliI);
        %
        %     if(dyn==1)
        %     Ch2D=reshape(kdata_GE,prod(gsize(kdata_GE,1:3)),nChannels);
        %     [~,S,sccmtx] = svd(Ch2D(1:1:end,:),'econ');
        %     sccmtxS(:,:,SliI)=sccmtx;
        %     else
        sccmtx=sccmtxS(:,:,SliI);
        %     end
        % clear Ch2D
        %     figure;plot(cumsum(diag(S).^2)./sum(diag(S).^2))
        
        ncc=31;
        kCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(kdata_GE,[4 1 2 3])),[4 1 2 3]);
        
        sig=permute(kCC,[1 2 5 4 6 7 3]);
        sigS(:,:,SliI,:,:,:,:,dyn)=sig;
        
        if(dyn==1)
            CurSensCCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2])),[3 1 2]);
            SensCSMap=perm43(CurSensCCC);
            SensCSMapS(:,:,SliI,:)=SensCSMap;
        end
    end
end
mask_sample=abs(sig(:,:,:,1,:,:,:,1))>0;
MskCC=perm75(repmat(mask_sample,[1 1 1 ncc 1 1 1 1 1]));
for dyn=1:4
    for SliI=1:nSlices
        tmp=sigS(:,:,SliI,:,:,:,:,dyn);
        SigMskd(:,SliI,dyn)=single(tmp(MskCC));
    end
end
clear tmp
save([OutP 'SigMskd_dyns.mat'],'SigMskd','MskCC','mask_sample');


Sz=[nRO nPE];
nEchos=size(sig,TS_Dim);
nCS=1;
clear kCC
save([OutP 'SensCSMapS.mat'],'SensCSMapS','sccmtxS');
%%
TEs_ms=(FirstTE_ms+(0:(nEchos-1))*ES_ms).';
NTEs=TEs_ms-TEs_ms(1);
NTEs=NTEs./NTEs(2);
NTEs=NTEs(:);
TEs_ms3=permute(TEs_ms,[1 3 2]);

EchoTimes_ms=TEs_ms.';
EchoTimes_ms3=permute32(EchoTimes_ms);
%
T=grepmat(gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]),nEchos,TS_Dim);
TD=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute(NTEs,[TS_Dim 1]);
TDx=gpermute(eye(nCS),[Ch_Dim CS_Dim 1 2]).*gpermute((EchoTimes_ms.'),[TS_Dim 1]);
M = fmac(ones(nEchos,1), T,Ch_Dim,[CS_Dim TS_Dim]);
P = fmac(NTEs(1:nEchos), 1i*T,Ch_Dim,[CS_Dim TS_Dim]);
B = fmac(NTEs(1:nEchos)/1000, -1i*2*pi*TDx/1000,Ch_Dim,[CS_Dim TS_Dim]);
TT = fmac(NTEs, -TDx,Ch_Dim,[CS_Dim TS_Dim]);
disp('ok operators');

ImSize=FillOnesTo16(Sz);
ImSize(TS_Dim)=nEchos;
ImSize(CS_Dim)=nCS;

ImSz16=FillOnesTo16(ImSize);

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
%% That's it. Now splitProx

%% To MLN
load([OutP 'SensCSMapS.mat'],'SensCSMapS','sccmtxS');
load([OutP 'SigMskd.mat'],'SigMskd','MskCC','mask_sample');
%%
BaseBaseOutP=[OutP 'MLN/'];
mkdir(BaseBaseOutP);
system(['chmod -R 777 ' BaseBaseOutP]);
disp(['Created ' BaseBaseOutP]);
%%
nCh=13;
SliI=7;
for SliI=1:nSlices
    disp(SliI);
    
    tmpSig=MskCC*0;
    tmpSig(MskCC)=SigMskd(:,SliI);
    clear SigX
    I1=[];I2=[];
    for i=1:nEchos
        tmpM=MskCC(:,:,1,1,i);
        [I1c, I2c]=find(tmpM);
        I1=[I1;I1c];
        I2=[I2;I2c];
        for j=1:nCh
            tmp=tmpSig(:,:,1,j,i);
            SigX(:,j,i)=tmp(tmpM);
        end
    end
    sig=CombineDims(SigX,[3 1]);
    sig=permute(sig,[3 1 4 2]);
    disp('Prepared sig for this slice');
    %
    
    
    
    Sensr=SensCSMapS(:,:,SliI,:); % 114   114     1    31
    %     sig=QQ.CurSig; % 1       47273           3          31
    
    % Sensr=imresize(Sensr,[116 116]);
    Sz128=gsize(Sensr,1:2);
    
    B0_Hz=UpdatedB0Map_HzSG(:,:,SliI);
    % B0_Hz=UpdatedB0Map_RS(:,:,WhichRSToUse);
    % B0_Hz=RefB0MrS(:,:,SliI);
    %     B0_Hz=UpdatedB0Map_RS_LLRS(:,:,2,SliI);
    B0_Hz=imresize(B0_Hz,Sz128);
    %
    BaseOutDir=[BaseBaseOutP 'S3Sli' num2str(SliI) filesep];
    mkdir(BaseOutDir);
    system(['chmod -R 777 ' BaseOutDir]);
    % batchSize=8;
    batchSize=16;
    nTS=15;
    SensCC=squeeze(Sensr(:,:,1:nCh)); % [X Y Ch]
    SensMsk=grmss(SensCC,3)>0.01; % [X Y]
    save([BaseOutDir 'SensCC1.mat'],'SensCC','SensMsk');
    
    % CurSig is 45501          13
    CurSig=squeeze(sig)*0.2;
    tmp=Row(CurSig);
    tmp2=[real(tmp) imag(tmp)];
    Data=tmp2;
    Data(batchSize,end)=0;
    save([BaseOutDir 'RealDataForNN.mat'],'Data');
    nTrajA=size(CurSig,1);
    TSBF=GetTSCoeffsByLinear(nTrajA,nTS).';
    
    Centers=linspacecenters(1,nTrajA,nEchos);
    nPointsPerEcho=nTrajA/nEchos;
    AA=TSBF(:,floor(Centers));
    TSBF2=kron(AA,ones(1,nPointsPerEcho));
    TSBF=TSBF2;
    
    save([BaseOutDir 'B0TS.mat'],'TSBF','B0_Hz');
    
    
    Sz128=gsize(SensCC,1:2);
    clear Trajm2
    Trajm2(1,:)=I1-Sz128(1)/2;
    Trajm2(2,:)=I2-Sz128(2)/2;
    [FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
    Kd=st.nufftStruct.Kd;
    SN=st.nufftStruct.sn;
    P=st.nufftStruct.p/sqrt(prod(Sz128));
    save([BaseOutDir 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');
    TotalAcqTime_ms=ES_ms*nEchos;
    TS_TimePoints=9+linspace(0,TotalAcqTime_ms,nTS);
    TSstr=strrep(num2str(TS_TimePoints,'%3.5f,'),' ','');
    TSstr=TSstr(1:end-1);
    % TSstr=['TimePoints_ms ' TSstr(1:end-1)];
    
    TS_TimePointsForRec=linspace(0,TotalAcqTime_ms,8);
    TS_TimePointsForRec_Str=strrep(num2str(TS_TimePointsForRec,'%3.5f,'),' ','');
    TS_TimePointsForRec_Str=TS_TimePointsForRec_Str(1:end-1);
    
    St=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForS10.txt');
    St.LabelsW=Sz128(1);
    St.LabelsH=Sz128(2);
    St.aDataW=Sz128(1);
    St.aDataH=Sz128(2);
    St.nTraj=size(CurSig,1);
    St.nccInData=nCh;
    St.RealDatancc=nCh;
    St.DataH=size(CurSig,1)*nCh*2;
    St.TimePoints_ms=TSstr;
    St.batch_size=batchSize;
    St.nToLoad=10000;
    St.TimePointsForRec_ms=TS_TimePointsForRec_Str;
    St.SessionNameBase=[OutP 'MLN/A3S' num2str(SliI)];
    St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
    St.BaseTSDataP=BaseOutDir;
    St.BaseNUFTDataP=BaseOutDir;
    % if(strcmp(HostName,'tiger'))
    % end
    % disp(TSstr)
    % ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForMES' num2str(SliI) '.txt'];
    ParamsOutFn=[BaseBaseOutP 'ParamsForME3S' num2str(SliI) '.txt'];
    Txt=gStruct2txt(St,ParamsOutFn);
    disp('ok')
    MSLines3{SliI}=['/usr/pubsw/bin/python3 /autofs/cluster/kawin/Gilad/TF/srezN/srez_main1.py ' ParamsOutFn];
end
%% End MLN






%%
load([OutP 'Sens.mat'],'CurSensSG','SensMskSG');
load([OutP 'GFit.mat'],'PDBaseSG','UpdatedB0Map_HzSG','UpdatedT2SMap_msSG','s_valsSG','PDBase0SG');
load([OutP 'RecAllEchosS.mat'],'RecAllEchosS');
%% CollectMLNResults
AllRes=cat(4,gflip(squeeze(RecAllEchosS(3:end-2,3:end-2,50,:)),2),squeeze(MLNRXS1(:,:,7,ROrd)),squeeze(MLNRXS3(:,:,7,ROrd)));
AllRes=AllRes./grms(AllRes,1:3);
%%
for i=1:nSlices
    load(['/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00884_FID32119_gSpi2d_T12_Dw11_d110_VD1_Cor/Maps_S' num2str(i) '.mat']);
    MapsS(:,:,:,i)=cat(3,Maps{:});
end
%%
SplitProxRes=squeeze(MapsS(:,:,1,ROrd).*exp(-1*42./MapsS(:,:,4,ROrd)));

AllRes=cat(4,gflip(squeeze(RecAllEchosS(3:end-2,3:end-2,50,:)),2),squeeze(MLNRXS1(:,:,7,ROrd)),squeeze(MLNRXS3(:,:,7,ROrd)),SplitProxRes);
AllRes=AllRes./grms(AllRes,1:3);
fgmontagex(gflip(AllRes(:,:,[5 12],:),1))
%%
% load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/meas_MID04675_FID21637_gre_te4_9/sTwixX.mat');
% GRE_Ref_FN='meas_MID04675_FID21637_gre_te4_9';
GRE_Ref_FN='meas_MID00856_FID32091_gre_4echo_32_22'; % 68x
load([directory_rawdata GRE_Ref_FN '/sTwixX.mat']);
RefLocs=load([directory_rawdata GRE_Ref_FN '/Locs.mat']);
RefLocs=RefLocs.RotatedLocs;

sTwixX.hdr.MeasYaps.sSliceArray.asSlice
sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition

FOVHere=[parameters.sSliceArray(1).dReadoutFOV parameters.sSliceArray(1).dPhaseFOV];
FOVRef=[sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV];
RefS=load([directory_rawdata GRE_Ref_FN '/B0T2S.mat']);

RefB0=gflip(rot90(padarray(imresize(padarray(RefS.B0M_Hz,[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefT2SBase=max(4,min(200,abs(RefS.UpdatedT2SMap_ms)));
RefT2S=gflip(rot90(padarray(imresize(padarray(RefT2SBase,[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefT2S=max(4,min(200,abs(RefT2S)));

Sz=gsize(PDBaseSG,1:2);


RefB0=gflip(rot90(padarray(imresize(padarray(RefS.B0M_Hz,[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
% RefB0=circshift(circshift(RefB0,2,1),4,2);
RefB0=circshift(circshift(RefB0,2,1),-8,2);
RefB0r=imresize(RefB0,Sz);

fgmontagex(s_valsSG(:,:,1,1:3:end));title('GRAPPA')
fgmontagex(RefSvalsr(:,:,3:3:end));title('Ref GRE-ME')

RefSvals=gflip(rot90(padarray(imresize(padarray(squeeze(RefS.s_vals(:,:,1,:)),[0 2],'Both'),FOVRef),(FOVHere-FOVRef)/2,'Both'),3),1);
RefSvals=circshift(circshift(RefSvals,2,1),-8,2);
RefSvalsr=imresize(RefSvals,Sz);
RefSvalsrX=(RefSvalsr(:,:,1:2:end)+RefSvalsr(:,:,2:2:end))/2;
fgmontagex(RefSvalsrX(:,:,1:3:end));title('RefSvalsrX')


B0MapToUse_Hz=RefB0r;
% B0MapToUse_Hz=imresize(RefB0,Sz);
B0MapToUse_Hz=UpdatedB0Map_HzSG;
% B0MapToUse_Hz=gflip(UpdatedB0Map_HzCseS,2);
% B0MapToUse_Hz=SmoothBySlices(gflip(UpdatedB0Map_HzCseS,2),[10 10],1);

Thickness_mm=parameters.sSliceArray(1).dThickness;
% DistBetweenSlices_mm=mean(diff(SlbLoc(3,:)));
DistBetweenSlices_mm=mean(diff(SlbLoc(2,:)));
%%
ExRefEcho=(RefS.PDBase0(:,:,SliI*2-1)).*exp(-20./RefS.UpdatedT2SMap_ms(:,:,SliI*2-1));
fgmontagex(gflip(rot90(ExRefEcho,3),1))
%%
dB0dx=symD(B0MapToUse_Hz,1);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dz=symD(B0MapToUse_Hz,3)*Thickness_mm/DistBetweenSlices_mm;
% dB0dz=symD(B0MapToUse_Hz,3);

ThroughPlanePhaseDiffS=2*pi.*EchoTimes_ms3.*perm43(dB0dz)/1000;
ThroughPlaneDecayS=sinc(ThroughPlanePhaseDiffS/(2*pi));
%%
SliI=7;
figure;
subplot(2,2,1);imagesc(B0MapToUse_Hz(:,:,SliI),[-200 200]);removeTicks;
subplot(2,2,2);imagesc(dB0dx(:,:,SliI),[-20 20]);removeTicks;
subplot(2,2,4);imagesc(dB0dz(:,:,SliI),[-20 20]);removeTicks;colormap gray
subplot(2,2,3);imagesc(dB0dy(:,:,SliI),[-20 20]);removeTicks;
h=colorbar('westoutside');
%%
ChosenP=[35 60];
% ChosenP=[64 92];
% ChosenP=[42 55];
DemoVal_Hz=dB0dz(ChosenP(2),ChosenP(1),SliI)
% sinc Sin(pi*x)/(pi*x) function.
% DemoVal_Hz=20;
RadsV=2*pi*EchoTimes_ms*DemoVal_Hz/1000;
% figure;plot(EchoTimes_ms,Rads);
SincV=sinc(RadsV/(2*pi));
SigV=abs(squeeze(RecAllEchos(ChosenP(2),ChosenP(1),:))).';
SigV=SigV*SincV(21)./SigV(21);
figure;plot(EchoTimes_ms,SincV,'k');
hold on;
plot(EchoTimes_ms,SigV,'r');
title(DemoVal_Hz);
%%
















%%
CurPrefix='qweqw';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
BaseSP=['/autofs/space/daisy_002/users/Gilad/' CurPrefix];
BaseFP=['/tmp/' CurPrefix];

LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);
% ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
% BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
% SigToUse=sig;
% disp('got ksp_adj');

% ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
% MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
% WriteLinopToFile(LS_ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;

BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap mask_sample};
%%












% figure;plot(EchoTimes_ms,SigV./SincV,'k');
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

OperatorTest2(Aop,ImSize,DataSize);
%%
% LS_ScriptFN=[pwd filesep 'Cart2CS_ITS.txt'];
% % SensCSMap is 0
% % mask_sampleP is 1
% % Copy with sens/CS mask, sum over CS, FFT and sample mask
% ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
% WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);
% 
% Ax_LS=bart(['linopScript ' LS_ScriptFN],ImSz16,x,SensCSMap,mask_sample);
% grmss(Ax-Ax_LS)
%%






%%
CurPrefix='qweqw';
ToBARTP=['/autofs/space/daisy_002/users/Gilad/gUM/' CurPrefix];
BaseSP=['/autofs/space/daisy_002/users/Gilad/' CurPrefix];
BaseFP=['/tmp/' CurPrefix];

LS_ScriptFN=[ToBARTP 'Cart_mCS_ITS_TSC.txt'];

ITS_Cmnds={['fmac 0 ' num2str(CS_Flag)],'fftc 3','fmac 1 0'};
WriteLinopToFile(LS_ScriptFN,ITS_Cmnds);
% ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
% BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);
% SigToUse=sig;
% disp('got ksp_adj');

% ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
% MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
% WriteLinopToFile(LS_ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;

BARTS_Aopx=struct();
BARTS_Aopx.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aopx.ImSz16=ImSz16;
BARTS_Aopx.Others={SensCSMap mask_sample};

% save('ForCall.mat','ToBARTP','BARTS_Aopx');

BARTS_Aop=WriteBARTStructToFiles(BARTS_Aopx,BaseFP);

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

disp('Prepared for BART_Aop');
%%
PDBase0=PDBase0SG(:,:,SliI);
UpdatedB0Map_Hz=UpdatedB0Map_HzSG(:,:,SliI);
UpdatedT2SMap_ms=UpdatedT2SMap_msSG(:,:,SliI);

% UpdatedB0Map_Hz=RefB0r(:,:,SliI);

% UpdatedB0Map_Hz=imresize(RefB0(:,:,SliI),Sz);
% UpdatedB0Map_Hz=circshift(circshift(UpdatedB0Map_Hz,2,1),4,2); % positive is right down
% UpdatedB0Map_Hz=circshift(circshift(UpdatedB0Map_Hz,2,1),2,2); % positive is right down
%%
% [PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
%     FitToModel_MPBD1CSf(RecAllEchos,WhichEchosToUse,ES_ms,FirstTE_ms);
m0=abs(PDBase0(:,:,1));
p0=angle(PDBase0(:,:,1));
m0=(min(median(m0(:))*10,m0));
b0=UpdatedB0Map_Hz;
t0=max(5,min(200,abs(UpdatedT2SMap_ms)));

Mm0=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Est0=Mm0.*expPp0.*expBb0.*expTt0;
Est0X=squeeze(sum(Est0,CS_Dim));

AHA_Est0=bart(['linopScript -N ' LS_ScriptFN],BARTS_Aop.ImSz16,Est0X,BARTS_Aop.Others{:});

AEst0=bart(['linopScript ' LS_ScriptFN],ImSz16,Est0X,SensCSMap,mask_sample);

AFac=grmss(sig)./grmss(AEst0);
m0=m0*AFac;

Elem0={m0,p0,b0,t0};

disp('Produced AEst0 and sig');
%%
ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
ElemL={M,P,B,TT};
Elem0={m0,p0,b0,t0};
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};

save([OutP 'ForMPBT_Sli' num2str(SliI) '.mat'],'Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');
disp(['Saved ' OutP 'ForMPBT_Sli' num2str(SliI) '.mat']);
%%

save('ForMPBT.mat','Elem0','ElemL','ElemsL','ksp_adj','BARTS_Aopx','sig','LS_ScriptFN','Sz','ThroughPlaneDecay');

%
ElemsAlphas=[0.0323, 1.5366e-07, 1.3036e-05, 0.0015];
ElemsAlphas=[0.323, 1e-03, 1e-1, 1e-0];
ElemsLambda=[0.015,0.05,0.05,0.0003]; % Very good for each one separately
ElemsLambda=[0.02,0.05,0.05,0.05];

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
end

for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[Sz 1 1 1 1 1]));
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
ElementTypes=[1 2 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[0 0 0 0];
for i=1:nElements
    ninneriterBART(i)=2;
end
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%% Continue?
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
end
%%
QQ=bart(['splitProx -i 500 -s 60 -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
%%
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
figure;plot(ErrVec)
%%
for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end

Mm0=M*Maps{1};
expPp0 = exp(P * Maps{2});
expBb0 = exp(B * Maps{3});
expTt0 = exp(TT * (1./Maps{4}));
Rec0=Mm0.*expPp0.*expBb0.*expTt0;
Rec0X=squeeze(sum(Rec0,CS_Dim));

disp('Loaded maps');
%%
ElemRanges={[0 10],[-pi pi],[-400 400],[0 100]};
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
    if(i==1), title('Optimized'); end
end
%%
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Elem0{i},ElemRanges{i});removeTicks;
    if(i==1), title('Start point'); end
end
%%
fgmontage(Maps{3}-Elem0{3},[-5 5])
%%
EchoIToShow=30;
% fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec_LLRX(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)),'Size',[1 3]);
fgmontagex(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow)));
fgmontagex(angle(cat(3,RecAllEchos(:,:,EchoIToShow),Rec0X(:,:,EchoIToShow))));

fgmontagex(RecAllEchos(:,:,EchoIToShow));
fgmontagex(Rec0X(:,:,EchoIToShow));
%%






















%% LLR
%% T2* histogram
%% Using components
T2svalues_ms=linspace(5,300,200);
Decays=exp(-EchoTimes_ms./(T2svalues_ms.'));

[Ud,Sd,Vd]=svd(Decays,'econ');
%%
figure;
for i=1:4
    subplot(2,2,i);
    plot(Vd(:,i),'--');hold on
end
legend({'T2* decays','W T2* decays'});
%%
% [PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
%     FitToModel_MPBD1CSf(RecAllEchos,WhichEchosToUse,ES_ms,FirstTE_ms);
% T2S_ms=40+0*max(5,min(200,abs(UpdatedT2SMap_ms)));

B0M_Hz=-UpdatedB0Map_Hz;
nComponentsToUse=2;

ScriptFN_CompgBo=[BaseSP 'CartCompgB0_N.txt'];
TSCOnlyM=exp(-EchoTimes_ms3./T2S_ms);
TSCOnlyP=exp(1i*2*pi.*B0M_Hz*(1e-3).*EchoTimes_ms3);
TSC=TSCOnlyP;
% TSC=TSCOnlyP.*ThroughPlaneDecay;

Sz16CompgB0=FillOnesTo16([Sz 1 1 1 nComponentsToUse]);
% CompsP=permute(WVd(:,1:nComponentsToUse),[7:-1:3 2 1]);
CompsP=permute(Vd(:,1:nComponentsToUse),[7:-1:3 2 1]);
SensP=SensCSMap;
MskP=mask_sample;
TSCP=perm73(TSC);
% # Img is [x y z 1 1 Comp]
% # file 0 is sensitivity maps [x y z Ch Maps]  
% # file 1 is sampling pattern/Trajectory [kx ky kz 1 1 1 TS]
% # file 2 is TSC [x y z 1 1 1 TS]
% # file 3 is Components [1 1 1 1 1 Comp TS]

WriteLinopToFile(ScriptFN_CompgBo,{'fmac 3 32','fmac 2 0','fmac 0 16','fftc 3','fmac 1 0'});
% THLR_lambda=0.1;
% LLR_lambda=10;
LLR_lambda=0.1;
% RhoStr='';
RhoStr=[' -u ' num2str(1e-3) ' '];
BlkSz=4;
disp('Prepared LLR');
%% By rep (160sec per rep)
nccToUse=11;
WhichEchosToUse=1:70;

Rec_LLR=bart(['picsS -m ' RhoStr ' -w 1 -b ' num2str(BlkSz) ' -R L:3:3:' num2str(LLR_lambda) ' ' ScriptFN_CompgBo],...
    Sz16CompgB0,sig(:,:,1,1:nccToUse,1,1,WhichEchosToUse),SensP(:,:,1,1:nccToUse),MskP(:,:,1,1,1,1,WhichEchosToUse),...
    TSCP(:,:,1,1,1,1,WhichEchosToUse),CompsP(:,:,1,1,1,:,WhichEchosToUse));
Rec_LLRX=sum(Rec_LLR.*CompsP,6);