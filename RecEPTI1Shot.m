%% GRAPPA, only 30 echos (1 zag)
for SliI=1:nSlices
disp([SliI SliI]);
    kdata_GE=kdata_GES(:,:,:,:,SliI);
    kdata_GE=permute(kdata_GE,[2 3 1 4]);
    kdata_GE=gflip(kdata_GE,2);
    
    load([CalibP,filename_calib,'_Slice',num2str(SliI),'.mat'],'kdata_calib');
    kdata_calib=double(kdata_calib);
    
    Sz=gsize(kdata_GE,1:2);
    nChannels=size(kdata_GE,4);
    
    kdata_calibP=permute(kdata_calib,[2 3 1 4]);

    Data=squeeze(sum(kdata_GE(:,:,1:30,:),3));
    Datab=squeeze(sum(kdata_GE(:,:,30+(1:30),:),3));
    Calib=squeeze(kdata_calibP(:,:,1,:));
    
%     res = GRAPPA(Data,Calib,[5 5],0.01);
%     resS(:,:,:,SliI)=res;
    resb = GRAPPA(Datab,Calib,[5 5],0.01);
    resbS(:,:,:,SliI)=resb;
end
disp('Finished GRAPPA');
iresS=ifft2cg(resS);
iresbS=ifft2cg(resbS);
% fgmontage(grmss(iresS,3))
%%
save([OutP 'resGRAPPASquashed.mat'],'resS','resbS');
%%
load([OutP 'resGRAPPASquashed.mat'],'resS','resbS');
%%
%% Run on all slices
slice=4;
for slice=[1:nSlices]
    disp(slice);
    im_recon_forSens=squeeze(iresS(:,:,:,slice));
    % disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
    CurSens=RunESPIRiTForSensMapsMultiMap(im_recon_forSens,0);
    CurSens=CurSens(:,:,:,1);
    SensMsk=grmss(CurSens,3)>0.01;
    CurSensSG(:,:,:,slice)=CurSens;
    SensMskSG(:,:,slice)=SensMsk;
    CurSens=CurSensSG(:,:,:,slice);
end
save([OutP 'Sens.mat'],'CurSensSG','SensMskSG');
%%
Echo1=sum(iresS.*conj(CurSensSG),3);
Echo2=sum(iresbS.*conj(CurSensSG),3);

    %
%     RecAllEchosChannels=ifft2c(recon);
%     RecAllEchos=sum(RecAllEchosChannels.*conj(perm43(CurSens)),4);
%     RecAllEchosS(:,:,:,slice)=RecAllEchos;
    %
%     WhichEchosToUse=15:55;
%     % WhichEchosToUse=10:14;
%     [PDBaseSG(:,:,slice), UpdatedB0Map_HzSG(:,:,slice), UpdatedT2SMap_msSG(:,:,slice), s_valsSG(:,:,:,slice),...
%         FittedSG(:,:,:,slice), PDBase0SG(:,:,slice)]=...
%         FitToModel_MPBD1CSf(RecAllEchosS(:,:,:,slice),WhichEchosToUse,ES_ms,FirstTE_ms);
% end

%%
directory_rawdata_fully = '/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
FNBaseFully='meas_MID00868_FID32103_ep2d_ge_sms1_EPTI_1p9_fully';
load([directory_rawdata_fully FNBaseFully filesep 'kdata_Slice' num2str(4) '.mat'],'kdata');
idata4=ifft2cg(kdata);

load([directory_rawdata_fully FNBaseFully filesep 'kdata_Slice' num2str(11) '.mat'],'kdata');
idata11=ifft2cg(kdata);

save([directory_rawdata_fully FNBaseFully filesep 'Recs411.mat'],'idata4','idata11');
fgmontagex(grmss(idata(:,:,[1,15,45],:),4))
GRecs=cat(3,grmss(iresS(:,:,:,4),3),grmss(iresbS(:,:,:,4),3));
FRecs=grmss(idata4(:,1:120,[15 45],:),4);
GFRecs=cat(4,GRecs,FRecs);
GFRecs=GFRecs./grms(GFRecs,1:3);

fgmontagex(GFRecs,[0 3]);
Ttls={'GRAPPA 1shot','Fully sampled'};
AddTopTtls(Ttls,[2 2]);

%%

%%
CS_Dim=5;
Ch_Dim=4;
TS_Dim=7;

CS_Flag=2^(CS_Dim-1);
Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);
% kdata_GE
%
% kdata_calib
%%
SliI=11;
for SliI=1:nSlices
    disp(SliI);
    kdata_GE=kdata_GES(:,:,:,:,SliI);
    kdata_GE=permute(kdata_GE,[2 3 1 4]);
    kdata_GE=gflip(kdata_GE,2);
    
    load([CalibP,filename_calib,'_Slice',num2str(SliI),'.mat'],'kdata_calib');
    kdata_calib=double(kdata_calib);
    
    Sz=gsize(kdata_GE,1:2);
    nChannels=size(kdata_GE,4);
    
    kdata_calibP=permute(kdata_calib,[2 3 1 4]);
    idata_calibP=ifft2cg(kdata_calibP);
    
    CurSens=RunESPIRiTForSensMapsMultiMap(squeeze(idata_calibP(:,:,1,:)),0);
    CurSens=CurSens(:,:,:,1);
    SensMsk=grmss(CurSens,3)>0.01;
    CurSensA=CurSens;
    CurSens=imresize(CurSens,Sz);
    CalibCCombined=sum(idata_calibP.*conj(perm43(CurSensA)),4);
    CalibCCombinedr=imresize(CalibCCombined,Sz);
    ES_ms=parameters.iEffectiveEpiEchoSpacing/1e3;
    FirstTE_ms=9;
    WhichEchosToUse=1:23;
    [PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
        FitToModel_MPBD1CSf(CalibCCombined,WhichEchosToUse,ES_ms,FirstTE_ms);
    
    B0r=imresize(UpdatedB0Map_Hz,Sz);
    B0rS(:,:,SliI)=B0r;
    
    
    Ch2D=reshape(kdata_GE,prod(gsize(kdata_GE,1:3)),nChannels);
    [~,S,sccmtx] = svd(Ch2D(1:1:end,:),'econ');
    sccmtxS(:,:,SliI)=sccmtx;
    
    ncc=31;
    kCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(kdata_GE,[4 1 2 3])),[4 1 2 3]);
    
    sig=permute(kCC,[1 2 5 4 6 7 3]);
    sigS(:,:,SliI,:,:,:,:)=sig;
    
    CurSensCCC=ipermute(MultMatTensor(sccmtx(:,1:ncc).',permute(CurSens,[3 1 2])),[3 1 2]);
    SensCSMap=perm43(CurSensCCC);
    SensCSMapS(:,:,SliI,:)=SensCSMap;
    
    mask_sample=abs(sig(:,:,:,1,:,:,:))>0;
    nEchos=size(sig,TS_Dim);
    nCS=1;
    clear kCC
    save([OutP 'SensCSMapS.mat'],'SensCSMapS','sccmtxS');
    MskCC=perm75(repmat(mask_sample,[1 1 1 ncc 1 1 1 1 1]));
    
    
    tmp=sigS(:,:,SliI,:,:,:,:,:,:,:);
    SigMskd(:,SliI)=single(tmp(MskCC));
end
disp('Finished preparation');
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
    
    
    B0_Hz=B0rS(:,:,SliI);
%     B0_Hz=UpdatedB0Map_HzSG(:,:,SliI);
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
    St.nNeighbors=70;
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
    St.SessionNameBase=[OutP 'MLN/S' num2str(SliI)];
    St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
    St.BaseTSDataP=BaseOutDir;
    St.BaseNUFTDataP=BaseOutDir;
    % if(strcmp(HostName,'tiger'))
    % end
    % disp(TSstr)
    % ParamsOutFn=['/autofs/cluster/kawin/Gilad/TF/ParamsForMES' num2str(SliI) '.txt'];
    ParamsOutFn=[BaseBaseOutP 'ParamsForMES' num2str(SliI) '.txt'];
    Txt=gStruct2txt(St,ParamsOutFn);
    disp('ok')
    MSLines{SliI}=['/usr/pubsw/bin/python3 /autofs/cluster/kawin/Gilad/TF/srezN/srez_main1.py ' ParamsOutFn];
end
%% End MLN
