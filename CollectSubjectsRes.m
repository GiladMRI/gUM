%% Collect MLN Results
WhichSubj=[1 2 3 4];
CurSubj=4;
for CurSubj=WhichSubj
    BaseP='/media/a/DATA/ASLSubjData/';
    CurSubjP=[BaseP 'S' num2str(CurSubj,'%02d') filesep];
    D=dir(CurSubjP);
    MeasFN=D(strhas({D.name},'_U19')).name;
    CurSubjMP=[CurSubjP MeasFN filesep];
    D=dir(CurSubjMP);
    D=D([D.isdir]);
    D=D(strhas({D.name},'train'));
    for SliI=1:12
%         SliI=5;
        CurSliP=[CurSubjMP D(strhas({D.name},['_Sli' num2str(SliI,'%02d') '__'])).name filesep];
        [ScrNC,BatchNC,MinNC,LastFN]=GraphOptFromFolderf(CurSliP);
        X=imread([CurSliP LastFN]);
        Y=X(1:128,128*5+(1:128),1);
        P=double(X(1:128,128*7+(1:128),1));
        Y=double(Y).*exp(1i*2*pi*P/256+1i*pi);
        
        MLNOut(:,:,SliI,CurSubj)=Y;
    end
    DM=dir([CurSubjMP '*.mat']);
    L1ESP=DM(strhas({DM.name},'^L1ESPI'));
    L1ESP=L1ESP(strhas({L1ESP.name},'_CC.mat'));
    X=load([CurSubjMP L1ESP.name]);
    resL1ESPIRiTCCS1(:,:,:,CurSubj)=X.resL1ESPIRiTCCS1;
    
    BARTR=DM(strhas({DM.name},'^BARTRecon_A'));
    X=load([CurSubjMP BARTR.name]);
    RecS(:,:,:,CurSubj)=X.RecS;
    
    BARTR2=DM(strhas({DM.name},'^BARTRecon2Ma'));
    X=load([CurSubjMP BARTR2.name]);
    RecSMM(:,:,:,CurSubj)=X.RecSMM;
    
    X=load([CurSubjMP 'im_resS.mat']);
    SPARSEMRI(:,:,:,CurSubj)=X.im_resS;
end
AllRec=cat(5,double(MLNOut),resL1ESPIRiTCCS1,SPARSEMRI,RecS,RecSMM);
for i=1:size(AllRec,3)
    for j=1:size(AllRec,4)
        for k=1:size(AllRec,5)
            CurI=AllRec(:,:,i,j,k);
            Msk=CurI~=0;
            AllRecN(:,:,i,j,k)=rot90(CurI/grmss(CurI(Msk)));
        end
    end
end
%%
% fgmontage(MLNOut)
% fgmontage(resL1ESPIRiTCCS1)
% fgmontage(RecS)
% fgmontage(RecSMM)
% fgmontage(SPARSEMRI)
% L1ESPIRiT_B0_lam1e-05_CC
% BARTRecon_AllS_NoB0_W1e-05
% BARTRecon2Maps_AllS_NoB0_W1e-05
% im_resS.mat