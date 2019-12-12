nSlices=16;

SetName='68xAx';
% SetName='68xCor';
if(strcmp(SetName,'68xAx'))
    EPTIP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns/';
    SpiP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00860_FID32095_gSpi2d_T10_Dw11_d110_VD1/';
    Sz=[114 114];
    
    MLNRXS1=CollectMLNResultsf('MID00860','ME_S',nSlices);
    MLNRXS3=CollectMLNResultsf('MID00860','ME_3S',nSlices);
end
if(strcmp(SetName,'68xCor'))
    EPTIP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00903_FID32138_ep2d_ge_EPTI_1p9_3shot_4dyns_Cor/';
    SpiP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00884_FID32119_gSpi2d_T12_Dw11_d110_VD1_Cor/';
    Sz=[116 116];
    
    MLNRXS1=CollectMLNResultsf('MID00884','ME_S',nSlices);
    MLNRXS3=CollectMLNResultsf('MID00884','ME_3S',nSlices);
end

load([EPTIP 'RecAllEchosS.mat'],'RecAllEchosS');
RecAllEchosS=gflip(RecAllEchosS,2);


OrdA=[2:2:nSlices 1:2:nSlices];
[~,ROrdA]=sort(OrdA);

ESz=[120 120];
EPTIRefS=zeros([ESz,80,nSlices]);
EMapsX=zeros([ESz,4,nSlices]);
for SliI=1:nSlices
    try
        EPTIRefPrefix=[EPTIP 'Sli' num2str(SliI) '_'];
        EPTIRefElemL3=readcfl([EPTIRefPrefix 'ElemsL_3']);
        
        EPTIRefPrefix=[EPTIP 'Sli' num2str(SliI) '_dyn2_'];
        EPTIRefElem0=readcfl([EPTIRefPrefix 'Elem0']);
        EPTIRefElem1=readcfl([EPTIRefPrefix 'Elem1']);
        EPTIRefElem2=readcfl([EPTIRefPrefix 'Elem2']);
        EPTIRefElem3=readcfl([EPTIRefPrefix 'Elem3']);
        
        EMapsX(:,:,1,SliI)=EPTIRefElem0;
        EMapsX(:,:,2,SliI)=EPTIRefElem1;
        EMapsX(:,:,3,SliI)=EPTIRefElem2;
        EMapsX(:,:,4,SliI)=EPTIRefElem3;
        
        EPTIRefS(:,:,:,SliI)=gflip(squeeze(EPTIRefElem0.*exp(EPTIRefElemL3./EPTIRefElem3)),2);
        disp(SliI);
    catch
    end
end
EMapsX=gflip(EMapsX,2);

RepSets={1, 1:2, 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:36};

RepSets={1, [1 13 25], 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:36};

SpiRefSR=zeros([Sz,21,nSlices,numel(RepSets)]);
MapsX=zeros([Sz,4,nSlices,numel(RepSets)]);
for rs=1:numel(RepSets)
    WhichRepsToUse=RepSets{rs};
    RepsStr=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((WhichRepsToUse).'),1),' ','')),'[','A'),']','A');
    
    for SliI=1:nSlices
        try
            SpiRefPrefix=[SpiP 'Sli' num2str(SliI) '_R' RepsStr '_'];
            SpiRefPrefix=[SpiP 'gB0x' 'Sli' num2str(SliI) '_R' RepsStr '_'];
            SpiRefElemL3=readcfl([SpiRefPrefix 'ElemsL_3']);
            SpiRefPrefix=[SpiP 'gB0x' 'Sli' num2str(SliI) '_R' RepsStr '_Lm45b_'];
            SpiRefElem0=readcfl([SpiRefPrefix 'Elem0']);
            SpiRefElem1=readcfl([SpiRefPrefix 'Elem1']);
            SpiRefElem2=readcfl([SpiRefPrefix 'Elem2']);
            SpiRefElem3=readcfl([SpiRefPrefix 'Elem3']);
            SpiRefSR(:,:,:,OrdA(SliI),rs)=squeeze(SpiRefElem0.*exp(SpiRefElemL3./SpiRefElem3));
            MapsX(:,:,1,OrdA(SliI),rs)=SpiRefElem0;
            MapsX(:,:,2,OrdA(SliI),rs)=SpiRefElem1;
            MapsX(:,:,3,OrdA(SliI),rs)=SpiRefElem2;
            MapsX(:,:,4,OrdA(SliI),rs)=SpiRefElem3;
            disp([rs SliI]);
        catch
        end
    end
end

RepSetsX={1, [1 13 25], 1:3, 1:4, 1:5, 1:6 1:7 1:8, 1:9, 1:15, 1:MaxRepsToUse};
rs=2;
RepsStrX=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((RepSetsX{rs}).'),1),' ','')),'[','A'),']','A');
for SliI=1:nSlices
    for rs=1:2
        RepsStrX=strrep(strrep(GroupToStr( strrep(gmat2cell(num2str((RepSetsX{rs}).'),1),' ','')),'[','A'),']','A');
        try
            CurSRPrefix=[mainP filesep 'gB0Sli' num2str(SliI) '_R' RepsStrX '_'];
            
            SpiRefElemL3gB0=readcfl([CurSRPrefix 'ElemsL_3']);
            SpiRefElem0gB0=readcfl([CurSRPrefix 'Elem0']);
            SpiRefElem1gB0=readcfl([CurSRPrefix 'Elem1']);
            SpiRefElem2gB0=readcfl([CurSRPrefix 'Elem2']);
            SpiRefElem3gB0=readcfl([CurSRPrefix 'Elem3']);
            SpiRefSRgB0(:,:,:,OrdA(SliI),rs)=squeeze(SpiRefElem0gB0.*exp(SpiRefElemL3gB0./SpiRefElem3gB0));
            MapsXgB0(:,:,1,OrdA(SliI),rs)=SpiRefElem0gB0;
            MapsXgB0(:,:,2,OrdA(SliI),rs)=SpiRefElem1gB0;
            MapsXgB0(:,:,3,OrdA(SliI),rs)=SpiRefElem2gB0;
            MapsXgB0(:,:,4,OrdA(SliI),rs)=SpiRefElem3gB0;
            disp([rs SliI]);
        catch
        end
    end
end
%%
EPTIRefS
%%
if(strcmp(SetName,'68xCor'))
    SpiRefSRx=gflip(padarray(SpiRefSR,[2 2],'Both'),1);
    RecAllEchosS=gflip(RecAllEchosS,1);
    EPTIRefS=gflip(EPTIRefS,1);
    
    MapsX=gflip(padarray(MapsX,[2 2],'Both'),1);
    EMapsX=gflip(EMapsX,1);
    
    MLNRXS3=gflip(double(padarray(MLNRXS3,(gsize(EPTIRefS,1:2)-gsize(MLNRXS3,1:2))/2,'Both')),1);
    MLNRXS1=gflip(double(padarray(MLNRXS1,(gsize(EPTIRefS,1:2)-gsize(MLNRXS1,1:2))/2,'Both')),1);
end

if(strcmp(SetName,'68xAx'))
    SpiRefSRx=gflip(padarray(SpiRefSR,[3 3],'Both'),[]);
    MapsX=padarray(MapsX,[3 3],'Both');
    
    MLNRXS3=double(padarray(MLNRXS3,(gsize(EPTIRefS,1:2)-gsize(MLNRXS3,1:2))/2,'Both'));
    MLNRXS1=double(padarray(MLNRXS1,(gsize(EPTIRefS,1:2)-gsize(MLNRXS1,1:2))/2,'Both'));
end

WhichTSToUs_MLN=1:8;
for i=1:nSlices
    disp(i);
    [~, UpdatedB0Map_MLN3(:,:,i), UpdatedT2SMap_msMLN3(:,:,i), s_valsMLN3(:,:,:,i), ~, PDBase0MLN3(:,:,i)]=...
        FitToModel_MPBD1CSf(squeeze(MLNRXS3(:,:,:,i)),WhichTSToUs_MLN,50/8,0+2.38);
    [~, UpdatedB0Map_MLN1(:,:,i), UpdatedT2SMap_msMLN1(:,:,i), s_valsMLN1(:,:,:,i), ~, PDBase0MLN1(:,:,i)]=...
        FitToModel_MPBD1CSf(squeeze(MLNRXS1(:,:,:,i)),WhichTSToUs_MLN,50/8,0+2.38);
end
disp('ok');
%%
SpiRefSRx=padarray(SpiRefSRgB0,[3 3],'Both');
SpiRefSRx=padarray(SpiRefSR,[3 3],'Both');
SliI=4;

fgmontagex(SpiRefSR(:,:,2:5:end,SliI,3));title('Spi splitProx');
fgmontagex(RecAllEchosS(:,:,9:17:end-5,SliI));title('EPTI GRAPPA');
fgmontagex(EPTIRefS(:,:,9:17:end-5,SliI));title('EPTI splitProx');

AllRecs=cat(4,RecAllEchosS(:,:,9:17:end-5,SliI),EPTIRefS(:,:,9:17:end-5,SliI),SpiRefSRx(:,:,2:5:end,SliI,3),SpiRefSRx(:,:,2:5:end,SliI,5));
AllRecs=AllRecs./grms(AllRecs,1:3);

AllRecs=cat(4,RecAllEchosS(:,:,9:17:end-5,SliI),EPTIRefS(:,:,9:17:end-5,SliI),SpiRefSRx(:,:,2:5:end,SliI,1),SpiRefSRx(:,:,2:5:end,SliI,2));
AllRecs=AllRecs./grms(AllRecs,1:3);
%%
AllRecs=cat(4,RecAllEchosS(:,:,9:17:end-5,SliI),EPTIRefS(:,:,9:17:end-5,SliI),SpiRefSRx(:,:,2:5:end,SliI,3),SpiRefSRx(:,:,2:5:end,SliI,5),...
    MLNRXS3(:,:,2:2:end,SliI),MLNRXS1(:,:,2:2:end,SliI));
AllRecs=AllRecs./grms(AllRecs,1:3);

fgmontagex(AllRecs(:,:,[1 3],:));caxis(caxis/1.5);
%%
SliI=9;
AllRecs=cat(4,RecAllEchosS(:,:,9:17:end-5,SliI),EPTIRefS(:,:,9:17:end-5,SliI),SpiRefSRx(:,:,2:5:end,SliI,3),SpiRefSRx(:,:,2:5:end,SliI,5));
% AllRecs=cat(4,RecAllEchosS(:,:,9:17:end-5,SliI),SpiRefSRx(:,:,2:5:end,SliI,3),SpiRefSRx(:,:,2:5:end,SliI,5));
AllRecs=AllRecs./grms(AllRecs,1:3);

fgmontagex(AllRecs(:,:,[1 3],:));Ttls={'EPTI3 GRAPPA','EPTI3 SplitProx','Spi3 SplitProx','Spi5 SplitProx'};
MSz=size(get(get(gca,'Children'),'CData'));
for i=1:numel(Ttls)
    h=text(MSz(2)/4*i-MSz(2)/4/2,7,Ttls{i},'FontSize',20,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',0);
end
ylabel(['Slice #' num2str(SliI)]);caxis(caxis/1.5);

fgmontagex(SpiRefSRx(:,:,15,SliI,:));title('Spirals 1     2     3     4     5     6     7     8     9    15    36');ylabel(['Slice #' num2str(SliI)]);
%%
fgmontagex(SpiRefSRx(:,:,15,:,2));
AllRecs2=perm43(cat(3,EPTIRefS(:,:,40,:),SpiRefSRx(:,:,10,:,3)));

AllRecs2=perm43(cat(3,EPTIRefS(:,:,40,:),SpiRefSRx(:,:,10,:,3),MLNRXS3(:,:,4,:),MLNRXS1(:,:,4,:)));
AllRecs2=AllRecs2./grms(AllRecs2,1:3);
% fgmontagex(perm43(AllRecs2(:,:,[1 4 7 10 14 16],:)));caxis(caxis/1.5);
fgmontagex(perm43(AllRecs2(:,:,[1 4 7 10 14 15],:)));caxis(caxis/1.5);
MSz=size(get(get(gca,'Children'),'CData'));
Ttls={'EPTI3 splitProx','Spiral3 splitProx','Spiral3 MLN','Spiral1 MLN'};
for i=1:numel(Ttls)
    h=text(7,MSz(1)/4*i-MSz(1)/4/2,Ttls{i},'FontSize',14,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',90);
end

h=text(7,60,'EPTI-3','FontSize',20,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',90);
h=text(7,180,'Spiral-3','FontSize',20,'Color',[1 0 0.5],'HorizontalAlignment','center');set(h,'Rotation',90);

save([SetName '_Res.mat'],'SpiRefSRx','EPTIRefS','RecAllEchosS','MapsX','EMapsX');
%%
AllT=perm43(cat(3,EMapsX(:,:,4,:),MapsX(:,:,4,:,3)));
AllT=perm43(cat(3,EMapsX(:,:,4,:),MapsX(:,:,4,:,3),perm43(UpdatedT2SMap_msMLN3),perm43(UpdatedT2SMap_msMLN1)));
fgmontagex(perm43(AllT(:,:,[1 4 7 10 14 15],:)),[0 200]);colormap hot
%%
fgmontagex(SpiRefSx(:,:,2:5:end,5)); title('Spi splitProx');
fgmontagex(gflip(RecAllEchosS(:,:,9:17:end-5,5),2)); title('EPTI GRAPPA');
fgmontagex(gflip(EPTIRefS(:,:,9:17:end-5,5),2)); title('EPTI splitProx');

fgmontagex(SpiRefS(:,:,12,:))
