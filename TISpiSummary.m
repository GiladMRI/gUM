ScanP='/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/';
BaseFN='meas_MID01088_FID09952_gSpi2d_T12_d110_Dw11_40rep_VD123';
RefFldMapP=[ScanP 'meas_MID01104_FID09968_gre_te7_40' filesep];

BaseFNs={'meas_MID01090_FID09954_gSpi2d_T13_d110_Dw11_40rep_VD123',...
    'meas_MID01092_FID09956_gSpi2d_T14_d110_Dw11_40rep_VD123',...
    'meas_MID01094_FID09958_gSpi2d_T15_d110_Dw11_40rep_VD123',...
    'meas_MID01102_FID09966_gSpi2d_T12_d110_Dw22_40rep_VD123',...
    'meas_MID01106_FID09970_gSpi2d_T12_d110_Dw11_40rep_VD123_Thin',...
    'meas_MID01108_FID09972_gSpi2d_T14_d110_Dw11_40rep_VD123_Thin',...
    'meas_MID01110_FID09974_gSpi2d_T15_d110_Dw11_40rep_VD123_Thin'};

for fIdx=1:numel(BaseFNs)
    BaseFN=BaseFNs{fIdx};
%%
mainP=[ScanP BaseFN];
%%
clear DataForSlice
load([mainP filesep 'DataForSlice.mat']);
%%
nSlices=10;
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

TrajType=str2num(BaseFN(32:33));
ResType=floor((TrajType-10)/2)+1; % 2,1.5,1
TimingType=mod(TrajType-10,2)+1; % 5ms, 10ms
load('GAll5ms10ms.mat');
GTrajaCBase=GAll(:,TimingType,ResType);
if(TimingType==1)
    nInnerShots=10;
else
    nInnerShots=5;
end
ResStrs={'xxx','1.9mm','1.3mm'};
InnerShotLengthStrs={'5ms','10ms'};

ResStr=ResStrs{ResType};
InnerShotLengthStr=InnerShotLengthStrs{TimingType};

ThinStr='';
if(strcmp(BaseFN(end-3:end),'Thin'))
    ThinStr='Thin_';
end
FNBase=[ScanP BaseFN(6:14) ResStr '_' InnerShotLengthStr '_' ThinStr];
%%
ncc=31;

clear UpdatedB0MapS UpdatedB0Map0S UpdatedT2SMap_msS UpdatedT2SMap_ms0S Rec1ccpMtS Rec1ccpMtInS Rec1ccpMtOutAndInS THLRMultiShotS
% SliI=2;
for SliI=1:nSlices
    UpdatedB0Map=DataForSlice{SliI}{1};
    UpdatedT2SMap_ms=DataForSlice{SliI}{2};
    UpdatedB0Map0=DataForSlice{SliI}{3};
    UpdatedT2SMap_ms0=DataForSlice{SliI}{4};
    THLRMultiShot=DataForSlice{SliI}{5};
    Rec1ccpMtOutIn=DataForSlice{SliI}{6};
    Rec1ccpMt=DataForSlice{SliI}{7};
    Rec1ccpMtIn=DataForSlice{SliI}{8};
    SelfSens1=DataForSlice{SliI}{9};
    sccmtx=DataForSlice{SliI}{10};

    SensCC=permute(sum(SelfSens1.*permute(sccmtx(:,1:ncc),[3 4 1 5 6 7 8 9 2]),3),[1:2 9 3:8]);
    SensCC=permute43(SensCC);
    disp('ok SensCC');


    Rec1ccpMtOutAndIn=cat(4,Rec1ccpMt,Rec1ccpMtIn);
%     Rec1ccpMtOutAndIn=CombineDims(Rec1ccpMtOutAndIn,[3 4]);
    UpdatedB0MapS(:,:,SliI)=UpdatedB0Map;
    UpdatedB0Map0S(:,:,SliI)=UpdatedB0Map0;
    UpdatedT2SMap_msS(:,:,SliI)=UpdatedT2SMap_ms;
    UpdatedT2SMap_ms0S(:,:,SliI)=UpdatedT2SMap_ms0;
    Rec1ccpMtS(:,:,:,SliI)=Rec1ccpMt;
    Rec1ccpMtInS(:,:,:,SliI)=Rec1ccpMtIn;
    Rec1ccpMtOutAndInS(:,:,:,:,SliI)=Rec1ccpMtOutAndIn;
    THLRMultiShotS(:,:,:,SliI)=squeeze(THLRMultiShot);
end
%%
Rec1ccpMtOutAndInSN=Rec1ccpMtOutAndInS./grms(Rec1ccpMtOutAndInS,[1 2 3]);
Rec1ccpMtOutAndInSN=CombineDims(Rec1ccpMtOutAndInSN,[3 4]);
OutAndInGifFN=[FNBase 'SpiralOutAndIn.gif'];
delete(OutAndInGifFN);
figure;pause(4);gmontage(Rec1ccpMtOutAndInSN(:,:,1,ROrd),[0 max(abs(Rec1ccpMtOutAndInSN(:)))/2]);
pause(4);gif(OutAndInGifFN);pause(4);
for i=2:size(Rec1ccpMtOutAndInSN,3)
    gmontage(Rec1ccpMtOutAndInSN(:,:,i,ROrd),[0 max(abs(Rec1ccpMtOutAndInSN(:)))/2]);
    gif
end
disp(OutAndInGifFN);
%%
Rec1ccpMtOutAndInSN=Rec1ccpMtOutAndInS(:,:,:,:,ROrd)./grms(Rec1ccpMtOutAndInS(:,:,:,:,ROrd),[1 2 3]);
Rec1ccpMtOutAndInSN=CombineDims(Rec1ccpMtOutAndInSN,[5 4]);
OutAndInGifFN=[FNBase 'SpiralOutAndIn_SideBySide.gif'];
delete(OutAndInGifFN);
figure;gmontage(Rec1ccpMtOutAndInSN(:,:,1,:),[0 max(abs(Rec1ccpMtOutAndInSN(:)))/2]);
pause(4);gif(OutAndInGifFN);pause(4);
for i=2:size(Rec1ccpMtOutAndInSN,3)
    gmontage(Rec1ccpMtOutAndInSN(:,:,i,:),[0 max(abs(Rec1ccpMtOutAndInSN(:)))/2]);
    gif
end
disp(OutAndInGifFN);
%% Single image
fgmontage(Rec1ccpMtOutAndInSN(:,:,3,:));
gprint([FNBase 'SingleInnerShot_MultiShot.png']);
%%
fgmontage(UpdatedB0Map0S(:,:,ROrd),[-100 100]);
gprint([FNBase 'B0FromSpiralOuts.png']);
fgmontage(UpdatedT2SMap_ms0S(:,:,ROrd),[0 100]);
gprint([FNBase 'T2SFromSpiralOuts.png']);
fgmontage(UpdatedB0MapS(:,:,ROrd),[-100 100]);
gprint([FNBase 'B0FromTHLRMS.png']);
fgmontage(UpdatedT2SMap_msS(:,:,ROrd),[0 100]);
gprint([FNBase 'T2SFromTHLRMS.png']);
%%
THLRMultiShotSN=THLRMultiShotS./grms(THLRMultiShotS,[1 2 3]);
fgmontage(THLRMultiShotSN(:,:,:,ROrd))
gprint([FNBase 'THLRMS.png']);
%%
close all
end