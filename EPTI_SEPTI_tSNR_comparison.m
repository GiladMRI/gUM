mainP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00864_FID32099_gSpi2d_T10_Dw11_d120_VD1';
load([mainP filesep 'aRecMXSR.mat']);
aSEPTI=aRecMXSR;
aSEPTI=aSEPTI(:,:,:,ROrd,:);

OutP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/meas_MID00876_FID32111_ep2d_ge_EPTI_1p9_3shot_4dyns/';
load([OutP 'aRecMXDS.mat']);
aEPTI=aRecMXDS;
aEPTI=gflip(aEPTI,2);
aEPTI=aEPTI(4:end-3,4:end-3,:,:,:);
%%
nEchos_EPTI=size(aEPTI,3);
nEchos_SEPTI=size(aSEPTI,3);

EchoForMsk_EPTI=80;
EchoForMsk_SEPTI=21;

maEPTI=mean(aEPTI,5);
saEPTI=std(aEPTI,[],5);
taEPTI=maEPTI./saEPTI;

se = strel('disk',5);

ForMaskE=squeeze(mean(maEPTI,3));
ForMaskE=ForMaskE./grmss(ForMaskE);
[~,MskEPTI,MskNEPTI]=CalcSlicesSNR(ForMaskE,false,9);
% MskNEPTI=imfillholesBySlices(~MskNEPTI)>0.5;
MskNEPTI=imfillholesBySlices(MskEPTI)>0.5;
for i=1:nSlices
    MskNSEPTI(:,:,i)=imopen(MskNSEPTI(:,:,i),se);
    MskNEPTI(:,:,i)=getLargestComponent(MskNEPTI(:,:,i));
end

MskNEPTI=imfillholesBySlices(getLargestComponent(ForMaskE>0.7))>0.5;

maSEPTI=mean(aSEPTI,5);
maSEPTI(:,:,:,9)=mean(aSEPTI(:,:,:,9,[1:5 7:12]),5);
saSEPTI=std(aSEPTI,[],5);
saSEPTI(:,:,:,9)=std(aSEPTI(:,:,:,9,[1:5 7:12]),[],5);
taSEPTI=maSEPTI./saSEPTI;

ForMaskS=squeeze(mean(maSEPTI,3));
ForMaskS=ForMaskS./grmss(ForMaskS);
[~,MskSEPTI,MskNSEPTI]=CalcSlicesSNR(ForMaskS,false,9);
% MskNSEPTI=imfillholesBySlices(~MskNSEPTI)>0.5;
MskNSEPTI=imfillholesBySlices(MskSEPTI)>0.5;
for i=1:nSlices
    MskNSEPTI(:,:,i)=imopen(MskNSEPTI(:,:,i),se);
    MskNSEPTI(:,:,i)=getLargestComponent(MskNSEPTI(:,:,i));
end

MskNSEPTI=imfillholesBySlices(getLargestComponent(ForMaskS>0.7))>0.5;

disp('ok');
%%
fgmontagex(squeeze(maEPTI(:,:,EchoForMsk_EPTI,:)));title('EPTI');
fgmontagex(squeeze(maSEPTI(:,:,EchoForMsk_SEPTI,:)));title('SEPTI');

fgmontagex(MskNEPTI);
fgmontagex(MskNSEPTI);

fgmontagex((squeeze(maEPTI(:,:,floor(linspace(1,nEchos_EPTI,8)),8))));title('EPTI');caxis(caxis/2);
fgmontagex((squeeze(maSEPTI(:,:,floor(linspace(1,nEchos_SEPTI,8)),8))));title('SEPTI');caxis(caxis/2);
%%
BinEdges=0:5:200;
BinCenters=(BinEdges(1:end-1)+BinEdges(2:end))/2;

EchoFortSNR=[40 11];
figure;
tmp=squeeze(taEPTI(:,:,EchoFortSNR(1,1),:));
subplot(1,2,1);histogram(tmp(MskNEPTI),BinEdges);title('tSNR EPTI 3-shot')
tmp=squeeze(taSEPTI(:,:,EchoFortSNR(1,2),:));
subplot(1,2,2);histogram(tmp(MskNSEPTI),BinEdges);title('tSNR SKEPTIC 3-shot')
close all
%%
EchoFortSNR=[1 1; 40 11; 80 21];
figure;

for i=1:size(EchoFortSNR,1)
    tmp=squeeze(taEPTI(:,:,EchoFortSNR(i,1),:));
    he=histcounts(tmp(MskNEPTI),BinEdges);
    he=he./sum(he);
    tmp=squeeze(taSEPTI(:,:,EchoFortSNR(i,2),:));
    hs=histcounts(tmp(MskNSEPTI),BinEdges);
    hs=hs./sum(hs);
    subplot(1,size(EchoFortSNR,1),i);
    plot(BinCenters,he,'r');hold on;plot(BinCenters,hs,'k');
    title(['tSNR, Echo ' num2str(i) '/21']);
    if(i==1)
        legend({'EPTI','SKEPTIC'});
    end
end
%%
maEPTIx=maEPTI(:,:,floor(linspace(1,80,21)),:);
Both=cat(5,maEPTIx,maSEPTI);
Both=Both./grms(Both,1:3);

fgmontagex(perm43(squeeze(Both(:,:,[2 7 14 21],3,:))));ylabel([PadStringWithBlanks('SKEPTIC',40) 'EPTI'],'FontSize',20)
fgmontagex(perm43(squeeze(Both(:,:,[2 7 14 21],8,:))));ylabel([PadStringWithBlanks('SKEPTIC',40) 'EPTI'],'FontSize',20)
fgmontagex(perm43(squeeze(Both(:,:,[2 7 14 21],13,:))));ylabel([PadStringWithBlanks('SKEPTIC',40) 'EPTI'],'FontSize',20)