ScanP='/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/';
BaseFN='meas_MID01088_FID09952_gSpi2d_T12_d110_Dw11_40rep_VD123';
% BaseFN='meas_MID01106_FID09970_gSpi2d_T12_d110_Dw11_40rep_VD123_Thin';

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
FNBaseX=[BaseFN(6:14) ResStr '_' InnerShotLengthStr '_' ThinStr];
FNBase=[ScanP FNBaseX];


load([ScanP BaseFN filesep 'PerSlice_gB0.mat']);
% load([ScanP BaseFN filesep 'PerSliceRec.mat']);
nSlices=10;
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

for SliI=1:nSlices
    Rec_CompgB0_MX=squeeze(sum(Rec_CompgB0_MS{SliI}.*CompsPS{SliI},6));
    Rec_CompgB0_MXS(:,:,:,:,SliI)=Rec_CompgB0_MX;
end
%%
HankelTemporalLen=2;
%%
Rec_CompgB0_MXSN=Rec_CompgB0_MXS./grms(Rec_CompgB0_MXS,[1 2 3 4]);
% XX=grmss(Rec_CompgB0_MXSN(:,:,:,10:30,ROrd),4);
% XX=grmss(Rec_CompgB0_MXSN(:,:,:,22:end,ROrd),4);
XX=grmss(Rec_CompgB0_MXSN(:,:,:,10:end,ROrd),4);

% YY=Rec_CompgB0_MXSN(:,:,:,22:end,ROrd);
% nHpoints=size(YY,4);
% [~, ~, ~,H]=ghankel(nHpoints,HankelTemporalLen,gsize(YY,1:2));
% 
% s_LLR=zeros([gsize(YY,1:2) 2 size(YY,3) nSlices]);
% V_LLR=zeros([gsize(YY,1:2) 2 2 size(YY,3) nSlices]);
% for i=1:size(YY,3)
%     for s=1:nSlices
%         disp([i s]);
%         [ ~, s_LLR(:,:,:,i,s), V_LLR(:,:,:,:,i,s) ] = batch_svd(H*squeeze(YY(:,:,i,:,s)));
%     end
% end
% XX=squeeze(s_LLR(:,:,1,:,:));

sSig=std(XX,0,3);
mSig=mean(XX,3);
tSig=mSig./sSig;

Mx=max(mSig(:))/1.5;
fgmontage(mSig,[0 Mx]);title('rms over central echos, mean over shots');
gprint([ScanP 'SubspaceLLR_' FNBaseX 'meanOverShots']);
fgmontage(sSig,[0 Mx/10]);title('rms over central echos, std over shots, 10x');
gprint([ScanP 'SubspaceLLR_' FNBaseX 'stdOverShots']);
fgmontage(tSig);title('m/s');colorbar;
gprint([ScanP 'SubspaceLLR_' FNBaseX 'tSNROverShots']);

fgmontage(Rec_CompgB0_MXSN(:,:,13,:,ROrd));
gprint([ScanP 'SubspaceLLR_' FNBaseX 'SingleShot_AllEchos_AllSlices']);
fgmontage(Rec_CompgB0_MXSN(:,:,13,:,5));
gprint([ScanP 'SubspaceLLR_' FNBaseX 'SingleShot_AllEchos_1Slices']);
fgmontage(Rec_CompgB0_MXSN(:,:,13,:,1));
gprint([ScanP 'SubspaceLLR_' FNBaseX 'SingleShot_AllEchos_1Slicesx']);
fgmontage(XX(:,:,1:8:end,2:3:end));title('rms over central echos, some shots, some slices');
gprint([ScanP 'SubspaceLLR_' FNBaseX 'SingleShot_SomeEchos_SomeSlices']);
%%
close all
%%


% XX=grmss(Rec_CompgB0_MXSN(:,:,:,10:end,ROrd),4);
% fgmontage(Rec_CompgB0_MXSN(:,:,13,5:15,ROrd(5)))
% XX=grmss(Rec_CompgB0_MXSN(:,:,:,22:end,ROrd),4);
% fgmontage(XX(:,:,1:8:end,2:3:end));title('rms over central echos, some shots, some slices');
% 
% mSig=mean(XX,3);
% Mx=max(mSig(:))/1.5;
% 
% fgmontage(mSig,[0 Mx]);title('rms over central echos, mean over shots');

