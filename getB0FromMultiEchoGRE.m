BaseP='/media/a/DATA/2018_01_25/';

FN='meas_MID135_BP_gre_Sri_2Dsingleslice_FatSat_FID279';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg');
Data=sTwix.image();
%%
D=gpermute(Data,[3 1 5 2 8 9 4 6 7]);
D=double(sum(D,6));
nSli=size(D,3);
% Ord=[2:2:nSli 1:2:nSli];
Ord=[1:2:nSli 2:2:nSli];
[~,ROrd]=sort(Ord);
I=fft2cg(D(:,:,ROrd,:,:));
I=flip(I,2);
nEchos=size(I,5);

%%
BadChannelsI=[6 11 29];
NCha=sTwix.image.NCha;
GoodChannelsI=setdiff(1:NCha,BadChannelsI);

setenv('TOOLBOX_PATH','/home/a/bart-gpu_tensor_NUFMAC')
for s=1:nSli
    Sens(:,:,:,s)=RunESPIRiTForSensMaps(squeeze(I(:,:,s,GoodChannelsI,1)),30);
end
for s=1:nSli
    for e=1:nEchos
        I1(:,:,s,e)=CalcSENSE1f(squeeze(I(:,:,s,GoodChannelsI,e)),Sens(:,:,:,s));
    end
end
%% No FATSAT part
FN_NoFS='meas_MID137_BP_gre_Sri_2Dsingleslice_NoFatSat_FID281';


sTwix_NoFS = mapVBVD([BaseP FN_NoFS '.dat'],'removeOS','ignoreSeg');
Data_NoFS=sTwix_NoFS.image();
%
D_NoFS=gpermute(Data_NoFS,[3 1 5 2 8 9]);
D_NoFS=double(sum(D_NoFS,6));
I_NoFS=fft2cg(D_NoFS(:,:,ROrd,:,:));
%
for s=1:nSli
    for e=1:nEchos
        I1_NoFS(:,:,s,e)=CalcSENSE1f(squeeze(I_NoFS(:,:,s,:,e)),Sens(:,:,:,s));
    end
end
%
Both=cat(4,I1(:,:,:,1),I1_NoFS(:,:,:,1));
DBoth=abs(Both(:,:,:,2))-abs(Both(:,:,:,1));
%%
save([BaseP FN '.mat'],'I1','sTwix','Sens');

% Xs=80:187;
Xs=1:size(I1,2);
Ix=I1(:,Xs,:,:);
Ix1=Ix(:,:,:,1);
%%

WhichTwo=[1 2];
% WhichTwo=[2 3];
M=Ix(:,:,:,WhichTwo(1))./Ix(:,:,:,WhichTwo(2));
Mag=grmss(Ix,4);
Msk=imfillholesBySlices(Mag>2e-4);

Combined=sum(Mag.*exp(1i*angle(M)),4);
gammaMHz=42.5774806;
gammaHz=gammaMHz*1e6;

TEs_us=[sTwix.hdr.Phoenix.alTE{:}];
deltaTE_us=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));
scanFreq_Hz=sTwix.hdr.Config.ScanFrequency; % SystemFrequency

dAngle=double(angle(Combined.*exp(1i*2*pi*0*scanFreq_Hz*deltaTE_us/1e6)));
B0_Hz=dAngle/(2*pi*deltaTE_us/1e6);
B0_Hz(~Msk)=NaN;
fgmontage(B0_Hz);colorbar
%%
[PhiCostantini] = cunwrap(dAngle(:,:,s), struct('weight',Mag(:,:,s),'RoundK',false,'maxblocksize',60));

B0_HzU=PhiCostantini/(2*pi*deltaTE_us/1e6);

fgmontage(PhiCostantini)