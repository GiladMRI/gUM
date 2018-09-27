MLN=rot90(loadniidata('MLN_mc.nii'));
MLNm=mean(MLN,4);
MLN_E=MLN(:,:,:,2:2:end);
MLN_O=MLN(:,:,:,1:2:end);
MLN_Perf=MLN_E-MLN_O;
MLN_Perfm=mean(MLN_Perf,4);
MLN_Perfs=std(MLN_Perf,0,4);
MLN_Perf_t=MLN_Perfm./MLN_Perfs;

MLN_Msk=rot90(loadniidata('MLN_ref_brain_mask_hdrfix.nii'));
% MLN_tSNR=
ESP=rot90(loadniidata('resL1ESPIRiTCCS1A_mc.nii'));
ESPm=mean(ESP,4);
ESPm=ESPm.*grmss(MLNm)./grmss(ESPm);
ESP_E=ESP(:,:,:,2:2:end);
ESP_O=ESP(:,:,:,1:2:end);
ESP_Perf=ESP_E-ESP_O;
ESP_Perfm=mean(ESP_Perf,4);
ESP_Perfs=std(ESP_Perf,0,4);
ESP_Perf_t=ESP_Perfm./ESP_Perfs;
ESP_Msk=rot90(loadniidata('resL1ESPIRiTCCS1A_ref_brain_mask_hdrfix.nii'));
%%
Fore=MLN_Perf_t.*MLN_Msk;
Back=MLNm;
PStr='MLN';

% Fore=ESP_Perf_t.*ESP_Msk;
% Back=ESPm;
% PStr='ESP';

ClimF=[0.05 1.8];
ClimB=[0 1];
BThreshForF=0.2;
WhichSlicesShow=1:12;

Fore(Back<BThreshForF)=0;
A=(min(max(Fore,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
for s=1:numel(WhichSlicesShow)
    RGB(:,:,:,s)=ind2rgb(round(squeeze(A(:,:,WhichSlicesShow(s)))*255)+1,hot(256));
end
Msk=permute(Fore(:,:,WhichSlicesShow)>ClimF(1),[1 2 4 3]);

B=(min(max(Back,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+permute(B(:,:,WhichSlicesShow),[1 2 4 3]).*(1-Msk);

% figure;imshow(X)
%
Z=PartitionDim(X,4,3);
Z=CombineDims(Z,[4 2]);
Z=CombineDims(Z,[4 1]);
figure;imshow(Z);title(PStr)