N=64;
EstSize=40;
SigCenterSize=8;

Sz=[N N];

nEchos=size(X,3);
X=imresizeBySlices(X,Sz);

FullSig=fft2cg(X);
% EstSig=crop(FullSig,EstSize,EstSize,nEchos);
% EstX=ifft2cg(EstSig);
% EstXPadded=imresizeBySlices(EstX,Sz)/(N/EstSize);
EstXPadded=SmoothBySlices(X,[20 20],2);

TotalOSRatio=1.5;
AccPerEcho=nEchos/TotalOSRatio;
sAccPerEcho=sqrt(AccPerEcho);
disp('ok sig');
%%
clear Msk
for i=1:nEchos
    disp(i);
    Msk(:,:,i)=squeeze(bart(['poisson -Y ' num2str(N) ' -Z ' num2str(N) ' -y ' num2str(sAccPerEcho) ' -z ' num2str(sAccPerEcho) ' -C ' num2str(SigCenterSize) ' -s ' num2str(rand*100000)]));
end

disp('ok Msk');
% Msk=permute(squeeze(mask_sample(1:nEchos

% Msk=Msk*0+1;
% Msk(:,:,5)=0;
%%
Sig=FullSig.*Msk;
disp(['Effective OS ratio: ' num2str(gsum(Msk)/prod(Sz))]);

% fgmontage(sum(Msk,3))
% fgmontage(sum(Msk,3)>0)
%% Now with sens
ManySens=load('/autofs/cluster/kawin/Gilad/CCSensMaps.mat');
Sens=perm43(imresizeBySlices(perm31(ManySens.SensCC(:,:,:,15)),Sz));
disp('Loaded sens');
%%
FullSig=fft2cg(X.*Sens);
Sig=FullSig.*Msk;
SensMsk=grmss(Sens,3:4)>0.05;
disp('OK sig with sens');
%%
X_Estr=imresizeBySlices(X_Est,Sz);

B0_Hzr=imresizeBySlices(B0_Hz,Sz);

B0_Hz_Estr=imresizeBySlices(B0_Hz_Est,Sz);

SB0Var_Hzr=imresizeBySlices(SB0Var_Hz,Sz);

T2S_msr=imresizeBySlices(T2S_ms,Sz);
R2Sr=1./T2S_msr;

T2S_ms_Estr=imresizeBySlices(T2S_ms_Est,Sz);
R2S_Estr=1./T2S_ms_Estr;
R2S_Estr=max(1e-5,min(0.1,R2S_Estr));
disp('OK est resize');