if(~exist('Ma7TMPBD','var'))
    Ma7TMPBD=load('/autofs/cluster/kawin/Gilad/All_Orientation-0x.mat');
end
%%
CurSli=permute(Ma7TMPBD.CurSetAll(:,:,:,134),[2 3 1 4]);

C=CurSli(:,:,1).*exp(1i*CurSli(:,:,2));
B0_Hz=CurSli(:,:,4);
B0_Hz(~isfinite(B0_Hz))=0;
% D=CurSli(:,:,3);
T2S_ms=CurSli(:,:,3);
T2S_ms=max(4,min(T2S_ms,300));

nEchos=30;
TimeBetweenEchos_ms=2;
EchoTimes_ms=TimeBetweenEchos_ms*(0:nEchos-1);
EchoTimes_ms3=permute32(EchoTimes_ms);
R2S=1./T2S_ms;

TSC=exp(-EchoTimes_ms3./T2S_ms).*exp(-1i*2*pi*EchoTimes_ms3.*B0_Hz/1000);

X=C.*TSC;
X(isnan(X))=0;
disp('OK');
%%
SmRange=[15 15];
SmSig=2;
Sm_X=SmoothBySlices(X,SmRange,SmSig);

SB0Var_Hz=GenerateRandomB0Variation(gsize(X,1:2),2,.1,3);

S_TSC=exp(-1i*2*pi*EchoTimes_ms3.*SB0Var_Hz/1000);

X_Est=Sm_X.*S_TSC;

B0_Hz_Est=SmoothBySlices(B0_Hz+SB0Var_Hz,SmRange,SmSig);

T2S_ms_Est=SmoothBySlices(T2S_ms,SmRange,SmSig);
disp('ok est');