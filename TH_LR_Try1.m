% Temporal Hankelization Low Rank
% "Image" size is [X Y Echos]
% Signal size is the same
% Image to signal operator is FT on X,Y and masking

ImSz=[Sz nEchos];
ImSzF=FillOnesTo16(ImSz);

Main_Ops={'fftc 3','fmac 0 0'};
% TH_Ops='hankel 3 4 2';

TH_LR_FN='TH_LR.txt';
WriteLinopToFile(TH_LR_FN,Main_Ops);
% WriteLinopToFile(TH_LR_FN,{Main_Ops,TH_Ops});

THLRLambda=0.003;
% TH_LR_Rec=bart(['picsS -m -R K:4:3:' num2str(THLRLambda) ':2:1:0:2 ' TH_LR_FN],ImSzF,Sig,Msk); % 11,2
TH_LR_Rec=bart(['picsS -m -R K:8:3:' num2str(THLRLambda) ':2:1:0:2 ' TH_LR_FN],ImSzF,Sig,Msk); % 2,11
Fac1=168.5870;
WarmSartFN='WarmStart';
writecfl(WarmSartFN,X/Fac1);
% TH_LR_Rec=bart(['picsS -m  -W ' WarmSartFN ' ' TH_LR_FN],ImSzF,Sig,Msk); % 2,11
TH_LR_Rec=bart(['picsS -m -R K:8:3:' num2str(THLRLambda) ':2:1:0:2 -W ' WarmSartFN ' ' TH_LR_FN],ImSzF,Sig,Msk); % 2,11

WLambda=0.003;
WTH_LR_Rec=bart(['picsS -m -R K:8:3:' num2str(THLRLambda) ':2:1:0:2 -R W:3:0:' num2str(WLambda) ' ' TH_LR_FN],ImSzF,Sig,Msk); % 2,11

THLLRLambda=0.003;
THLLR_LocalSize=5;
TH_LLR_Rec=bart(['picsS -m -R K:8:3:' num2str(THLLRLambda) ':3:' num2str(THLLR_LocalSize) ':0:2 ' TH_LR_FN],ImSzF,Sig,Msk); % 2 704

TVLambda=0.1;
TV_Rec=bart(['picsS -m -R T:3:3:' num2str(TVLambda) ' ' TH_LR_FN],ImSzF,Sig,Msk);

NTH_LR_Rec=TH_LR_Rec.*grmss(X)./grmss(TH_LR_Rec);
NWTH_LR_Rec=WTH_LR_Rec.*grmss(X)./grmss(WTH_LR_Rec);
NTH_LLR_Rec=TH_LLR_Rec.*grmss(X)./grmss(TH_LLR_Rec);
NTV_Rec=TV_Rec.*grmss(X)./grmss(TV_Rec);

THLRs=ssimBySlices(NTH_LR_Rec,X);
WTHLRs=ssimBySlices(NWTH_LR_Rec,X);
THLLRs=ssimBySlices(NTH_LLR_Rec,X);
TVs=ssimBySlices(NTV_Rec,X);

figure;plot(TVs,'k');hold on;plot(THLRs,'r');plot(WTHLRs,'r--');plot(THLLRs,'g');