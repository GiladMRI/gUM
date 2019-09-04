load('maps128x128x8.mat');
nChannles=size(maps,3);

I=phantom(128);
IWithSens=I.*maps;

F=fft2cg(I);
M=ones(size(I));

FS=fft2cg(IWithSens);
%%
MM=conj(maps).*permute(maps,[1 2 6 5 4 3]);
MMI=MM-permute(eye(nChannles),[3 4 1 5 6 2]);

writecfl('/autofs/space/daisy_002/users/Gilad/gUM/MMI',MMI);

grmss(sum(IWithSens.*MMI,3))
%%
ET36=permute(eye(nChannles),[3 4 1 5 6 2]);
writecfl('/autofs/space/daisy_002/users/Gilad/gUM/ET36',ET36);
%%
BaseSP='/autofs/space/daisy_002/users/Gilad/';
ScriptFN=[BaseSP 'NoC.txt'];
%% Single channel fully sampled
SzI=size(I);
ImSize16=FillOnesTo16(SzI);

RecM=bart(['picsS -R T:3:3:0.01 -m ' ScriptFN],ImSize16,F,M);
%%
SzC=size(IWithSens);
ImSizeC16=FillOnesTo16(SzC);
RecC=bart(['picsS -R T:3:3:0.01 -m ' ScriptFN],ImSizeC16,FS,M);
%%
RecC=bart(['picsS -R C:0.01 -m ' ScriptFN],ImSizeC16,FS,M);
%%
QQ=bart(['linopScript ' ScriptFN],ImSize16,m,T,PD,SensCSMap,mask_sampleP);

RecM=bart(['picsS -R W:3:0:0.05 -m ' ScriptFN],ImSize16,y,T,PD,SensCSMap,mask_sampleP);
