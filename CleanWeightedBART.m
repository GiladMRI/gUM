BaseSP='/autofs/space/daisy_002/users/Gilad/';

FmacWOp={'fmac 0 0'};
FmacWScriptFN=[BaseSP 'FmacW.txt'];
WriteLinopToFile(FmacWScriptFN,FmacWOp);

tmpImSz16=FillOnesTo16(size(m0));
%%
fgmontage(b0,[-200 200]);
%%
% Clean=bart(['picsS -m -u 0.1 -w 1 -R W:3:0:100 ' FmacWScriptFN],tmpImSz16,b0.*m0,m0);
WW=m0.^2;
% Clean=bart(['picsS -m -u 10 -w 1 -R W:3:0:1000 ' FmacWScriptFN],tmpImSz16,b0.*WW,WW);
Clean=bart(['picsS -m -u 1e2 -w 1 -R D:3:0:10000 ' FmacWScriptFN],tmpImSz16,b0.*WW,WW);
fgmontage(Clean,[-200 200]);