I1=phantom(128);
II=bart('linopScript /autofs/space/daisy_002/users/Gilad/gUM/dblszS.txt',FillOnesTo16(size(I1)),I1);

% STraj=linspace(0,50,10000).*exp(1i*linspace(0,2*pi*50,10000));
STraj=kC.';
STraj3=CTo3Rows(STraj);

% SnufftStruct = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),size(I1)), size(I1), [6 6], size(I1)*2, size(I1)/2); % st.om
SnufftStruct = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(STraj),size(I1)), size(I1), [6 6], size(I1)*2); % st.om

% Kern=NUFFT_to_Toep_2blocks(SnufftStruct,ones(numel(STraj),1));
Kern=NUFFT_to_Toep_2blocks(SnufftStruct,TSB(:,4));
NUFTB=bart('linopScript /autofs/space/daisy_002/users/Gilad/gUM/nuftScript.txt',FillOnesTo16(size(I1)),I1,STraj3,Kern);
NUFTB=NUFTB.*(TSB(:,4).');
NUFTB_NA=bart('linopScript -A /autofs/space/daisy_002/users/Gilad/gUM/nuftScript.txt',FillOnesTo16(size(I1)),NUFTB,STraj3,Kern);
ShowAbsAngle((NUFTB_NA))

NUFTB_N=bart('linopScript -N /autofs/space/daisy_002/users/Gilad/gUM/nuftScript.txt',FillOnesTo16(size(I1)),I1,STraj3);
ShowAbsAngle((NUFTB_N))

NUFTB_N2=bart('linopScript -N /autofs/space/daisy_002/users/Gilad/gUM/nuftScriptN.txt',FillOnesTo16(size(I1)),I1,STraj3,(Kern));
ShowAbsAngle((NUFTB_N2))

NUF=nufft(I1,SnufftStruct);
NNUF=nufft_adj(NUF,SnufftStruct);
F=fft2(I1,size(I1,1)*2,size(I1,2)*2);
FK=F.*Kern;
IFK=ifft2(FK);
CIFK=IFK(1:size(I1,1),1:size(I1,2));
%%


ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_Script.txt';
RecrDecayMagLLRC=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP);
RecrDecayMagLLRCT=sum(RecrDecayMagLLRC.*TSB_EchosP,7);
%% with Teoplitz
ScriptFN_Toep='/autofs/space/daisy_002/users/Gilad/gUM/NuftTS_ScriptN.txt';

TSKern=NUFFT_to_Toep_2blocks(nufftStructR, TSBk);
TSKernP=permute(TSKern,[1 2 7 6 5 4 3]);

RecrDecayMagLLRC_Toep=bart(['picsS -m -u 0.010 -b 3 -R L:3:3:0.00001 ' ScriptFN_Toep],[Sz 1 1 1 1 nTS ones(1,9)],DataWithDecay_MagOnlyP,Sensr,CTo3Rows(kC(1:end,1).'),TSBkP,TSKernP);
RecrDecayMagLLRCT_Toep=sum(RecrDecayMagLLRC_Toep.*TSB_EchosP,7);
