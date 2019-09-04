BaseSPENP='/autofs/cluster/kawin/Gilad/SPEN_data/';
SR=load([BaseSPENP 'Amat.mat']);
SR=SR.FinalAFull; % 120 120
SPENData=load([BaseSPENP 'CmplxDataEO.mat']);
SPENData=SPENData.CmplxDataEO; % 120 120 5
Sens=load([BaseSPENP 'SensMaps.mat']);
Sens=Sens.Sens; % 120 120 5
% [SPEN RO Channel Aux]
%% Basic recon
FRO=fft1cg(SPENData,2);
BasicRec=MultMatTensor(SR',FRO);
BasicRecC=sum(BasicRec.*conj(Sens),3);
figure;subplot(1,2,1);imagesc(abs(BasicRecC));subplot(1,2,2);imagesc(angle(BasicRecC));colormap gray;
%% Linop script : Data format is [SPEN_AXIS_Img RO/kRO Ch SPEN_AXIS_Signal]
ImgSz=size(BasicRecC);
DataP=permute(SPENData,[4 2 3 1]); % [1 120 5 120] : [1 RO Ch SPENAUX]
SRP=permute(SR,[2 3 4 1]); % [120 1 1 120] : [SPEN 1 1 SPENAUX]
SPENScriptFN='SPEN.txt';
SPEN_Ops={'fmac 0 0 Apply sensitivities' 'ifftc 2 Fourier along RO','fmac 1 1 Apply SPEN'};
WriteLinopToFile(SPENScriptFN,SPEN_Ops);
% Test
Test1=bart(['linopScript ' SPENScriptFN],FillOnesTo16(size(BasicRecC)),BasicRecC,Sens,SRP);
fgmontage(cat(4,permute(Test1/3,[4 2 3 1]),permute(DataP,[4 2 3 1])));ylabel('Data             Op*BasicRecC','FontSize',20);
%% Now recon
RegularizationStr='-R T:3:3:0.001';
Rec=bart(['picsS -m ' RegularizationStr  ' ' SPENScriptFN],FillOnesTo16(ImgSz),DataP,Sens,SRP);
figure;
subplot(2,2,1);imagesc(abs(BasicRecC));title('Basic rec (SR'')');
subplot(2,2,3);imagesc(angle(BasicRecC));
subplot(2,2,2);imagesc(abs(Rec));title('BART (linopScript) Rec');
subplot(2,2,4);imagesc(angle(Rec));title(RegularizationStr);
colormap gray