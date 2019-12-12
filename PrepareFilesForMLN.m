
SensCCA=SensCC;
%%
BaseBaseOutP='/autofs/cluster/kawin/Gilad/TF/MEx_Ma5x5/';
mkdir(BaseBaseOutP);
system(['chmod -R 777 ' BaseBaseOutP]);
%%
Sensr=SensCCA;
% sig=QQ.CurSig;
sig=DataCCP;
% DataCCP(:,TrajPartMed,CurReps,1:nccToUse)
% B0_Hz=-UpdatedB0Map1;
B0_Hz=UpdatedB0Map1;
%%
WhichRep=1;
TrajPartToUseMLN=1:12500;
SliI

AmpFac=300;
MLNParamsPrefix='/autofs/cluster/kawin/Gilad/TF/ParamsForME5x5S';

nCh=13;
batchSize=16;
nTS=15;
nTSOut=8;
%%
BaseOutDir=[BaseBaseOutP 'Sli' num2str(SliI) filesep];
mkdir(BaseOutDir);
system(['chmod -R 777 ' BaseOutDir]);
SensCC=squeeze(Sensr(:,:,1:nCh)); % [X Y Ch]
SensMsk=grmss(SensCC,3)>0.01; % [X Y]
save([BaseOutDir 'SensCC1.mat'],'SensCC','SensMsk');

% AcqDwellTime_us=2*1.1;
% AcqDwellTime_us=3*1.1;
% TrajPartToUse=1:24000;
% TrajPartToUse=1:45501;
% CurSig=     sig(1,TrajPartToUse(1:2:end),WhichRep,1:nCh)+...
%             sig(1,TrajPartToUse(2:2:end),WhichRep,1:nCh);
CurSig=     sig(1,TrajPartToUseMLN,WhichRep,1:nCh);
CurSig=squeeze(CurSig);
tmp=Row(CurSig);
tmp2=[real(tmp) imag(tmp)]*AmpFac;
Data=tmp2;
Data(batchSize,end)=0;
save([BaseOutDir 'RealDataForNN.mat'],'Data');

nTrajA=size(CurSig,1);
TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TimePoints_ms3=permute(TimePoints_ms,[1 3 2]);
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TS_TimePoints3=permute(TS_TimePoints,[1 3 2]);
TSBF=GetTSCoeffsByLinear(nTrajA,nTS).';
% WhichRSToUse=1;
% TSC=exp(-TS_TimePoints3./UpdatedT2SMap_ms_RS(:,:,WhichRSToUse)).*exp(-1i*2*pi*UpdatedB0Map_RS(:,:,WhichRSToUse).*TS_TimePoints3/1e3);
% TSBF: [15×7162 double]
% TSC: [128×128×15 double]
% B0_Hz=UpdatedB0Map_RS(:,:,WhichRSToUse);
% B0_Hz=RefB0MrS(:,:,SliI);
save([BaseOutDir 'B0TS.mat'],'TSBF','B0_Hz');
%
% Traj=(TrajM(WhichRep,TrajPartToUse(1:2:end))+TrajM(WhichRep,TrajPartToUse(2:2:end)))/2;
% Traj=(TrajM(WhichRep,TrajPartToUse(1:3:end))+TrajM(WhichRep,TrajPartToUse(2:3:end))+TrajM(WhichRep,TrajPartToUse(3:3:end)))/3;
Traj=TrajM(WhichRep,TrajPartToUseMLN(1:1:end));
Sz128=gsize(SensCC,1:2);
clear Trajm2
Trajm2(1,:)=real(Traj);
Trajm2(2,:)=imag(Traj);
[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
P=st.nufftStruct.p/sqrt(prod(Sz128));
save([BaseOutDir 'TrajForNUFT.mat'],'Trajm2','SN','Kd','P');

TimePoints_ms=(1:nTrajA)*AcqDwellTime_us/1000;
TS_TimePoints=linspace(0,TimePoints_ms(end),nTS);
TSstr=strrep(num2str(TS_TimePoints,'%3.5f,'),' ','');
TSstr=TSstr(1:end-1);
% TSstr=['TimePoints_ms ' TSstr(1:end-1)];

St=getParamsStructFromFN('/autofs/cluster/kawin/Gilad/TF/ParamsForS10.txt');
St.TimePoints_ms=TSstr;
TS_TimePointsForRec=linspace(0,TimePoints_ms(end),nTSOut);
TS_TimePointsForRec_Str=strrep(num2str(TS_TimePointsForRec,'%3.5f,'),' ','');
TS_TimePointsForRec_Str=TS_TimePointsForRec_Str(1:end-1);
St.TimePointsForRec_ms=TS_TimePointsForRec_Str;
St.SessionNameBase=['ME_S' num2str(SliI)];
St.RealDataFN=[BaseOutDir 'RealDataForNN.mat'];
St.BaseTSDataP=BaseOutDir;
St.BaseNUFTDataP=BaseOutDir;
% disp(TSstr)
ParamsOutFn=[MLNParamsPrefix num2str(SliI) '.txt'];
Txt=gStruct2txt(St,ParamsOutFn);
disp(['wrote ' ParamsOutFn]);
SensCC=SensCCA;