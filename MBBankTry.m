%% Call TF
% CurOurPSli='/media/a/DATA/MBBankTry/Sli05/';
% ParamsSDefaults=getParamsStructFromFN('/media/a/DATA/ASLSubjDataCopy/S04/meas_MID149_gBP_VD11_U19_G35S155_FID23846/RegridTry3C2_7TS_S04_Sli01__2018-07-19_01-47-17_train/');
ParamsSDefaults=getParamsStructFromFN('/home/deni/TF/srez/RegridTry3C2_7TS_S02_Sli01__2018-07-23_09-40-59_train/');
ParamsS=ParamsSDefaults;
% MIDStr='XXX';
ParamsS.SessionNameBase=['RegridTry3C2_7TS_' MIDStr '_Sli' num2str(SliI,'%02d')];

ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];

ParamsS.BaseTSDataP=CurOurPSli;
ParamsS.BaseNUFTDataP=[ScanP BaseFN filesep];
% ParamsS.BaseNUFTDataP='/media/a/DATA/MBBankTry/Base/';

% ParamsS.DatasetMatFN='/media/a/H1/HCPData_256x256_int16.mat';
ParamsS.DatasetMatFN='/home/deni/HCPData_256x256_int16.mat';
ParamsS.nToLoad=10000;

% Parames for MB:
% nccX=16;
ParamsS.nccInData=ncc;
ParamsS.InputMode='RegridTry3FMB';
ParamsS.NetMode='RegridTry3C2_TS_MB';

% To change and decide
ParamsS.nccToUse=16;    %13
ParamsS.nNeighbors=12;  %12
% ParamsS.nTimeSegments=15; %11
ParamsS.nTimeSegments=11;

% ParamsS.nTimeSegments=5; %11

ParamsS.learning_rate_start=0.002;
ParamsS.learning_rate_half_life=30;
ParamsS.train_time=180;
BatchSizeInNN=6;
ParamsS.batch_size=BatchSizeInNN;
ParamsS.WL2_Lambda=0;

ParamsS.RandomPhaseLinearFac=3;
ParamsS.RandomPhaseQuadraticFac=0.05;
ParamsS.RandomPhaseScaleFac=2;

ParamsS.checkpoint_period=20;

ParamsS.InitForRFN='None';
ParamsS.InitForLFN='None';

% ParamsS.InitForRFN='/home/deni/TF/srez/RegridTry3C2_7TS_MID448_Sli06__2018-09-20_09-17-39_train/TrainSummary_014491.mat';

% ParamsS.BankSize=780;
ParamsS.BankSize=0; % Not use a signal bank
ParamsS.BankSize=1440;
ParamsS.BankK=7;

% ParamsS.nTraj=5118;
ParamsS.DataH=size(Data,2);
% ParamsS.achannelsIn=ParamsS.nNeighbors*ParamsS.nccToUse*2;

ParamsS.kMax=ceil(MaxK);

ParamsS.ShowRealData=1;
% ParamsS.nccToUse*ParamsS.nNeighbors

% DataX=Data;
Data=DataX(1:ParamsS.batch_size,:);
ParamsS.RealDataFN=[CurOurPSli 'RealDataForNNx.mat'];
save(ParamsS.RealDataFN,'Data')
clc

% disp(['Each batch 12, 163776, 1, 1 ' num2str(12*163776*4/1e6,'%.2f') 'MB']);
disp(['Constant images: ' num2str(ParamsS.nToLoad*256*256*2/1e9,'%.2f') 'GB']);
disp(['Regridded batch data ' num2str(128*128*ParamsS.nNeighbors*ParamsS.nccToUse*4*2*ParamsS.batch_size/1e6,'%.2f') 'MB']);
disp(['Kside ' num2str(ParamsS.batch_size*128*128*ParamsS.nNeighbors*ParamsS.nccToUse*ParamsS.nTimeSegments*2*4*2/1e9,'%.2f') 'GB']);
disp(['Signal Bank ' num2str(ParamsS.DataH*4*ParamsS.BankSize*2/1e9,'%.2f') 'GB']);
disp(['Image Bank ' num2str(128*128*4*2*ParamsS.BankSize*2/1e6,'%.2f') 'MB']);
disp(num2str(ParamsS.BankSize/ParamsS.batch_size));

Txt=gStruct2txt(ParamsS,'~/TF/Params.txt');

system('sudo /home/deni/RunTFForMatlabx.sh');

% Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
% system('sudo -H -u a /media/a/H2/home/a/RunTFForMatlabx.sh');