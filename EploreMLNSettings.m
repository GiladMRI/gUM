CurOurPSli=[ScanP BaseFN filesep 'Sli' num2str(SliI,'%02d') filesep];

%     ParamsSDefaults=getParamsStructFromFN('/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry3C2_7TS_RL_S3__2018-07-16_15-19-07_train/');
    ParamsSDefaults=getParamsStructFromFN('/home/deni/TF/srez/RegridTry3C2_7TS_S02_Sli01__2018-07-23_09-40-59_train/');
    ParamsS=ParamsSDefaults;
    ParamsS.SessionNameBase=['RegridTry3C2_7TS_' ScanP(end-3:end-1) '_Sli' num2str(SliI,'%02d')];
    ParamsS.RealDataFN=[CurOurPSli 'RealDataForNN.mat'];
    ParamsS.BaseTSDataP=CurOurPSli;
    ParamsS.BaseNUFTDataP=[ScanP BaseFN filesep];
%     Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    Txt=gStruct2txt(ParamsS,'~/TF/Params.txt');
    
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
    system('sudo /home/deni/RunTFForMatlabx.sh');
    
    %%
    ParamsS.nNeighbors=15;
    ParamsS.BaseNUFTDataP=[ScanP BaseFN filesep];
%     Txt=gStruct2txt(ParamsS,'~/HomeA/TF/Params.txt');
    Txt=gStruct2txt(ParamsS,'~/TF/Params.txt');
    
%     system('sudo -H -u a /media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/RunTFForMatlabx.sh');
    system('sudo /home/deni/RunTFForMatlabx.sh');
    