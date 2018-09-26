ParamSet2=UsedParams;
ParamSet2.BankSize=0;

ParamSet3=UsedParams;
ParamSet3.nTimeSegments=5;

ParamSet4=UsedParams;
ParamSet4.learning_rate_start=0.003;

ParameterSets={ParamSet2 ParamSet3 ParamSet4};

ResOnAllSets=zeros(1,ParameterSets);
for SetIdx=1:numel(ParameterSets)
    try
        CurSet=ParameterSets{SetIdx};
        Txt=gStruct2txt(CurSet,[BaseTFFolder 'Params.txt']);
        system('sudo /home/deni/RunTFForMatlabx.sh');
        ResOnAllSets(SetIdx)=1;
    catch
        ResOnAllSets(SetIdx)=-1;
    end
end