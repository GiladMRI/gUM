function Out=getMostRecentDirWithBase(TBaseP,CurTrainDirbase)

D=dir([TBaseP CurTrainDirbase '*']);
D=D([D.isdir]);
DName={D.name}';
D=D(strhas(DName,'_train'));
DNameFull={D.name}';
DName=cell(0);
for i=1:numel(DNameFull)
    DName{i}=DNameFull{i}((numel(CurTrainDirbase)+1):end-6);
end
DT=datenum(DName,'yyyy-mm-dd_HH-MM-SS');
[~,MI]=max(DT);
Out=DNameFull{MI};