function St=getParamsStructFromFN(FN)
if(FN(end)==filesep)
    FN=[FN 'ParamsUsed.txt'];
end
Lines=getLines(FN);
LinesC=regexp(Lines,' ','split');
St=struct();
for i=1:numel(LinesC)
    St.(LinesC{i}{1})=LinesC{i}{2};
end