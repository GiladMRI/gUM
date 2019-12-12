function BARTS_Out=WriteBARTStructToFiles(BARTStruct,BaseFP)
Flds=setdiff(fieldnames(BARTStruct),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTStruct.(Flds{i})))
        for j=1:numel(BARTStruct.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTStruct.(Flds{i}){j});
            BARTStruct.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTStruct.(Flds{i}));
        BARTStruct.(Flds{i})=tmpFN;
    end
end
BARTS_Out=BARTStruct;