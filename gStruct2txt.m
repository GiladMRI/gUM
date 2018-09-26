function Txt=gStruct2txt(S,FN)
fid=fopen(FN,'w+');
flds=fieldnames(S);
for i=1:numel(flds)
    CurLine='';
    CurFld=flds{i};
    CurVal=S.(CurFld);
    CurStr=CurVal;
    if(isnumeric(CurVal))
        CurStr=num2str(CurVal);
    end
    CurLine=sprintf('%s',[CurFld ' ' CurStr]);
    Txt{i}=CurLine;
    fprintf(fid,'%s',CurLine);
    if(i<numel(flds))
        fprintf(fid,'\r\n');
    end
end
fclose(fid);
[status,msg,msgID] = fileattrib(FN,'+w','a');
Txt=Txt';