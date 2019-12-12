function Out=readFSLMatFile(FN)
QQ=getLines(FN);
for i=1:4
    Out(i,:)=str2num(QQ{i});
end