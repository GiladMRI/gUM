function S=gLines2Struct(L)
S=struct();
for i=1:numel(L)
    CurLine=L{i};
    X=regexp(CurLine,' ','split');
    Fld=X{1};
    Val=X{2};
    ValN=str2num(Val);
    if(~isempty(ValN))
        Val=ValN;
    end
    S.(Fld)=Val;
end