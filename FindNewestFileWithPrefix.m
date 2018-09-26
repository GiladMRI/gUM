function [FN, FFN]=FindNewestFileWithPrefix(P,Pre,Ext)
if(P(end)~=filesep)
    P=[P filesep];
end
if(nargin<3)
    Ext='*';
end
D=dir([P Pre '*.' Ext]);
if(isempty(D))
    FN='';
    FFN='';
    return;
end
[~,MI]=max([D.datenum]);
FN=D(MI).name;
FFN=[P FN];