function [FN, FFN]=FindNewestFolderWithPrefix(P,Pre,Str)
D=dir([P Pre '*.*']);
D2=dir([P Pre '*.']);
D3=dir([P Pre '*']);
D=[D; D2; D3];
D=D([D.isdir]);
if(nargin>2)
    DN={D.name};
    D=D(strhas(DN,Str));
end
if(isempty(D))
    FN='';
    FFN='';
    return;
end
[~,MI]=max([D.datenum]);
FN=D(MI).name;
FFN=[P FN];