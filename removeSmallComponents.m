function Out=removeSmallComponents(Msk)
[L,n]=bwlabel(Msk);
H=histcounts(L(:),0:(n+1));
if(numel(H)==1)
    Out=Msk;
    return;
end
[~,MI]=max(H(2:end));
Out=(L==MI);