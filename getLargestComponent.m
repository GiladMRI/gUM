function Out=getLargestComponent(Msk)
Out=false(size(Msk));
for i=1:prod(gsize(Msk,3:ndims(Msk)))
    [L,n]=bwlabel(Msk(:,:,i),4);
    H=histcounts(L(:),0:(n+1));
    if(numel(H)==1)
        Out(:,:,i)=Msk(:,:,i);
        continue;
    end
    [~,MI]=max(H(2:end));
    Out(:,:,i)=(L==MI);
end