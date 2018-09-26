function Out=imfillholesBySlices(In)
Out=zeros(size(In));
for i=1:prod(gsize(In,3:20))
    Out(:,:,i)=imfill(In(:,:,i),'holes');
end