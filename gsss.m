function Out=gsss(Data,Dims)
if(nargin<2)
    Dims=1:ndims(Data);
end
Out=squeeze(gss(Data,Dims));