function In=gifftshift(In,Dims)
if(nargin<2)
    Dims=1:ndims(In);
end
for i=1:numel(Dims)
    In=ifftshift(In,Dims(i));
end