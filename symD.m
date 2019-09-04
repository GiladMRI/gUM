function Out=symD(In,Dim)
N=size(In,Dim);
D=diff(In,1,Dim);
D1=SliceFromDim(D,1,Dim);
Dend=SliceFromDim(D,N-1,Dim);

DA=cat(Dim,D1,D);
DB=cat(Dim,D,Dend);

Out=(DA+DB)/2;