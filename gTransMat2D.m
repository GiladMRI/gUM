function Out=gTransMat2D(T,T2)
if(nargin>1)
    T=[T, T2];
end
Out=eye(3);
Out(3,1:2)=T;