function Out=OperatorTest2(A,ImSize,DataSize)
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = A*x;
Aty = A'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
if(nargin<5)
    disp(['Operator conj test: ' num2str(Out)]);
else
    Err=Out;
    Out=Out<1e-10;
    if(~Out)
%         error('Failed OperatorTest');
        disp(['Operator conj test: ' num2str(Err)]);
    end
end