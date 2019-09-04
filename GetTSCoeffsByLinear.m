function M=GetTSCoeffsByLinear(N,n)
% N=5000;
% n=5;
if(n==1)
    M=ones(N,1);
    return;
end
M=zeros(N,n);
Ttimes=linspace(0,1,n);
for i=1:n
    tmp=zeros(1,n);
    tmp(i)=1;
    M(:,i)=interp1(Ttimes,tmp,linspace(0,1,N));
end