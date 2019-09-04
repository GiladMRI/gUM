function [M, dT, TimePointsR] =GetTSCoeffsByLinearWithPlateau(N,n)
% N=5000;
% n=5;
if(n==1)
    M=ones(N,1);
    return;
end
M=zeros(N,n);
% Ttimes=linspace(0,1,2*n+1);
Ttimes=linspace(-0.5,N+0.5,2*n+1)/N;
Ttimes(3:2:end-1)=[];
OutTimes=linspace(0,1,N);

dT=Ttimes(3)-Ttimes(2);
TimePointsR=Ttimes(2:end-1);

tmp=zeros(1,n+2);
tmp(1:2)=1;
M(:,1)=interp1(Ttimes,tmp,OutTimes);

tmp=zeros(1,n+2);
tmp(end-1:end)=1;
M(:,end)=interp1(Ttimes,tmp,OutTimes);

for i=2:n-1
    tmp=zeros(1,n+2);
    tmp(i+1)=1;
    M(:,i)=interp1(Ttimes,tmp,OutTimes);
end
%
% figure;plot(M,'-*');hold on;
% for i=1:numel(Ttimes)
%     plot([Ttimes(i) Ttimes(i)]*(N-1)+1,[0 1],'b');
% end