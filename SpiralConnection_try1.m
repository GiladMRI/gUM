N=42;
N2=3;
LocX0=1500;
LocXm1=1450;
LocYm1=-80;
MaxAcc=70;
MaxV=200;
%
Vy0=-LocYm1;
ay=2*Vy0/(N+1);
Vy=Vy0-(1:N)*ay;
Y=cumsum(Vy);

Vx0=LocX0-LocXm1;
ax=(Vx0*N+LocX0)/(N2*(N2+1)/2+N2*(N-N2))
Vx=Vx0-[(1:N2) ones(1,N-N2)*N2]*ax;
X=LocX0+cumsum(Vx);

% X=linspace(LocX0,0,N);

figure;
plot(LocXm1,LocYm1,'k+');hold on
plot(LocX0,0,'k+');
plot(X,Y,'r.');