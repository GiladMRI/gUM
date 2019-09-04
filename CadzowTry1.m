% M=(randn(5,2)+1i*randn(5,2))*10;
M=(rand(5,2)+1i*rand(5,2))*10;
Mb=M+randn(5,2)*4;

[U0,S0,V0]=svd(M,0);
grmss(M-U0*S0*V0')
D=grmss(M-U0(:,1)*S0(1,1)*V0(:,1)')
U=U0(:,1)*sqrt(S0(1,1));
V=V0(:,1)*sqrt(S0(1,1));
U2=U0(:,2)*sqrt(S0(2,2));
V2=V0(:,2)*sqrt(S0(2,2));
grmss(M-U*V'-U2*V2')

[U0b,S0b,V0b]=svd(Mb,0);
Db=grmss(Mb-U0b(:,1)*S0b(1,1)*V0b(:,1)')
% Db=grmss(Mb-U0b(:,2)*S0b(2,2)*V0b(:,2)')
% Cadzow
K=20;
for i=1:K
%     U=Mb/V';
%     V=(U\Mb)';
    U=Mb*V/gsss(V);
    V=(U'*Mb)'/gsss(U);
end

Dbx=grmss(Mb-U*V')
%%
% v1: [1 a ab]
% v2: [1 c cd]
% have a balanced linear combination, w, 1-w:
% aw+c(1-w), abw+cd(1-w)
% (aw+c(1-w))^2=abw+cd(1-w)
syms a b c d w
eqn=(a*w+c*(1-w))^2==a*b*w+c*d*(1-w);
sol=solve(eqn,w);
qs=a*sol+c*(1-sol);
simplify(qs)
% w1=  ((a^2*b^2 - 4*a^2*b*c + 4*a^2*c*d + 4*a*b*c^2 - 2*a*b*c*d - 4*a*c^2*d + c^2*d^2)^(1/2) + a*b - 2*a*c - c*d + 2*c^2)/(2*(a - c)^2)
% w2= -((a^2*b^2 - 4*a^2*b*c + 4*a^2*c*d + 4*a*b*c^2 - 2*a*b*c*d - 4*a*c^2*d + c^2*d^2)^(1/2) - a*b + 2*a*c + c*d - 2*c^2)/(2*(a - c)^2)
% q1=a*w1+c*(1-w1);
% q2=a*w2+c*(1-w2);
% simplify(q1)
% simplify(q2)