% [x,time]=socp2sdp(f,A,b,N,x,z,Nu,abs_tol,rel_tol,max_iter)

function [x,time]=socp2sdp(f,A,b,N,x,z,Nu,abs_tol,rel_tol,max_iter)

L=length(N);
n=length(f);

k=0;
Z=[];
F=[];
for i=1:L,
	Fi0=[b(k+N(i))*eye(N(i)-1) b(k+1:k+N(i)-1);
		b(k+1:k+N(i)-1)' b(k+N(i))];
	Fi=Fi0(:);
	for j=1:n,
		Fij=[A(k+N(i),j)*eye(N(i)-1) A(k+1:k+N(i)-1,j);
			A(k+1:k+N(i)-1,j)' A(k+N(i),j)];
		Fi=[Fi Fij(:)];
	end;
	F=[F; Fi];

	Zi=[z(k+N(i))/N(i)*eye(N(i)-1) z(k+1:k+N(i)-1)/2;
		z(k+1:k+N(i)-1)'/2 z(k+N(i))/N(i)];
	Z=[Z; Zi(:)];
	k=k+N(i);
end;

%keyboard;

[x,Z,ul,info,time] = ...
             sp(F,N,f,x,Z,Nu,abs_tol,rel_tol,0,max_iter);
