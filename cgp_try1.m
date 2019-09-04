tic;
S=rand(5);
S=S+S';
x0=zeros(5,1);
b=rand(5,1);
[x, t] = cgp(x0, S, speye(1), b, 3000, 10^-8, @(Z, o) Z*o, @(Z, o) o);
toc
%%
% I1
SN=nufftStruct.sn;
P=nufftStruct.p.G;

Fac=N1*2*N2*2;
forw=@(x) P*reshape(fft2(padarray(x.*SN,[N1,N2],'post')),[],1)*Fac;
back=@(x) subsref(ifft2(reshape(P'*x,[N1*2,N2*2])),struct('type','()','subs',{{1:N1,1:N2}})).*conj(SN);

ToepH=NUFFT_to_Toep(nufftStruct,1);

toep=@(x) subsref(ifft2(fft2(padarray(x,[N1,N2],'post')).*ToepH),struct('type','()','subs',{{1:N1,1:N2}}));
% Padv=v;
% Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N1);


Sig=forw(I1);

AHA_I=back(Sig);
Toep_I=toep(I1);
grmss(AHA_I)/grmss(AHA_I-Toep_I)
%% Just cg
MaxIter=15;
tic;
x0=zeros(prod(size(I1)),1);
[x, ItersUsed] = cgp(x0, 0, 0, reshape(AHA_I,[],1), MaxIter, 10^-8, @(Z, x) reshape(toep(reshape(x,[128,128])),[],1), @(Z, o) o);
toc
ShowAbsAngle(reshape(x,Sz128))
%% Add L2 regularization
Lambda=1e6;

Dx = @(x) diff(x,1,2);
Dy = @(x) diff(x,1,1);
Dxc=@(x) [Dx(x) zeros(128,1)]-[zeros(128,1) Dx(x)];
Dyc=@(x) [Dy(x); zeros(1,128)]-[zeros(1,128); Dy(x)];
% Dxc=@(x) [zeros(128,1) Dx(x)];
% Dyc=@(x) [zeros(1,128); Dy(x)];

toepAndReg=@(x) toep(x)-Lambda*(Dxc(x)+Dyc(x));
tic;
x0=zeros(prod(size(I1)),1);
[x, ItersUsed] = cgp(x0, 0, 0, reshape(AHA_I,[],1), MaxIter, 10^-8, @(Z, x) reshape(toepAndReg(reshape(x,[128,128])),[],1), @(Z, o) o);
toc
ShowAbsAngle(reshape(x,Sz128))
%% Add L2 regularization, 2D
Lambda=1e6;

Dx = @(x) diff(x,1,2);
Dy = @(x) diff(x,1,1);
Dxc=@(x) [Dx(x) zeros(128,1)]-[zeros(128,1) Dx(x)];
Dyc=@(x) [Dy(x); zeros(1,128)]-[zeros(1,128); Dy(x)];
% Dxc=@(x) [zeros(128,1) Dx(x)];
% Dyc=@(x) [zeros(1,128); Dy(x)];

toepAndReg=@(x) toep(x)-Lambda*(Dxc(x)+Dyc(x));
tic;
x0=zeros(size(I1));
[x, ItersUsed] = cgp2D(x0, 0, 0, AHA_I, MaxIter, 10^-8, @(Z, x) toepAndReg(x), @(Z, o) o);
toc
ShowAbsAngle(reshape(x,Sz128))
