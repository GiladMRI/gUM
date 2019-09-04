% function Out=FitToModel_1CSf(In,H)
In=squeeze(UpdatedX0);
nCurEchos=size(In,3);
TH=H*In;
[ ~, ~, V_tmp] = batch_svd(TH);
% [ ~, ~, V_tmp] = batch_svd(H_AllTS4_3*In);
% [ U_tmp, s_tmp, V_tmp] = batch_svd(H_AllTS4_3*In);
% s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);
VH12=perm43(conj(V_tmp(:,:,:,1:2))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
MainVHFac=VH12(:,:,:,1);
VH12N=VH12./MainVHFac;
a=VH12N(:,:,1,2);
b=VH12N(:,:,1,3);
c=VH12N(:,:,2,2);
d=VH12N(:,:,2,3);
Denom=(2*(a - c));
BasePart=b-d;
Discriminant=4*a.^2.*d - 4*a.*b.*c - 4*a.*c.*d + b.^2 + 4*b.*c.^2 - 2*b.*d + d.^2;
Sol1=(BasePart+sqrt(Discriminant))./Denom;
Sol2=(BasePart-sqrt(Discriminant))./Denom;
Sol=cat(3,Sol1,Sol2);
Sol=min(1,abs(Sol)).*exp(1i*angle(Sol));

R1TSC=perm43(Sol).^(perm32(0:(nCurEchos-1)));

ad=sum(abs(R1TSC).^2,3);
bc=sum(R1TSC(:,:,:,1).*conj(R1TSC(:,:,:,2)),3);
Det=ad(:,:,1,1).*ad(:,:,1,2) - bc.*conj(bc); % a*d - b*c
IMat=cat(3,cat(4,ad(:,:,1,2),-(bc)),cat(4,-conj(bc),ad(:,:,1,1)))./Det;
% IMat=cat(3,cat(4,ad(:,:,1,2),-conj(bc)),cat(4,-(bc),ad(:,:,1,1)))./Det;
IMat5=permute(IMat,[1 2 5 3 4]);
R1TSCIMat=squeeze(sum(conj(R1TSC).*IMat5,4));
MFac=sum(In.*R1TSCIMat,3);
Fitted=sum(MFac.*R1TSC,4);
Fitted(isnan(Fitted))=0;
%%
syms a b c d w
eqn=(a.*w+c.*(1-w))^2==b.*w+d.*(1-w);
sol=solve(eqn,w);
qs=a*sol+c*(1-sol);
sqs=simplify(qs)
%   (b - d + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
%  -(d - b + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
(b - d + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
Denom=(2*(a - c))
BasePart=b-d
Discriminant=4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2

Sol1=(BasePart+sqrt(Discriminant))/Denom;
Sol2=(BasePart-sqrt(Discriminant))/Denom;

Sol=cat(3,Sol1,Sol2);
Sol=min(1,abs(Sol)).*exp(1i*angle(Sol));

R1TSC=perm43(Sol).^(perm32(0:(nCurEchos-1)));

ad=sum(abs(R1TSC).^2,3);
bc=sum(R1TSC(:,:,:,1).*conj(R1TSC(:,:,:,2)),3);
Det=ad(:,:,1,1).*ad(:,:,1,2) - bc.*conj(bc); % a*d - b*c
IMat=cat(3,cat(4,ad(:,:,1,2),-conj(bc)),cat(4,-(bc),ad(:,:,1,1)))./Det;
IMat5=permute(IMat,[1 2 5 3 4]);
R1TSCIMat=squeeze(sum(conj(R1TSC).*IMat5,4));
% MFac=
%%
% U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
VH1=perm43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
% Normalize VH1
MainVH1Fac=VH1(:,:,1,1);
% U1=U1.*MainVH1Fac;
VH1N=VH1./MainVH1Fac;
% now VH1N 2 is the R1
VH1N(:,:,:,2)=min(abs(VH1N(:,:,:,1)),abs(VH1N(:,:,:,2))).*exp(1i*angle(VH1N(:,:,:,2)));
R1_tmp=VH1N(:,:,:,2);
% TSC vec from 
R1TSC=R1_tmp.^(perm32(0:(nCurEchos-1)));
MFac=sum(In.*conj(R1TSC),3)./gsss(R1TSC,3);
Out=R1TSC.*MFac;
%%
syms a b c d
M=[a b; c d]
% [ a, b]
% [ c, d]
inv(M)
% [  d/(a*d - b*c), -b/(a*d - b*c)]
% [ -c/(a*d - b*c),  a/(a*d - b*c)]
%%
a=[1 3 5+3i]
b=[1 4 9; 4 5 6]
x=a/b;
b(1,:)*x(1)+b(2,:)*x(2)

ad=sum(abs(b).^2,2)
bc=sum(b(1,:).*conj(b(2,:)))
Det=ad(1)*ad(2)-bc*conj(bc)
IMat=[ad(2) -bc; -conj(bc) ad(1)]./Det
bIMat=b'*IMat

[a*bIMat ; x]