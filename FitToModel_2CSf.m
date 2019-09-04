function Out=FitToModel_2CSf(In,H)
% function Out=FitToModel_1CSf(In,H)
% In=squeeze(UpdatedX0);
nCurEchos=size(In,3);
% TH=H*In;
% [ ~, ~, V_tmp] = batch_svd(TH);
[ ~, ~, V_tmp] = batch_svd(H*In);
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
Out=sum(MFac.*R1TSC,4);
Out(isnan(Out))=0;