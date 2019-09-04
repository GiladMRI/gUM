function Out=FitToModel_1CSfx(In,H,N)
% In=squeeze(UpdatedX0);
nCurEchos=size(In,3);
[ ~, ~, V_tmp] = batch_svd(H*In);
% [ U_tmp, s_tmp, V_tmp] = batch_svd(H_AllTS4_2*In);
% s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);

% U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
VH1=perm43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));
% Normalize VH1
MainVH1Fac=VH1(:,:,1,1);
% U1=U1.*MainVH1Fac;
VH1N=VH1./MainVH1Fac;
% now VH1N 2 is the R1
% VH1N(:,:,:,2)=min(abs(VH1N(:,:,:,1)),abs(VH1N(:,:,:,2))).*exp(1i*angle(VH1N(:,:,:,2)));
VH1N(:,:,:,2)=min(abs(VH1N),[],4).*exp(1i*angle(VH1N(:,:,:,2)));
R1_tmp=VH1N(:,:,:,2);
% TSC vec from 
R1TSC=R1_tmp.^(perm32(0:(nCurEchos-1)));
MFac=sum(In.*conj(R1TSC),3)./gsss(R1TSC,3);

BMFac=abs(MFac)>0 & isfinite(MFac);
Medval=median(abs(MFac(BMFac)));
MFac=min(abs(MFac),Medval*6).*exp(1i*angle(MFac));

R1TSCx=R1_tmp.^(perm32( (0:(N-1))+nCurEchos-N ));
Out=R1TSCx.*MFac;
Out(isnan(Out))=0;