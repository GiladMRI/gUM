N1=12;
PadMtx=[eye(N1); zeros(N1)];
F1=dftmtx(N1*2);
F1=fftshift(fftshift(F1,1),2);
% P1=rand(5,N1*2);
P1=zeros(5,N1*2);
% for i=1:5
%     R=rand(1,3);
%     R=R./sum(R);
% %     P1(i,randi(18)+(0:2))=R;
% %     P1(i,randi(18)+(0:2))=[0.25 0.5 0.25];
%     P1(i,i+(0:2))=[0.25 0.5 0.25];
% end
P1(1,4)=1;
% P1(1,4:6)=[0.25 0.5 0.25]/norm([0.25 0.5 0.25]);
P1(2,14)=1;
P1(3,8)=1;
P1(4,7)=1;
P1(5,12)=1;
% PF1=P1*F1*PadMtx;
PF1=P1*F1;
T1=PF1'*PF1;
ShowAbsAngle(T1)
%%
PF1=F1([2,3,5,7,8,9,11:16,19],:);
T1=PF1'*PF1;
T1
%%
N1=16;
nTraj=13;

Sz128=[N1,N1];
Trajm2=(rand(2,nTraj)-0.5)*N1;
% Trajm2(2,:)=0;
[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
% P=st.nufftStruct.p/sqrt(prod(Sz128));
P=st.nufftStruct.p.G/sqrt(prod(Sz128));
%
% N1=48;
DSN=diag(SN(:));
PadRow=kron(kron(eye(N1),[1 0]'),eye(N1));
PadRow(N1*N1*4,1)=0;

% PadMtx=[eye(N1); zeros(N1)];
% FTBase=dftmtx(N1*2)*PadMtx;
FTBase=dftmtx(N1*2);
% FTBase=fftshift(fftshift(FTBase,1),2);
% FTBase=fftshift(FTBase,2);

FT1=kron(eye(N1*2),FTBase);
FT2=kron(FTBase,eye(N1*2));
% F2D=FT1*FT2;
F2D=FT2*FT1*PadRow*DSN;

% FT1=kron(eye(N1),FTBase);
% FT2=kron(FTBase,eye(N1*2));
% % F2D=FT1*FT2;
% F2D=FT2*FT1;

PF1=P*F2D;
T1=PF1'*PF1;
ShowAbsAngle(T1)
%% real
w1=sin(linspace(0,3*pi,nTraj)).'+1i*sin(pi/4+linspace(0,pi,nTraj)).';

PFreal=RepDotMult(w1,PF1);
Treal=PFreal'*PFreal;
ShowAbsAngle(Treal)
%% complex
wc=sin(linspace(0,3*pi,nTraj)).'+1i*sin(pi/4+linspace(0,pi,nTraj)).';

PFc=RepDotMult(wc,PF1);
Tc=PFc'*PFc;
ShowAbsAngle(Tc)
%% 2 stuff
w1=sin(linspace(0,3*pi,nTraj)).'+1i*sin(pi/4+linspace(0,pi,nTraj)).';
w2=sin(linspace(0,pi,nTraj)).'+1i*sin(linspace(0,2*pi,nTraj)).';

PFc1=RepDotMult(w1,PF1);
PFc2=RepDotMult(w2,PF1);
Tc12=PFc2'*PFc1;
ShowAbsAngle(Tc12)
%%
M=Tc12;
MM=M;
for i=1:N1
    CurIdx=(1:N1)+(i-1)*N1;
    MM(CurIdx,:)=circshift(MM(CurIdx,:),-(i-1)*N1,2);
end
% ShowAbsAngle(MM)
for i=1:N1
    CurIdx=(1:N1)+(i-1)*N1;
    CurCol=MM(:,CurIdx);
    CurCol3=reshape(CurCol.',[N1 N1 N1]);
    Test1=CurCol3(:,:,1:end-(i-1))-repmat(CurCol3(:,:,1),[1 1 N1-i+1]);
    if(isempty(Test1))
        Test1Scr=0;
    else
        Test1Scr=grmss(Test1)./grmss(CurCol3(:,:,1));
    end
    Test2=CurCol3(:,:,end-(i-2):end)-repmat(CurCol3(:,:,end),[1 1 i-1]);
    if(isempty(Test2))
        Test2Scr=0;
    else
        Test2Scr=grmss(Test2)./grmss(CurCol3(:,:,end));
    end
    Test1Scrs(i)=Test1Scr;
    Test2Scrs(i)=Test2Scr;
end
[Test1Scrs; Test2Scrs]
%%
N1=6;
nTraj=13;

Sz128=[N1,N1];
Trajm2=(rand(2,nTraj)-0.5)*N1;

N2=N1;
Nd=[N1 N2];
% [G1 G2] = Gnufft_gram_setup2(arg);
[G1,st1] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Trajm2m=Trajm2;
Trajm2m(1,:)=-Trajm2m(1,:);
[G2,st2] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2m,Sz128),Sz128);

% block1 = reshape(reuse.G1' * real(wi), Nd); % kludge
% block2 = reshape(reuse.G2' * real(wi), Nd);
wi=ones(nTraj,1);
% block1 = reshape(G1' * real(wi), Nd); % kludge
% block2 = reshape(G2' * real(wi), Nd);
block1 = G1' * real(wi); % kludge
block2 = G2' * real(wi);


z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
	[block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
	zeros(1,2*N2);
	[flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
]; % [(2Nd)]
% kern = dsft_gram_hermitify(kern, true); % force Hermitian symmetry
fftkern = fftn(kern);
%%
% FirstRow=T1(1,:);
% FirstRowLastBlock=fliplr(T1(end-N1+1,:));
% 
% h11=reshape(FirstRow,[N1,N1]).';
% h12=conj(reshape(FirstRowLastBlock,[N1,N1]));
% h12(:,end)=0;
% h12=circshift(h12,1,2);
% h11x=h11;
% h11x(1,:)=0;
% h11x(:,1)=0;
% h2=[h11 h12; h12' h11x'];
% % h2=[h11 h12; h12' h11x.'];
% % h2=[h11 h12; h12' h11x];
% % h2=[h11 h12; h12' conj(h11x)];
% H2=fft2(h2);
% % fgmontage(h2)
H2=fftkern;
v=rand(N1,N1);
v=v+v';
Tv=T1*v(:);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=H2.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N1);

Tv=Tm*v;
GGv=G1'*(G1*v);
figure;plot(abs(Tv(:)),'r');hold on;plot(abs(CFHFPadv(:)),'g.');
plot(abs(GGv(:)),'m.');
%%
rmpath('/media/g/aa384808-2266-4edb-87a7-637bc772dc2f/SPENOnline/CS/ESPIRiT/nufft_files')

N1=50;
nTraj=31;
N2=N1;
Nd=[N1 N2];
Sz128=[N1,N1];
Trajm2=(rand(2,nTraj)-0.5)*N1;

% ig = image_geom('nx', N1, 'ny', N2, 'dx', 1, 'offsets', 'dsp');
ig = image_geom('nx', N1, 'ny', N2, 'dx', 1);

% ig.mask = ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
N = ig.dim;
ig.mask = true(N);

% [kspace, omega, wi] = mri_trajectory('spiral0', {}, N, ig.fov, {'voronoi'});
% warn 'todo: crudely fixing dcf at edge of spiral'
% wi(350:end) = wi(350);

kspace=BART2Fes_NUFT_Idxs(Trajm2,Sz128);
% nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
nufft_args = {N, [6 6], 2*N, N/2, 'minmax:kb'};

% Gm = Gmri(kspace/2/pi, ig.mask, 'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);
Gm = Gmri(kspace/2/pi, ig.mask, 'fov', ig.fov, 'basis', {'dirac'}, 'nufft', nufft_args);
% Gm = Gmri(kspace, ig.mask, 'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);

Tm = build_gram(Gm, 1);
fftkern=Tm.fftkern;
H2=fftkern;

[G1,st1] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);

v=rand(N1,N1);
% v=v+v';

Tv=Tm*v;
GGv=G1'*(G1*v)*N1*N1;

GGv2=Gm.Gnufft'*(Gm.Gnufft*v);
GGv2=Gm'*(Gm*v);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=H2.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N1);

figure;
% plot(angle(Tv(:)),'r');hold on;
% plot(angle(CFHFPadv(:)),'g.');hold on;
% plot(angle(GGv(:)),'m.');
% plot(angle(GGv2(:)),'k.');
plot(real(Tv(:)),'r');hold on;
plot(real(CFHFPadv(:)),'go');hold on;
plot(real(GGv(:)),'m+');
plot(real(GGv2(:)),'k*');
legend({'Tm','CFHFPadv','GGv','GGv2'});

% abs(GGv(:)-CFHFPadv(:))
rDiff=abs(GGv(:)-CFHFPadv(:))./abs(GGv(:));
[max(rDiff) mean(rDiff)]*100
rDiff2=abs(GGv(:)-CFHFPadv(:))/mean(abs(GGv(:)));
[max(rDiff2) mean(rDiff2)]*100
%%
wi=ones(nTraj,1);

Trajm2m=Trajm2;
Trajm2m(1,:)=-Trajm2m(1,:);

st=Gm.Gnufft.st;

% nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128,[5 5],Sz128*2);
% nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), st.Nd, st.Jd, st.Kd); % , [0 0] st.om
% nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), st.Nd, st.Jd, st.Kd); % , [0 0] st.om
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), st.Nd, [6 6], st.Kd); % , [0 0] st.om
nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), st.Nd, [6 6], st.Kd); % , [0 0] st.om

block1=nufft_adj(wi, nufftStruct);
block2=nufft_adj(wi, nufftStruct2);

% [G1,st1] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
% [G2,st2] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2m,Sz128),Sz128);
% st1.n_shift = 0 * st.n_shift;
% st2.n_shift = 0 * st.n_shift;
% 
% block1a = (G1' * real(wi))*N1; % kludge
% block2a = (G2' * real(wi))*N1;
% block1a=gfftshift(block1a,1:2);
% block2a=gfftshift(block2a,1:2);
% kspace2=kspace;
% kspace2(:,1) = -kspace2(:,1);
% Gm2 = Gmri(kspace2/2/pi, ig.mask, 'fov', ig.fov, 'basis', {'dirac'}, 'nufft', nufft_args);

% 
% G1 = Gnufft(ig.mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0]});
% st.om(:,1) = -st.om(:,1);
% G2 = Gnufft(ig.mask, {st.om, st.Nd, st.Jd, st.Kd, [0 0]});

% % Working:
% st = rmfield(st, 'phase_shift'); % trick: eliminate phase shift
% st.n_shift = 0 * st.n_shift;
% % build two Gnufft objects, one with negative om1.
% G1 = Gnufft(ig.mask, st);
% st.om(:,1) = -st.om(:,1);
% G2 = Gnufft(ig.mask, st);
% block1 = reshape(G1' * real(wi), Nd); % kludge
% block2 = reshape(G2' * real(wi), Nd);


% block1 = reshape(reuse.G1' * real(wi), Nd); % kludge
% block2 = reshape(reuse.G2' * real(wi), Nd);
% block1 = G1' * real(wi); % kludge
% block2 = G2' * real(wi);
% block1 = Gm' * real(wi); % kludge
% block2 = Gm2' * real(wi);
% block1 = reshape(Gm' * real(wi), Nd); % kludge
% block2 = reshape(Gm2' * real(wi), Nd);

z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
	[block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
	zeros(1,2*N2);
	[flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
]; % [(2Nd)]
% kern = dsft_gram_hermitify(kern, true); % force Hermitian symmetry
fftkern = fftn(kern);
%
H2=fftkern;


[G1,st1] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);

v=rand(N1,N1);
% v=v+v';

Tv=Tm*v;
GGv=G1'*(G1*v)*N1*N1;

GGv2=Gm.Gnufft'*(Gm.Gnufft*v);
GGv2=Gm'*(Gm*v);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=H2.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N1);

figure;
% plot(angle(Tv(:)),'r');hold on;
% plot(angle(CFHFPadv(:)),'g.');hold on;
% plot(angle(GGv(:)),'m.');
% plot(angle(GGv2(:)),'k.');
plot(real(Tv(:)),'r');hold on;
plot(real(CFHFPadv(:)),'go');hold on;
plot(real(GGv(:)),'m+');
plot(real(GGv2(:)),'k*');
legend({'Tm','CFHFPadv','GGv','GGv2'});

% abs(GGv(:)-CFHFPadv(:))
rDiff=abs(GGv(:)-CFHFPadv(:))./abs(GGv(:));
[max(rDiff) mean(rDiff)]*100
rDiff2=abs(GGv(:)-CFHFPadv(:))/mean(abs(GGv(:)));
[max(rDiff2) mean(rDiff2)]*100





%%
TrajForNUFT=load('TrajForNUFT.mat');
Sz128=size(TrajForNUFT.SN);
nTraj=size(TrajForNUFT.Trajm2,2);
Trajm2=TrajForNUFT.Trajm2;
N1=Sz128(1);
N2=Sz128(2);
wi=ones(nTraj,1);

Trajm2m=Trajm2;
Trajm2m(1,:)=-Trajm2m(1,:);

nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om

block1=nufft_adj(wi, nufftStruct);
block2=nufft_adj(wi, nufftStruct2);

fftkern=blocksToFFTKern(block1,block2);

fftkern=single(real(fftkern));


v=rand(Sz128);
% v=v+v';
GGv=nufft_adj(nufft(v,nufftStruct), nufftStruct);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N1);

figure;
plot(real(CFHFPadv(:)),'go');hold on;
plot(real(GGv(:)),'m+');
legend({'CFHFPadv','GGv'});

% abs(GGv(:)-CFHFPadv(:))
rDiff=abs(GGv(:)-CFHFPadv(:))./abs(GGv(:));
[max(rDiff) mean(rDiff)]*100
rDiff2=abs(GGv(:)-CFHFPadv(:))/mean(abs(GGv(:)));
[max(rDiff2) mean(rDiff2)]*100
%%
fftkern=single(real(fftkern));
save('fftkern.mat','fftkern')
%%
B0Data=load('B0TS.mat');
% TSBF=B0Data.TSBF.';
TSC=B0Data.TSC;
L=size(TSC,3);
TSBF=GetTSCoeffsByLinear(nTraj,L).';
v=rand(Sz128) + 1i*rand(Sz128);

SN=TrajForNUFT.SN;
P=TrajForNUFT.P;

P=nufftStruct.p.G;
SN=nufftStruct.sn;
Kd=nufftStruct.Kd;

save('TrajForNUFTx.mat','P','SN','Kd','Trajm2')

TSB=TSBF.';

save('B0Datax.mat','TSBF','TSC')
%%
for l=1:L
%     wi=abs(TSBF(:,t)).^2;
    wi=TSBF(l,:).';
    block1=nufft_adj(wi, nufftStruct);
    block2=nufft_adj(wi, nufftStruct2);
    
    z1 = zeros(N1,1);
    z2 = zeros(N1-1,1);
    kern = [
        [block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
        zeros(1,2*N2);
        [flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
        ]; % [(2Nd)]    
    N1 = size(kern,1);
    N2 = size(kern,2);
    n1 = 0:(N1-1);
    n2 = 0:(N2-1);
    [n1, n2] = ndgrid(n1, n2);
    tmp1 = kern(sub2ind([N1 N2], 1+n1, 1+n2));
    m1 = mod(-n1, N1); % circular symmetry
    m2 = mod(-n2, N2);
    tmp2 = kern(sub2ind([N1 N2], 1+m1, 1+m2));
    tmp2 = conj(tmp2);
    kern = (tmp1+tmp2)/2; % force it to be Hermitian
    N1=Sz128(1);
    N2=Sz128(2);

    fftkernTS(:,:,l) = fftn(kern);
end
fftkernTS=single(fftkernTS);
fftkernTS=real(fftkernTS);
save('fftkernTS.mat','fftkernTS')
%% test O(L)
% v=v+v';
vWithTSC=RepDotMult(v,TSC);
% vWithTSCWithSN=RepDotMult(vWithTSC,SN);
vWithTSCWithSN=RepDotMult(vWithTSC,nufftStruct.sn);
vWithTSCWithSNPadded=padarray(vWithTSCWithSN,[N1 N2 0],'post');
FvWithTSCWithSNPadded=fft2(vWithTSCWithSNPadded);
FvWithTSCWithSNPaddedReshape=reshape(FvWithTSCWithSNPadded,[N1*N2*4,L]);
% SigPerTS_DirectBeforeTSB=P*FvWithTSCWithSNPaddedReshape;
SigPerTS_DirectBeforeTSB=nufftStruct.p.G*FvWithTSCWithSNPaddedReshape;
SigPerTS_Direct=SigPerTS_DirectBeforeTSB.*(TSBF.');
for l=1:L
    CurSig=nufft(vWithTSC(:,:,l),nufftStruct);
    SigPerTS(:,l)=CurSig.*(TSBF(l,:).');
end
[grmss(SigPerTS) grmss(SigPerTS_Direct) grmss(SigPerTS-SigPerTS_Direct)]
SigF=sum(SigPerTS,2);
SigWithTSB=RepDotMult(SigF,conj(TSBF.'));
vPerTSDirect=nufftStruct.p.G'*SigWithTSB;
vPerTSDirectR=reshape(vPerTSDirect,[N1*2,N2*2,L]);
FvPerTSDirectR=ifft2(vPerTSDirectR)*N1*N2*4;
FvPerTSDirectRCropped=FvPerTSDirectR(1:N1,1:N2,:);
FvPerTSDirectRCroppedTSC=FvPerTSDirectRCropped.*conj(TSC);
FvPerTSDirectRCroppedTSCSN=RepDotMult(FvPerTSDirectRCroppedTSC,conj(nufftStruct.sn));
vHDirect=sum(FvPerTSDirectRCroppedTSCSN,3);
for l=1:L
%     SigWithTSB=SigF.*conj(TSBF(l,:).');
    Curv=nufft_adj(SigWithTSB(:,l),nufftStruct);
    CurvHwithTSC(:,:,l)=Curv.*conj(TSC(:,:,l));
end
vH=sum(CurvHwithTSC,3);
disp('finished using nufft');
for l=1:L
    Padv=padarray(vWithTSC(:,:,l),[N1 N2],'post');
    FPadv=fft2(Padv);
    HFPadv=fftkernTS(:,:,l).*FPadv;
    FHFPadv=ifft2(HFPadv);
    AfterTPerTS(:,:,l)=FHFPadv(1:N1,1:N1);
end
AfterTCPerTS=AfterTPerTS.*conj(TSC);
vHUsingT=sum(AfterTCPerTS,3);
disp('finished using T');

M1=mean(abs(vHUsingT(:)));
M2=mean(abs(vH(:)));
dM=mean(abs(vH(:)-vHUsingT(:)));
[M1,M2,dM,dM/max(M1,M2),max(M1,M2)/dM,M1/M2]
% [grmss(vH) grmss(vHDirect) grmss(vHUsingT) grmss(vH-vHUsingT)]
% grmss(vH)/grmss(vH-vHUsingT)
% 
% mean(abs(vHUsingT(:)))
% mean(abs(vH(:)))