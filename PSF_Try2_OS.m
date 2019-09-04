T2_ms=70;
T2S_ms=50;
B0_Hz=70;
deltat_ms=2.5e-3;
N=192;
HalfN=N/2;
[X, Y]=ndgrid(1:N,1:N);
RR=sqrt((X-HalfN).^2+(Y-HalfN).^2);
NN=N*N;
Sz=[N N];
t_ms=(0:(NN-1))*deltat_ms;

% Non-uniform
WndSz=5;
TableN=2^10;
Jd=[WndSz WndSz];
Jd32=int32(Jd);
Tbl32=int32([TableN TableN]);
Kd32=int32([N N]);

%
Acc=3;
nTraj=floor(NN/Acc);
%%
Idx=(X-1)*N+Y;
A=zeros(N);
A(Idx(:))=exp(-t_ms/T2S_ms);
F=fftshift(fft2(A));
%%
AR=zeros(N);
RP=randperm(NN);
AR(RP)=exp(-t_ms/T2S_ms);
FAR=fftshift(fft2(AR));

Msk=mod(Y,3)==0;
A3=zeros(N);
A3(Msk(:))=exp(-t_ms(1:NN/3)/T2S_ms);
FA3=fftshift(fft2(A3));
%% Cartesian Acc
CurTtl='Cartesian undersampled no decay';
RPerm=randperm(NN);
CurG=zeros(N);
CurG(RPerm(1:nTraj))=1;

ShowCurrentRegriddedPSF;
%%
CurTtl='Cartesian undersampled + decay';
CurG=zeros(N);
CurG(RPerm(1:nTraj))=exp(-t_ms(1:nTraj)/T2S_ms);

ShowCurrentRegriddedPSF;
%% Spiral
CurTtl='Spiral';

VD=0.6;
% VD=1;
nLoops=12;

CurTtl=[CurTtl ' VD=' num2str(VD) ', ' num2str(nLoops) ' loops'];

R=(linspace(0,1,nTraj).^VD)*(N/2);
Phi=linspace(0,2*pi*nLoops,nTraj);
CTraj=R.*exp(1i*Phi);
Traj=[real(CTraj); imag(CTraj)];
%%
nTraj=size(Traj,2);
w=exp(-t_ms(1:nTraj)/T2S_ms).';

nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Traj, Sz),Sz, [6 6], Sz*2,Sz/2); % , [0 0] st.om

SensMBTSC=imresizeBySlices(SensSB,Sz);
SensMB=SensMBTSC;
wCTS=w;
%
P=nufftStruct.p;
SN=nufftStruct.sn;

fftkerns=NUFFT_to_Toep(nufftStruct,w);
H=ifft2(fftkerns);

% fgmontage(fftshift(log(abs(H))))
%
MultSpTensor=@(Sp,x,Sz) reshape(Sp*reshape(x,size(Sp,2),[]),Sz);
fft2os=@(x) fft2(x,size(x,1)*2,size(x,2)*2);
ifft2os=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2]);
fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));

NUFT_forwSMBTS=@(x,SensMBTSC) MultSpTensor(P,fft2osN((x.*SN).*SensMBTSC),[size(P,1),gsize(SensMBTSC,3:4),1,1,gsize(SensMBTSC,6:9)]);
NUFT_adjSMBTS=@(x,SensMBTSC) sum(sum(ifft2osN(MultSpTensor(P',x,[gsize(SensMBTSC,1:2)*2,gsize(SensMBTSC,3:9)])).*conj(SN).*conj(SensMBTSC),3),6);

NUFT_forwSMBCTS=@(x,SensMBTSC,wa) sum(sum(NUFT_forwSMBTS(x,SensMBTSC).*wa,3),6);
NUFT_adjSMBCTS=@(x,SensMBTSC,wb) NUFT_adjSMBTS(x.*wb,SensMBTSC);

NUFT_adj_forwSMBCTS=@(x,SensMBTSC,wa,wb) NUFT_adjSMBCTS(NUFT_forwSMBCTS(x,SensMBTSC,wa),SensMBTSC,wb);

% That's L^2
NUFT_TSMBCTS_L2=@(x,SensMBTSC,fkerns) sum(sum(permute54(sum(ifft2osN(fkerns.*fft2osN(x.*SensMBTSC)),4)).*conj(SensMBTSC),3),6);

% That's O(L) + L^2 Hadamards
NUFT_TSMBCTS_L=@(x,SensMBTSC,fkerns) sum(sum(ifft2osN(permute54(sum(fkerns.*fft2osN(x.*SensMBTSC),4))).*conj(SensMBTSC),3),6);

% Verify operator normalization
ImSize=[gsize(SensMB,1:2) 1 size(SensMB,4)];
DataSize=[size(P,1) gsize(SensMB,3)];
% Op=OpFromFunc(@(x) NUFT_forwSMB(x,SensMB),@(x) NUFT_adjSMB(x,SensMB));
Op=OpFromFunc(@(x) NUFT_forwSMBCTS(x,SensMBTSC,wCTS),@(x) NUFT_adjSMBCTS(x,SensMBTSC,conj(wCTS)));
x = randn(ImSize) + 1j*randn(ImSize);
y = randn(DataSize) + 1j*randn(DataSize);
Ax = Op*x;
Aty = Op'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
disp(['Operator conj test: ' num2str(Out)]);

w1=interp1(1:50,bas(:,1),linspace(1,50,nTraj),'spline').';
w2=interp1(1:50,bas(:,2),linspace(1,50,nTraj),'spline').';

OnlyCenterP=zeros(N);
OnlyCenterP(N/2+1,N/2+1)=1;
SigX=NUFT_forwSMBCTS(OnlyCenterP,SensMBTSC,w1);
SigY=NUFT_forwSMBCTS(OnlyCenterP,SensMBTSC,w2);

ImXX=NUFT_adjSMBCTS(SigX,SensMBTSC,w1);
ImXY=NUFT_adjSMBCTS(SigX,SensMBTSC,w2);
ImYX=NUFT_adjSMBCTS(SigY,SensMBTSC,w1);
ImYY=NUFT_adjSMBCTS(SigY,SensMBTSC,w2);

ImX=cat(3,ImXX,ImXY);
ImY=cat(3,ImYX,ImYY);
Im=cat(4,ImX,ImY);
%
MxVecXX=getRadialMaxVec(abs(ImXX));
MxVecXY=getRadialMaxVec(abs(ImXY));
MxVecYX=getRadialMaxVec(abs(ImYX));
MxVecYY=getRadialMaxVec(abs(ImYY));

CTrajx=Traj(1,:)+1i*Traj(2,:);
figure;
subplot(3,4,1);
plot(Traj(1,:),Traj(2,:));
title(CurTtl);
xlabel(['#Traj=' num2str(nTraj)]);
setXaxis([-100 100]);
setYaxis([-100 100]);
removeTicks;
subplot(3,4,2:4);
plot((angle(CTrajx)+pi)/(2*pi),'r');hold on;plot(abs(CTrajx)/96,'k','LineWidth',2)
setXaxis([1 nTraj]);
subplot(3,4,[5 6 9 10]);
% gmontage(abs(Im),[0 (max(MxVecXX(2:end))+max(MxVecXX(1:end)))/2]);
gmontage(abs(Im),[0 max(MxVecXX(2:end))*1.2]);
% subplot(3,4,[2 3 4 6 7 8 10 11 12]);
subplot(3,4,[7 8 11 12]);
plot(MxVecXX,'b-','LineWidth',2);hold on
plot(MxVecYY,'k-','LineWidth',2);
plot(MxVecXY,'r-','LineWidth',2);
plot(MxVecYX,'g-','LineWidth',2);
setXaxis([0 50]);
legend({'1-1','2-2','1-2','2-1'});
%%
% ImA=cat(4,cat(3,ImXX,ImXY
% AHA_I=NUFT_adjSMBCTS(DataSB.',SensSB,wCTS)*1e7;
AHA_I=NUFT_TSMBCTS_L(OnlyCenterP,SensMBTSC,fftkerns);

ShowAbsAngle(AHA_I)
%%
















%%


TableNB=TableN*WndSz+1;
% kbw = kaiser(TableNB,2.34*WndSz);
kbw=C1d.';
% figure;plot(linspace(-2,2,TableNB),kbw);12

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%% Spiral nShots
CurTtl='Spiral ';
nInnerShots=17;
VD=0.6;
% VD=1;
CurTtl=[CurTtl num2str(nInnerShots) ' shots, VD=' num2str(VD)];
nLoops=4;
nPreShot=nTraj/nInnerShots;
R=(linspace(0,1,nPreShot).^VD)*(N/2);
Phi=linspace(0,2*pi*nLoops,nPreShot);
CTraj=(R.*exp(1i*Phi)).';

ShotPerm=randperm(nInnerShots);
CTrajN=CTraj.*exp(1i*ShotPerm*2*pi/nInnerShots);
CTraj=CTrajN(:).';
Traj=[real(CTraj); imag(CTraj)];
disp(['Made ' CurTtl]);
%%
TableNB=TableN*WndSz+1;
% kbw = kaiser(TableNB,2.34*WndSz);
kbw=C1d.';
% kbw=kbw*0;kbw(2561+(-510:510))=sqrt(10);
% figure;plot(linspace(-2,2,TableNB),kbw);12
w=exp(-t_ms(1:size(Traj,2))/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%% Full circle cartesian grid no decay
CurTtl='Full circle cartesian grid no decay';
CurG=RR<=HalfN;

ShowCurrentRegriddedPSF;
%%
CurTtl='Cartesian undersampled in circle';

CMsk=RR<=HalfN;
RPC=randperm(gsum(CMsk));
CIdx=Idx(CMsk);
CurG=zeros(N);
CurG(CIdx(RPC(1:nTraj)))=exp(-t_ms(1:nTraj)/T2S_ms);

ShowCurrentRegriddedPSF;
%%
CurTtl='Cartesian undersampled in circle';

CMsk=RR<=HalfN;
RPC=randperm(gsum(CMsk));
CIdx=Idx(CMsk);
tmp=zeros(N);
tmp(CIdx(RPC(1:nTraj)))=1;
[I1,I2]=find(tmp);
RP2=randperm(numel(I1));
Traj=[I1,I2].'-HalfN;
Traj=Traj(:,RP2);
w=exp(-t_ms(1:size(Traj,2))/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%%
CurTtl='Random VD radius, phase, VD=';

VDRR=0.5;
CurTtl=[CurTtl num2str(VDRR)];
CTraj=(rand(1,nTraj).^VDRR)*HalfN.*exp(1i*rand(1,nTraj)*2*pi);
Traj=[real(CTraj);imag(CTraj)];
w=exp(-t_ms(1:size(Traj,2))/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%% Realistic spiral
CurTtl='Realistic spiral ';
CurTtl=[CurTtl 'VD=' num2str(VD)];
Traj=BARTTrajxQ; %kTraj.'*FOVx/1000/2/pi;
%%
w=exp(-t_ms(1:size(Traj,2))/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%%
CurTtl='Rosette';
nRShots=4;
kmax=96;
w1=1.4;
w2=40;
t=linspace(0,40,nTraj/nRShots);
CBase=kmax*sin(w1*t).*exp(1i*w2*t);
C=CBase;
for s=1:nRShots
    C=[C CBase*exp(1i*s*2*pi/nRShots)];
end
Traj=[real(C);imag(C)];
% figure;plot(real(C),imag(C),'.')
w=exp(-t_ms(1:size(Traj,2))/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%%
CurTtl='Equispaced undersampled in circle';
nTrajx=nTraj*4/pi;
snTrajx=floor(sqrt(nTrajx));
Ve=linspace(-HalfN,HalfN,snTrajx);
[Xe, Ye]=meshgrid(Ve,Ve);
RRe=sqrt(Xe.^2+Ye.^2);
Mske=RRe<=HalfN;
[I1e, I2e]=find(Mske);
I1e=Ve(I1e);
I2e=Ve(I2e);
Traj=[I1e;I2e];
nTraje=numel(I1e);
RPe=randperm(nTraje);
Traj=Traj(:,RPe);

% kbw=kbw*0;kbw(2561+(-256:256))=sqrt(10);

w=exp(-t_ms(1:size(Traj,2))*0/T2S_ms)/10;

Reg=fftshift(reshape(interp2_table_adj_mex(complex(w.'), kbw, kbw, Jd32, Tbl32, Traj.', Kd32),Sz));
CurG=min(Reg,1);

ShowCurrentRegriddedPSF;
%%
% figure;imagesc(A)
% figure;imagesc(abs(FA3))
figure;plot(getRadialMaxVec(abs(FA3)))
%%
M=abs(FA3);
%%% figure;plot(linspace(-2,2,TableNB),kbw);

Traj=[real(CTraj); imag(CTraj)];
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Traj,Sz), Sz, [5 5], Sz); % , [0 0] st.om
TablsN=2^10;
WndSz=5;
TableNB=TablsN*WndSz+1;
% kbw = kaiser(TableNB,2.34*WndSz);
kbw=C1d.';
figure;plot(linspace(-2,2,TableNB),kbw);
Jd=[WndSz WndSz];
st2=nufft_init(BART2Fes_NUFT_Idxs(Traj,Sz), Sz, Jd, Sz, 'table', TablsN, 'minmax:kb');
st2.sn=ones(Sz);