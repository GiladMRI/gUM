N=120;
HalfN=N/2;
% get image
I=phantom(N);
I=imresize(cameraman,[N N]);
I=cameraman;
% Traj
Acc=2;
nTraj=floor(N*N/Acc);
R=linspace(0,HalfN/(1),nTraj);
nLoops=70;
Phi=linspace(0,2*pi*nLoops,nTraj);
C=R.*exp(1i*Phi);
kx=real(C);
ky=imag(C);
% figure;plot(kx,ky,'-o');
% Apply NUFFT
BARTTraj=[kx;ky];
BARTTraj(3,1)=0;
V=bart('nufft',BARTTraj,I);
% V=bart('nufft',BARTTraj,Ix.*SN);
% BART NUFFT recon
BR=bart('nufft -i -l 0.2',BARTTraj,V);
% BART PICS recon
BRP=bart('pics -S -m -R T:3:0:1e-4 -t ',permute(BARTTraj,[1 2 3 4]),permute(V,[3 2 4 1]),ones(N));
ShowAbsAngle(BRP)
%%
% Find Neighbors
disp('Find Neighbors');
nNeighbors=12;
clear NearIdxs Dists
FLocs=(-HalfN+1):HalfN;
[X Y]=ndgrid(FLocs,FLocs);
RR=sqrt(X.^2+Y.^2);
for i=1:N
    for j=1:N
        CurLoc=FLocs(i)+1i*FLocs(j);
        D=C-CurLoc;
        [~,Ord]=sort(abs(D),'ascend');
        NearIdxs(i,j,:)=Ord(1:nNeighbors);
        Dists(i,j,:)=D(Ord(1:nNeighbors));
    end
end
% Density compensation
disp('Density compensation');
[Vd,Cd] = voronoin([kx.' ky.']);
Vdx=min(Vd,HalfN);
DC=NaN(1,nTraj);
for i=1:nTraj
    CurX=Vdx(Cd{i},1);
    CurY=Vdx(Cd{i},2);
    DC(i)=polyarea(CurX,CurY);
%     DC(i)=polyarea(max(min(CurX,HalfN),-HalfN),max(min(CurY,HalfN),-HalfN));
end
DCM=NaN(N,N,nNeighbors);
for i=1:N
    for j=1:N
        DCM(i,j,:)=DC(squeeze(NearIdxs(i,j,:)));
    end
end
% get Coeffs
Coeffs=Dists*0;
KB_J=6;
for i=1:N
    if(mod(i,10)==1)
        disp(['Coeffs ' num2str(i)]);
    end
    for j=1:N
%         Coeffs(i,j,:)=kaiser_bessel(abs(squeeze(Dists(i,j,:))), 8);
        Coeffs(i,j,:)=kaiser_bessel(real(squeeze(Dists(i,j,:))), KB_J).*kaiser_bessel(imag(squeeze(Dists(i,j,:))), KB_J);
    end
end

CoeffsDC=Coeffs.*DCM;
SCoeffs=sum(CoeffsDC,3);
% Regrid
ValsM=NaN(N,N,nNeighbors);
for i=1:N
    for j=1:N
        CurVals=V(NearIdxs(i,j,:));
        ValsM(i,j,:)=CurVals;
    end
end

RG=sum(ValsM.*CoeffsDC,3);

SCoeffsN=SCoeffs;
% SCoeffsN(SCoeffs<10)=1;
RGN=RG./SCoeffsN;
RGN(RR>(HalfN-5))=0;
% Recon
RC=ifft2cg(RGN);
ShowAbsAngle(RC)
%%
a=kaiser_bessel(linspace(-HalfN,HalfN,N), KB_J);
Z=1./abs(fft1cg(a,2));
ZZ=Z'*Z;
%%
a=10;
A=1:(a*2+1);
F1=fft1cg(A,2);
BT=(-a:a)-0.5;
BT(2,:)=BT;BT(1,:)=-0.5;
BT(3,1)=0;
A3=repmat(A,[3 1]);
A3([1 3],:)=0;
BF=bart('nufft',BT,A3)*sqrt(3);
abs([F1;BF])
diff(angle([F1;BF])*pi,[],1)