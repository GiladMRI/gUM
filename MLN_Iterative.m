N=64;
nIm=1000;
NN=N*N;
Acc=4;
nCh=8;
nData=NN*nCh/Acc;
% AI=S
% XS=I
% S'X'=I'
% X = lsqminnorm(A,B) returns an array X that solves the linear equation AX = B an
I=single(rand(NN,nIm)+1i*rand(NN,nIm));
S=single(rand(nData,nIm)+1i*rand(nData,nIm));
tic
XT = lsqminnorm(S',I');
% N=64, nIm=10000 : 66.719552 sec
% N=64, nIm=1000 : 2.459254 sec
% N=128, nIm=100 : 3.27 sec
% N=128, nIm=400 : 16.85 sec
% N=64, nIm=100, t=0.22156 sec
t=toc;
disp(['N=' num2str(N) ', nIm=' num2str(nIm) ', t=' num2str(t) ' sec']);
tic
[U,E,V]=svd(XT);
t=toc;
disp(['N=' num2str(N) ', X is ' num2str(size(X)) ', t=' num2str(t) ' sec']);
X=XT';
% N=64, X is 8192  4096, t=24.8155 sec
% Suppose we know the TSC part
% X is [NN,nData]
% XT is [nData,NN]
% X is C*Y
% Y is [nTS, nData]
nTS=7;
C=rand(NN,nTS);
tic
Y=lsqminnorm(C,X);
t=toc;
disp(['% Y is ' num2str(size(Y)) ', t=' num2str(t) ' sec']);
% Y is 7  8192, t=0.9194 sec
%% Now iterate: Limit Y to only neighbors, and get C, then get the new Y, ...
% Iterate option 1: to have CY=X
% Iterate option 2: XS=I, CYS=I. Then given Y, solve for C, given C, solve
% for Y, clean, iterate...
% given Y: CYS=I, (S'Y')C'=I'
nIm=10000;
I=single(rand(NN,nIm)+1i*rand(NN,nIm));
S=single(rand(nData,nIm)+1i*rand(nData,nIm));
tic
tmp=(S')*(Y');
CT = lsqminnorm(tmp,I');
t=toc
disp(['% N=' num2str(N) ', nIm=' num2str(nIm) ', solving for CT takes t=' num2str(t) ' sec']);
% N=64, nIm=1000, solving for CT takes t=0.069641 sec
% N=64, nIm=10000, solving for CT takes t=0.68901 sec
% Given C: CYS=I
% or CY=X
%% Altrnative: advance towards Cnew, Ynew, not jumping there directly
%% Forgot F. Here's directly via regridding:
N=64;
nNeighbors=12;
nCh=8;
nDataG=N*N*nNeighbors*nCh;
nI=3000;
SzD=nDataG*nI;
SzD/1e9
SzI=N*N*nI;
SzI/1e9
Sd=load('/media/a/DATA/All32kImWithPhaseComplexSingleX128x128_Sd.mat');Sd=Sd.Sd;
%%
tic
VecD=single(rand(nNeighbors*nCh,nI)+1i*rand(nNeighbors*nCh,nI));
VecA=single(rand(1,nI)+1i*rand(nNeighbors*nCh,nI));
parfor i=1:N*N
Q=VecA/VecD;
end
toc
% For N=64, nNeighbors*nCh=96, nI=3000, single complex, using parfor with 6 workers: 14 sec
%% FFT for N=64 is 1.5 sec
tic
X=rand(N*2,N*2,nTS,nI);
Y=fft2cg(X);
toc

%%
N=64;
nNeighbors=12;
nCh=8;
nDataG=N*N*nNeighbors*nCh;
nI=3000;
SzD=nDataG*nI;
SzD/1e9
SzI=N*N*nI;
SzI/1e9
Sd=load('/media/a/DATA/All32kImWithPhaseComplexSingleX128x128_Sd.mat');Sd=Sd.Sd;

%%
N=64;
SzI=[N N];

load brain_8ch

DATA=crop(DATA,N,N,size(DATA,3));
I=ifft2c(DATA);
Sens=RunESPIRiTForSensMaps(I);
Im=CalcSENSE1f(I,Sens);
ImWithSens=Im.*Sens;
FData=fft2cg(ImWithSens);
FDataP=permute(FData,[1 2 4 3]);
SensP=permute(Sens,[1 2 4 3]);
nCh=size(Sens,3);
%%
% FTM=dftmtx(N);
% 
% BigFTM=zeros(N,N,N,N);
% for i=1:N
%     for j=1:N
%         BigFTM(i,j,:,:)=FTM(:,i)*FTM(j,:);
%     end
% end
% BigFTMC=repmat(BigFTM,[1 1 1 1 nCh]);
% BigFTMC=BigFTMC.*permute(Sens,[1 2 5 4 3]);
%% get data
Ni=2500;
A=load('/media/a/H1/First3kIm128x128MagSinglex.mat');
A=A.First3kIm128x128MagSingle;
A=permute(A,[2 3 1]);
A=A(:,:,getKrandomSamples(3000,Ni));
AR=imresizeBySlices(A,SzI);
disp('ok');
%% add phase
for i=1:Ni
    if(mod(i,100)), disp(i), end
    Out=GenerateRandomSinPhase(N,5,.1);
    AR(:,:,i)=AR(:,:,i).*single(Out);
end
disp('Added phase');
%%
% FAR=fft2cg(AR);
%% get mask
Acc1D=1.3;
CenterSize=4;
Msk=squeeze(bart(['poisson -Y ' num2str(N) ' -Z ' num2str(N) ' -y ' num2str(Acc1D) ' -z ' num2str(Acc1D) ' -v -e -C ' num2str(CenterSize)]));
ActAcc=(N*N)/gsum(Msk);
MskB=Msk>0;
[I1, I2]=find(MskB);
HalfN=N/2;
CurBartTraj=[I1 I2]'-HalfN;
nTrajX=gsum(MskB);
fgmontage(MskB);title(['ActAcc=' num2str(ActAcc)]);
xlabel(Acc1D);
disp('Generated mask');
%% get NMAP
C=linspace(1-HalfN,HalfN,N);
% Generate Idx mat of neighbors
nNeighbors=12;
NMap=NaN(N,N,nNeighbors);
for i=1:N
    for j=1:N
        CurLoc=C([i j]).';
        D=CurBartTraj(1:2,:)-CurLoc;
        R=grmss(D,1);
        [~, Idx]=sort(R);
        NMap(i,j,:)=Idx(1:nNeighbors);
    end
end

nChToUseInNN=8;
NMapC=RepDotMult(ones(size(NMap))*nTrajX,permute(0:nChToUseInNN-1,[1 4 3 2]));
NMapC=NMapC+repmat(NMap,[1 1 1 nChToUseInNN]);
NMapCX=CombineDims(NMapC,[3 4]);
disp('Generated NMap');
%% get Shat
Nx=nCh*nNeighbors;
CAR=AR.*permute(Sens,[1 2 4 3]);
FCAR=fft2cg(CAR);
DataFlat=single(NaN(nTrajX,Ni,nCh));
for i=1:Ni
    for c=1:nCh
        tmp=FCAR(:,:,i,c);
        DataFlat(:,i,c)=tmp(MskB);
    end
end
DataFlatP=permute(DataFlat,[2 1 3]);
DataFlatPX=reshape(DataFlatP,Ni,nTrajX*nCh);

DataM=single(NaN(N,N,Ni,Nx));
for x=1:N
    disp(x);
    for y=1:N
%         disp([x y]);
        CurInd=squeeze(NMapCX(x,y,:));
        for i=1:Ni
            DataM(x,y,i,:)=DataFlatPX(i,CurInd);
        end
    end
end
disp('Regridded training data');
%%
QQ=zeros(N);
QQ(MskB)=1;
x=28;y=17;
CurInd=squeeze(NMap(x,y,:));
for i=1:nNeighbors
    QQ(I1(CurInd(i)),I2(CurInd(i)))=2;
end
QQ(x,y)=3;
fgmontage(QQ);
%%
% I is data (many images, possibly multichannel)
% A is transformation (e.g. partial FT)
% S is AI, i.e. the signal
% Shat is the signal regridded (all the neighbors, channels concatenated per grid location.
% Out partial left inverse is LFR, where F is (inverse) FT
% L is a maps operation
% R is per grid location (different) kernel operation, to several
% 'channels' or 'time-sgements'
% Finally, I=LFRShat
% Shat is DataM
% I is AR
%% Initialize L
Iter=0;
nTS=7;
L=rand(N,N,nTS)-0.5+1i*(rand(N,N,nTS)-0.5);
L=single(L);
R=single(NaN(N,N,Nx,nTS));
RS=single(NaN(N,N,Ni,nTS));

Lerr=single(NaN(N,N));
Rerr=single(NaN(N,N));

ARest=single(NaN(size(AR)));
%% given L
for i=1:3
    Iter=Iter+1;
    disp('Given L');
    disp(datestr(now));
    Lm=(1./L) .* ((abs(L)./vecnorm(L,2,3)).^2);
    LmI=AR.*permute(Lm,[1 2 4 3]);
    FmLmI=fft2cg(LmI);
    % tmp=lsqminnorm(CurA,CurB);
    % X = lsqminnorm(A,B) returns an array X that solves the linear equation AX = B an
    % RS=I
    % S'R'=I'
    for x=1:N
        %     disp(x);
        for y=1:N
            CurA=squeeze(DataM(x,y,:,:));
            CurB=squeeze(FmLmI(x,y,:,:));
            R(x,y,:,:)=lsqminnorm(CurA,CurB);
        end
    end
    
    % for x=1:N
    % %     disp(x);
    %     for y=1:N
    %         CurA=squeeze(DataM(x,y,:,:));
    %         CurB=squeeze(FmLmI(x,y,:,:));
    %         Rerr(x,y)=norm(CurA*squeeze(R(x,y,:,:)) - CurB);
    %     end
    % end
    % MRerr=gmean(Rerr);
    disp(datestr(now));
    % Given R
    disp('Given R');
    disp(datestr(now));
    for x=1:N
        %     disp(x);
        for y=1:N
            CurSh=squeeze(DataM(x,y,:,:));
            CurR=squeeze(squeeze(R(x,y,:,:)));
            RS(x,y,:,:)=squeeze(sum(CurSh.*permute(CurR,[3 1 2]),2));
        end
    end
    disp(datestr(now));
    FRS=ifft2cg(RS); % this is Y
    LFRS=sum(FRS.*permute43(L),4);
    AfterRErr(Iter)=grmss(AR-LFRS);
    disp(['Iter #' num2str(Iter,'%02d') ', AfterRErr = ' num2str(AfterRErr(Iter),'%.6f')]);
    for x=1:N
        %     disp(x);
        for y=1:N
            CurB=squeeze(AR(x,y,:,:));
            CurA=squeeze(FRS(x,y,:,:));
            L(x,y,:)=squeeze(lsqminnorm(CurA,CurB));
        end
    end
    disp(datestr(now));
    LFRS=sum(FRS.*permute43(L),4);
    AfterLErr(Iter)=grmss(AR-LFRS);
    disp(['Iter #' num2str(Iter,'%02d') ', AfterLErr = ' num2str(AfterLErr(Iter),'%.6f')]);
    % for x=1:N
    % %     disp(x);
    %     for y=1:N
    %         CurA=squeeze(AR(x,y,:,:));
    %         CurB=squeeze(FRS(x,y,:,:));
    %         Lerr(x,y)=norm(CurA*(squeeze(L(x,y,:)).') - CurB);
    %     end
    % end
    % MLerr=gmean(Lerr);
    % disp(datestr(now));
    % disp(['MLerr = ' num2str(MLerr,'%.2f')]);
end
%%
ShowAbsAngle(cat(4,AR(:,:,[45 764 1232]),LFRS(:,:,[45 764 1232])))
%%
clear RDataFlat RDataM
for c=1:nCh
    tmp=FData(:,:,c);
    RDataFlat(:,c)=tmp(MskB);
end
RDataFlatVec=Row(RDataFlat);
for x=1:N
    for y=1:N
        CurInd=squeeze(NMapCX(x,y,:));
        RDataM(x,y,:)=RDataFlatVec(CurInd);
    end
end

RS_R=sum(RDataM.*R,3);
FRS_R=ifft2cg(RS_R); % this is Y
LFRS_R=sum(FRS_R.*permute43(L),4);

