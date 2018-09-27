N=128;
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
disp('ok');
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
% A=load('/media/a/H1/First3kIm128x128MagSinglex.mat');
A=load('/media/a/ec52f4a8-12c2-4a24-9e6f-d65aaa410529/Gilad/First3kIm128x128MagSinglex.mat');
A=A.First3kIm128x128MagSingle;
A=permute(A,[2 3 1]);
A=A(:,:,getKrandomSamples(3000,Ni));
if(size(A,1)~=SzI(1))
    AR=imresizeBySlices(A,SzI);
else
    AR=A;
end
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
Acc1D=1.4;
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
disp('Takes a few seconds..');
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
    disp([datestr(now) ' ' num2str(x)]);
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
disp('Initialized L');
%%
nIters=5;
% With N=128, Ni=2500, 12 Neighbors, 7TS, it took ~9min per iteration (on the server)
% Total memory used was <50GB
for i=1:nIters
    Iter=Iter+1;
    disp('Given L');
    StartIter=now;
    disp(datestr(now));
    % given L
    % This part 2sec
    Lm=(1./L) .* ((abs(L)./vecnorm(L,2,3)).^2);
    LmI=AR.*permute(Lm,[1 2 4 3]);
    FmLmI=fft2cg(LmI);
    disp(['L1 ' datestr(now)]);
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
    disp(['L2 ' datestr(now)]);
    % for x=1:N
    % %     disp(x);
    %     for y=1:N
    %         CurA=squeeze(DataM(x,y,:,:));
    %         CurB=squeeze(FmLmI(x,y,:,:));
    %         Rerr(x,y)=norm(CurA*squeeze(R(x,y,:,:)) - CurB);
    %     end
    % end
    % MRerr=gmean(Rerr);
    % Given R
    disp('Given R');
    for x=1:N
        %     disp(x);
        for y=1:N
            CurSh=squeeze(DataM(x,y,:,:));
            CurR=squeeze(squeeze(R(x,y,:,:)));
            RS(x,y,:,:)=squeeze(sum(CurSh.*permute(CurR,[3 1 2]),2));
        end
    end
    disp(['R1 ' datestr(now)]);
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
    disp(['R2 ' datestr(now)]);
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
disp('Finished iterations');
%% See results on training
ShowAbsAngle(cat(4,AR(:,:,[45 764 1232]),LFRS(:,:,[45 764 1232])))
%% This part reorders the data
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
% The recon using R and L
RS_R=sum(RDataM.*R,3);
FRS_R=ifft2cg(RS_R); % this is Y
LFRS_R=sum(FRS_R.*permute43(L),4);

ShowAbsAngle(cat(4,Im,LFRS_R))

ref=double(abs(Im));
mI=mean(abs(Im(:)));

NMLN=double(LFRS_R*mI./mean(abs(LFRS_R(:))));

[MLN_Scr,MLN_SI]=ssim(abs(NMLN),ref);
title(MLN_Scr);
%%
% ShowAbsAngle(L)

R3=repmat(grmss(R,3:4),[1 1 3])*2;
Clr=[1 0 0];
for i=1:3
    tmp=R3(:,:,i);
    tmp(MskB)=Clr(i);
    R3(:,:,i)=tmp;
end
figure;imshow(R3)