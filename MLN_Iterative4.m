load('Bart_Acc_test.mat','MskM');
% load('Bart_Acc_testX.mat','TIs');
% for i=1:numel(TIs)
%     TIsX(:,:,i)=double(TIs{i});
% end
Accs=1.4:0.05:3;
WhichAccs=[1 12 19 26 33];

nIters=5;
%%
N=200;
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
Ni=10000;
% A=load('/media/a/H1/First3kIm128x128MagSinglex.mat');
% A=load('/media/a/ec52f4a8-12c2-4a24-9e6f-d65aaa410529/Gilad/First3kIm128x128MagSinglex.mat');
% A=A.First3kIm128x128MagSingle;
% A=load('/media/a/ec52f4a8-12c2-4a24-9e6f-d65aaa410529/Gilad/First10kIm128x128MagSingle.mat');
% A=A.First10kIm128x128MagSingle;
A=load('/media/a/ec52f4a8-12c2-4a24-9e6f-d65aaa410529/Gilad/First10kIm256x256Magint16.mat');
A=A.I;
A=permute(A,[2 3 1]);
A=A(:,:,getKrandomSamples(size(A,3),Ni));
% Normalize
A=single(A);
Mx=max(1,gmax(A,1:2));
A=A./Mx;
if(size(A,1)~=SzI(1))
    AR=single(NaN([SzI Ni]));
    for i=1:Ni
        if(mod(i,1000)==0), disp(i), end
        AR(:,:,i)=imresize(A(:,:,i),SzI);
    end
else
    AR=A;
end
clear A
disp('ok');
%% add phase
for i=1:Ni
    if(mod(i,1000)==0), disp(i), end
    Out=GenerateRandomSinPhase(N,5,.1);
    AR(:,:,i)=AR(:,:,i).*single(Out);
end
disp('Added phase');

FCAR=fft2cg(AR.*permute(Sens,[1 2 4 3]));
disp('FTd');
%%
% FAR=fft2cg(AR);
%% get mask
% Acc1D=1.4;
% CenterSize=4;
% Msk=squeeze(bart(['poisson -Y ' num2str(N) ' -Z ' num2str(N) ' -y ' num2str(Acc1D) ' -z ' num2str(Acc1D) ' -v -e -C ' num2str(CenterSize)]));
% ActAcc=(N*N)/gsum(Msk);
% MskB=Msk>0;
% [I1, I2]=find(MskB);
% HalfN=N/2;
% CurBartTraj=[I1 I2]'-HalfN;
% nTrajX=gsum(MskB);
% fgmontage(MskB);title(['ActAcc=' num2str(ActAcc)]);
% xlabel(Acc1D);
% disp('Generated mask');
HalfN=N/2;
nNeighbors=12;
nChToUseInNN=8;
Nx=nCh*nNeighbors;
nTS=7;
DataM=single(NaN(N,N,Ni,Nx));
DataM(1)=1i*0.00001; % ~ 150GB for N=200, Ni=5000
disp('ok m');
%%
for CurAccI=1:numel(WhichAccs)
    Acc1D=Accs(WhichAccs(CurAccI));
    Msk=MskM(:,:,WhichAccs(CurAccI))>0;
    MskB=Msk>0;
    [I1, I2]=find(MskB);
    CurBartTraj=[I1 I2]'-HalfN;
    ActAcc=(N*N)/gsum(Msk);
    nTrajX=gsum(MskB);
    % fgmontage(MskB);title(['ActAcc=' num2str(ActAcc)]);
    % xlabel(Acc1D);
    disp(CurAccI);
    disp('Generated mask');
    
    % get NMAP
    C=linspace(1-HalfN,HalfN,N);
    % Generate Idx mat of neighbors
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
    
    NMapC=RepDotMult(ones(size(NMap))*nTrajX,permute(0:nChToUseInNN-1,[1 4 3 2]));
    NMapC=NMapC+repmat(NMap,[1 1 1 nChToUseInNN]);
    NMapCX=CombineDims(NMapC,[3 4]);
    disp('Generated NMap');
    % get Shat
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
    
    for x=1:N
        if(mod(x,100)<4)
            disp([datestr(now) ' ' num2str(x)]);
        end
        for y=1:N
            %         disp([x y]);
            CurInd=squeeze(NMapCX(x,y,:));
            DataM(x,y,:,:)=DataFlatPX(:,CurInd);
        end
    end
    disp('Regridded training data');
    %
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
    % Initialize L
    Iter=0;
    L=rand(N,N,nTS)-0.5+1i*(rand(N,N,nTS)-0.5);
    L=single(L);
    R=single(NaN(N,N,Nx,nTS));
    RS=single(NaN(N,N,Ni,nTS));
    
    Lerr=single(NaN(N,N));
    Rerr=single(NaN(N,N));
    
    ARest=single(NaN(size(AR)));
    disp('Initialized L');
    %
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
        % X = lsqminnorm(A,B) returns an array X that solves the linear equation AX = B an
        % RS=I
        % S'R'=I'
        for x=1:N
            for y=1:N
                CurA=squeeze(DataM(x,y,:,:));
                CurB=squeeze(FmLmI(x,y,:,:));
                R(x,y,:,:)=lsqminnorm(CurA,CurB);
            end
        end
        disp(['L2 ' datestr(now)]);
        % Given R
        disp('Given R');
        for x=1:N
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
        CurSFN=['MLN_Iter_RL_Acc' num2str(Acc1D) '_Ni' num2str(Ni) 'Iter' num2str(Iter) '.mat'];
        save(CurSFN,'R','L');
        disp(['Saved ' CurSFN]);
    end
    disp('Finished iterations');
end
%%
SFN=['MLN_Iter_RL_Acc' num2str(Acc1D) '_Ni' num2str(Ni) '.mat'];
save(SFN,'R','L');
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

ZF=sum(ifft2cg(FData.*MskB).*conj(Sens),3);
ShowAbsAngle(cat(4,Im,LFRS_R,ZF),1,[],'Size',[1 3])

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
%% on test images, all accs
% Iter=5;
% Ni=5000;
for CurAccI=1:numel(WhichAccs)
    disp(CurAccI);
    Acc1D=Accs(WhichAccs(CurAccI));
    
    Msk=MskM(:,:,WhichAccs(CurAccI))>0;
    MskB=Msk>0;
    [I1, I2]=find(MskB);
    CurBartTraj=[I1 I2]'-HalfN;
    ActAcc=(N*N)/gsum(Msk);
    nTrajX=gsum(MskB);
    disp('Generated mask');
    
    % get NMAP
    C=linspace(1-HalfN,HalfN,N);
    % Generate Idx mat of neighbors
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
    
    NMapC=RepDotMult(ones(size(NMap))*nTrajX,permute(0:nChToUseInNN-1,[1 4 3 2]));
    NMapC=NMapC+repmat(NMap,[1 1 1 nChToUseInNN]);
    NMapCX=CombineDims(NMapC,[3 4]);
    disp('Generated NMap');
    NMapCXC{CurAccI}=NMapCX;
    
    clear R L
    CurSFN=['MLN_Iter_RL_Acc' num2str(Acc1D) '_Ni' num2str(Ni) 'Iter' num2str(Iter) '.mat'];
    load(CurSFN,'R','L');
    disp('Loaded R,L');
    
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
    disp('Reconed');
    
    MLNI_Rec(:,:,CurAccI)=LFRS_R;
end
%%
for CurAccI=1:numel(WhichAccs)
    LFRS_R=MLNI_Rec(:,:,CurAccI);
    NMLN=double(LFRS_R*mI./mean(abs(LFRS_R(:))));
    [MLN_Scr,MLN_SI]=ssim(abs(NMLN),ref);
    MLN_Scrs(CurAccI)=MLN_Scr;
end
%%
load('Bart_Acc_testX.mat','TIs');
for i=1:numel(TIs)
    TIsX(:,:,i)=double(TIs{i});
end
nTIs=numel(TIs);

load('Bart_Acc_test.mat','SMV','MI','RecNM');
%%
MskS=grmss(Sens,3)>0.1;
for CurAccI=1:numel(WhichAccs)
    disp(CurAccI);
    Acc1D=Accs(WhichAccs(CurAccI));
    
    Msk=MskM(:,:,WhichAccs(CurAccI))>0;
    MskB=Msk>0;
    [I1, I2]=find(MskB);
    CurBartTraj=[I1 I2]'-HalfN;
    ActAcc=(N*N)/gsum(Msk);
    nTrajX=gsum(MskB);
    disp('Generated mask');
    
    NMapCX=NMapCXC{CurAccI};
    
    clear R L
    CurSFN=['MLN_Iter_RL_Acc' num2str(Acc1D) '_Ni' num2str(Ni) 'Iter' num2str(Iter) '.mat'];
    load(CurSFN,'R','L');
    disp('Loaded R,L');
    
    for t=1:nTIs
        disp([datestr(now) ' test image #' num2str(t)]);
        ImN=double(TIs{t}).*MskS;
        ImN=ImN./gmax(abs(ImN));
        ImWithSens=ImN.*Sens;
        FData=fft2cg(ImWithSens);
        
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
        disp('Reconed');
        
        MLNI_RecX(:,:,CurAccI,t)=LFRS_R;
    end
end
disp('Finished recon test images');
%%
for t=1:nTIs
    Im=double(TIs{t}).*MskS;
    Im=Im./gmax(abs(Im));
    ref=double(abs(Im));
    mI=mean(abs(Im(:)));
    for CurAccI=1:numel(WhichAccs)
        LFRS_R=MLNI_RecX(:,:,CurAccI,t);
        NMLN=double(LFRS_R*mI./mean(abs(LFRS_R(:))));

        [MLN_Scr,MLN_SI]=ssim(abs(NMLN),ref);
        MLN_ScrsX(CurAccI,t)=MLN_Scr;
    end
end
disp('ok');
%%
save('MLNI_Recons.mat','MLNI_Rec','MLNI_RecX','MLN_Scrs');