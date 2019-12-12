%% Eig test
N=16000;
R=4e4;
B=rand(R,N)+1i*rand(R,N);
A=B'*B;
disp(['Starting ' datestr(now)]);
tic
[V,D]=eig(A);
t=toc;
disp(['N=' num2str(N), ' eig time=' num2str(t) ' sec']);

% N=1000 eig time=0.13246 sec
% N=2000 eig time=1.1052 sec
% N=4000 eig time=9.8171 sec
% N=8000 eig time=69.8039 sec
% N=16000 eig time=623.1064 sec
%% SVD Test
nP=4000;
N=nP*2;
X=rand(N,nP)+1i*rand(N,nP);
X=single(X);
mX=mean(X,1);
X=X-mX;

disp(['Starting ' datestr(now)]);
tic
[U,S,V]=svd(X,'econ');
t=toc;
disp(['N,R=' num2str(nP) ',' num2str(N) ', eig time=' num2str(t) ' sec']);

% Double:
% N,R=1000,2000, eig time=0.70437 sec
% N,R=2000,4000, eig time=4.8478 sec
% N,R=4000,8000, eig time=38.0483 sec
% N,R=8000,16000, eig time=330.5892 sec
% Single:
% N,R=1000,2000, eig time=0.28412 sec
% N,R=2000,4000, eig time=2.6184 sec
% N,R=4000,8000, eig time=19.497 sec

%%
AllImWithPhaseComplexSingle=load('/media/a/H1/All32kImWithPhaseComplexSingleX128x128.mat');AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle.AllImWithPhaseComplexSingle;
AllImWithPhaseComplexSingle=reshape(AllImWithPhaseComplexSingle,[31248,128*128]);
mX=mean(AllImWithPhaseComplexSingle,1);
AllImWithPhaseComplexSingle=AllImWithPhaseComplexSingle-mX;

disp(['Starting ' datestr(now)]);
[U,S,V]=svd(AllImWithPhaseComplexSingle,'econ');
disp(['Finished ' datestr(now)]);
Sd=diag(S);
save('/media/a/H1/All32kImWithPhaseComplexSingleX128x128_mX.mat','mX')
save('/media/a/H1/All32kImWithPhaseComplexSingleX128x128_Sd.mat','Sd')
save('/media/a/H1/All32kImWithPhaseComplexSingleX128x128_V.mat','-v7.3','V')
%% SVD Eig Test
n=20;
p=3;
X=rand(n,p)+1i*rand(n,p);
mX=mean(X,1);
X=X-mX;
C=X'*X/(n-1);
[U,S,V]=svd(X,'econ');
[VV,DD]=eig(X'*X,'vector');

[d,ind] = sort(DD,'descend');
DDs = DD(ind);
VVs = VV(:,ind);
[sqrt(d) diag(S)]
% VVs
% V
[abs(V); abs(VVs);]