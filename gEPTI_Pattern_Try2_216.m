N=216;
PatternLen=72;

HalfN=N/2;
BaseRatio=N/PatternLen;

%%
nCovers=4;
MaxJump=10;

RotMx=30;
RotMx=floor(PatternLen*2/nCovers);

GR=(1+sqrt(5))/2;

RotMx/GR
% RotG=17;
RotG=23;

nPatterns=1000;
% nPatterns=10;

BigSeed=345454;


LToHave=PatternLen;

setenv('TOOLBOX_PATH','/home/a/bart-0.4.04b')
% LToHave=60;
rng(BigSeed)
% GR=(1+sqrt(5))/2;
% GR2=1-1/GR;
%
Seeds=rand(1,nPatterns)*1e7;
clear COrdX
for j=1:nPatterns
    disp(j);
%     Q=bart(['poisson -Y ' num2str(N) ' -y 2 -Z 1 -z 1. -s ' num2str(Seeds(j))]);
    Q=bart(['poisson -Y ' num2str(N) ' -y ' num2str(BaseRatio) ' -Z 1 -z 1. -s ' num2str(Seeds(j))]);
    CurAddI=0;
    AddOrd=[0 kron(1:10,[1 -1])];
    while gsum(Q)<LToHave
        CurAddI=CurAddI+1;
        Q(HalfN+ AddOrd(CurAddI))=1;
    end
    %
    F=find(Q);
    F=F(1:LToHave);
    RD=mod(randperm(LToHave),nCovers)+1;
    
    S=sparse(F,RD,ones(1,LToHave));
    FS=full(S);
    if(size(FS,1)<N)
        FS(N,1)=0;
    end
    %
    clear COrd
    for i=1:nCovers
        COrd{i}=find(FS(:,i));
    end
    COrd=cat(2,COrd{:});
    COrd(:,2:2:end)=flip(COrd(:,2:2:end),1);
    COrdX(j,:)=COrd(:).';
end
COrdX1=COrdX;
%%
COrdX=COrdX1;
gmax(abs(diff(COrdX,1,2)))

for i=1:nPatterns
    CurPtrn=COrdX(i,:);
    CurPtrn=circshift(CurPtrn,-mod(round(RotG*(i-1)),RotMx));
    while (max(abs(diff(CurPtrn)))>MaxJump)
        D=diff(CurPtrn);
        [~,MI]=max(abs(D));
        CurD=D(MI);
        CurPtrn=[CurPtrn(1:MI) CurPtrn(MI)+floor(CurD/2)  CurPtrn( (MI+1):end)];
    end
    COrdX(i,:)=CurPtrn(1:LToHave);
end

gmax(abs(diff(COrdX,1,2)))
figure;plot(COrdX(:,:).','*','LineWidth',1)
%%
save(['Patterns_' num2str(N) '_' num2str(PatternLen) '.mat']);
%%
WritePatternToFile("Pattern4.txt",COrdX,N);
