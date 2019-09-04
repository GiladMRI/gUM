N=128;

Accs=3:10;
sAccs=sqrt(Accs);
nAccs=numel(Accs);
nEchos=10;
CenterSz=8;

CurM=zeros(N,N,nAccs,nEchos);
for i=1:nAccs
    for j=1:nEchos
        tic
        CurM(:,:,i,j)=squeeze(bart(['poisson -Y ' num2str(N) ' -Z ' num2str(N) ' -y ' num2str(sAccs(i)) ' -z ' num2str(sAccs(i)) ' -C ' num2str(CenterSz) ' -s ' num2str(rand*10000)]));
        t=toc
    end
end

CurMI=int8(CurM);
save('PossionMaps3to10_10reps.mat');