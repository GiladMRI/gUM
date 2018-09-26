% figure;imagesc(I,[0 1]);
%
nPhantoms=30000;
AllPhantomWithPhaseComplexSingle=single(zeros(nPhantoms,128,128));
AllPhantomWithPhaseComplexSingle(1)=0+0.0001i;

Ne=100;
Val=rand(Ne,nPhantoms)*10-5;
% WH=randn(Ne,2)*2-1;
WH=randn(Ne,2,nPhantoms)/4;
XY=rand(Ne,2,nPhantoms)*2-1;
Phi=rand(Ne,nPhantoms)*360;

for i=1:30000
    if(mod(i,100)==1)
        disp(i);
    end
    CurE=[Val(:,i) WH(:,:,i) XY(:,:,i) Phi(:,i)];
    I=single(min(1,abs(phantom(CurE,128))/15));
    AllPhantomWithPhaseComplexSingle(i,:,:)=I.*squeeze(exp(1i*angle(AllImWithPhaseComplexSingle(i,:,:))));
end
%%
save('/media/a/H1/AllPhantomWithPhaseComplexSingle128x128.mat','-v7.3','AllPhantomWithPhaseComplexSingle');