nufftStructT = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), st.Nd, [6 6], st.Kd,'table', 2^10, 'minmax:kb'); % , [0 0] st.om
%%
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
%%
v=imresize(cameraman,Sz128);
v(1)=v(1)+1i*0.001;
% v=rand(Sz128);
% v=v+v';v
Gv=nufft(v,nufftStruct);
GGv=nufft_adj(nufft(v,nufftStruct), nufftStruct);
fgmontage(GGv)
%%
Trajm2Base=Trajm2;
%%
r=linspace(0,30,nTraj);
Phi=linspace(0,2*pi*35,nTraj);
Trajm2C=r.*exp(1i*Phi);
Trajm2(1,:)=real(Trajm2C);
Trajm2(2,:)=imag(Trajm2C);
figure;plot(Trajm2(1,:),Trajm2(2,:),'.');grid on;set(gca,'XTick',-50:50);set(gca,'YTick',-50:50)
%%
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
Gv=nufft(v,nufftStruct);

%%
[X Y]=ndgrid(-63:64);
% [X Y]=ndgrid(-30:31,-30:31);
vq = griddata(Trajm2(1,:),Trajm2(2,:),Gv,X,Y, 'cubic');
vq(~isfinite(vq))=0;
fgmontage(fftshift(ifft2cg(vq)))

%%
nNeighnors=12;
C=-64:63;
osN=128;
NMap=zeros(osN,osN,nNeighnors);
RMap=zeros(osN,osN,nNeighnors);
DMap=zeros(osN,osN,nNeighnors,2);
for i=1:osN
    for j=1:osN
        CurLoc=[C(i), C(j)];
        D=Trajm2.'-repmat(CurLoc,[nTraj 1]);
        R=grmss(D,2);
        D1=D(:,1);
        D2=D(:,2);
        [RSort,Idx]=sort(R,'ascend');
        NIdx=Idx(1:nNeighnors);
        NMap(i,j,:)=NIdx;
        RMap(i,j,:)=R(NIdx);
        DMap(i,j,:,1)=D1(NIdx);
        DMap(i,j,:,2)=D2(NIdx);
    end
end
%%
SigM=Gv(NMap);
%%
Table=nufftStructT.h{1};
% Table=abs(Table);
% Table=real(Table);
DMapx=DMap*nufftStructT.Ld(1)*2+(numel(Table)+1)/2;
DMapx=max(1,min(numel(Table),round(DMapx)));
WMap=Table(DMapx);
WMap=WMap(:,:,:,1).*WMap(:,:,:,2);
WMap(max(abs(DMap),[],4)>1.5)=0;
% WnMap=sum(WMap,3);
WnMap=grmss(WMap,3);
NWMap=RepDotMult(WMap,1./WnMap);
NWMap(~isfinite(NWMap))=0;

Reg=sum(SigM.*NWMap,3);
FReg=ifft2cg(Reg);
FReg=fftshift(FReg);
% fgmontage(FReg)
ShowAbsAngle(FReg)
%%