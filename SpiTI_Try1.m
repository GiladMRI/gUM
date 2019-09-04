% mex VDSpiralMex.cpp VDSpiral.cpp -DMATLAB

FOVx=200;
dFOV=FOVx/1000;
paramLongROSamples=1024*12;
AccR=10.2;
nInnerInterleaves=6;

% paramLongROSamples=1024*15;
% AccR=2.8;

GRAD_RASTER_TIME=10;

spBW=4e5;

TrajPointsPerGrad=GRAD_RASTER_TIME*spBW/1e6;
paramLongInterleaves=1;
VD=1.3;
paramLongSpGradAmp=35;
paramLongSpSlewRate=155;
% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,1]);

% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/2,spBW,AccR,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples/nInnerInterleaves,spBW,AccR,...
    paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,nInnerInterleaves,1]);

dkToGFac=FOVx/1000/2/pi;
BARTTrajx=kTraj.'*FOVx/1000/2/pi;
nTraj=size(BARTTrajx,2);
C=BARTTrajx(1,1:size(BARTTrajx,2)/2)+1i*BARTTrajx(2,1:size(BARTTrajx,2)/2);
%
% figure;plot(BARTTrajx(1,1:end),BARTTrajx(2,1:end),'.')
% figure;plot(BARTTrajx(1,1:end-40),BARTTrajx(2,1:end-40),'b.')
% hold on
% plot(BARTTrajx(1,end-39:end-19),BARTTrajx(2,end-39:end-19),'k.')
% plot(BARTTrajx(1,end-18:end),BARTTrajx(2,end-18:end),'r.')
%
MaxK=max(BARTTrajx(:));
Acc=ceil(MaxK*2).^2/nTraj;

figure;
% subplot(2,2,1);
plot3(1:size(BARTTrajx,2),BARTTrajx(1,:),BARTTrajx(2,:))
%
HnTraj=nTraj/nInnerInterleaves;
CLRs='rgbcmkrg';
% subplot(2,2,2);
plot(BARTTrajx(1,1:HnTraj),BARTTrajx(2,1:HnTraj),[CLRs(1) '-'],'LineWidth',2)
hold on
plot(BARTTrajx(1,1:HnTraj),BARTTrajx(2,1:HnTraj),[CLRs(1) 'o'],'LineWidth',2)
for i=2:nInnerInterleaves
    plot(BARTTrajx(1,HnTraj*(i-1)+(1:HnTraj)),BARTTrajx(2,HnTraj*(i-1)+(1:HnTraj)),[CLRs(i) '-'],'LineWidth',2)
    plot(BARTTrajx(1,HnTraj*(i-1)+(1:HnTraj)),BARTTrajx(2,HnTraj*(i-1)+(1:HnTraj)),[CLRs(i) 'o'],'LineWidth',1)
end
axis equal
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
%%
subplot(2,2,3);
Ts=(0:(size(GradBuf,1)-1))*10;
plot(Ts,GradBuf*MaxGrad*1000,'.-');title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
% setXaxis([1 nTraj]);
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel('time (ms)');
ylabel('mT/m');
%
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(Ts(1:end-1),SlewBuf*100);
MaxSlew=max(max(abs(SlewBuf(20:HnTraj,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
% setXaxis([1 nTraj]);
setXaxis([0 Ts(end)]);
set(gca,'XTick',10000:10000:70000);
XT=get(gca,'XTick');
set(gca,'XTickLabel',XT/1000)
xlabel('time (ms)');
ylabel('mT/m/s');
%%
EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;

clear kTrajQ
kTrajQ(:,1) = interp1(1:size(kTraj,1),kTraj(:,1),1:1e5/spBW:(size(kTraj,1)-0.01));
kTrajQ(:,2) = interp1(1:size(kTraj,1),kTraj(:,2),1:1e5/spBW:(size(kTraj,1)-0.01));

BARTTrajx=kTrajQ.'*FOVx/1000/2/pi;
BARTTrajx(3,end)=0;

BARTTrajMS=BARTTrajx;
BARTTrajxC=BARTTrajx(1,:)+1i*BARTTrajx(2,:);
for i=2:paramLongInterleaves
    CurRotTrajC=BARTTrajxC*exp(1i*2*pi*(i-1)/paramLongInterleaves);
    BARTTrajMS(1,:,i)=real(CurRotTrajC);
    BARTTrajMS(2,:,i)=imag(CurRotTrajC);
end
BARTTrajMS(3,end)=0;
disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
XResmm=dFOV*1000/(2*MaxK);
%
Trajm2=BARTTrajMS(1:2,1:end-2);
nTraj=size(Trajm2,2);

Sz128=[128 128];

[FesNUFTOp,st] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2,Sz128),Sz128);
Kd=st.nufftStruct.Kd;
SN=st.nufftStruct.sn;
try
P=st.nufftStruct.p/sqrt(prod(Sz128));
catch
P=st.nufftStruct.p.G/sqrt(prod(Sz128));
end

Trajm2m=Trajm2;
Trajm2m(1,:)=-Trajm2m(1,:);

[FesNUFTOpm,stm] = nuFTOperator(BART2Fes_NUFT_Idxs(Trajm2m,Sz128),Sz128);
% Pm=stm.nufftStruct.p/sqrt(prod(Sz128));
Pm=stm.nufftStruct.p.G/sqrt(prod(Sz128));
%%
save(['TrajForNUFTm_6TI.mat'],'Trajm2','SN','Kd','P','Pm');
disp('Saved TrajForNUFTm_6TI.mat');
%%
QQ=reshape(full(P(100,:)),[256 256]);
B=reshape(QQ(abs(QQ)>0),[5 5]);
%%
wi=ones(nTraj,1);
block1=FesNUFTOp'*wi;
block2=FesNUFTOpm'*wi;

fftkern1=blocksToFFTKern(block1,block2);
fftkern = single(real(fftkern));
disp('fftnern ok');
%% Test
I1=imresize(cameraman,Sz);
I1=I1./gmax(I1);
I1(1)=I1(1)+1i*0.00001;
%%
nufftStruct = nufft_init(BART2Fes_NUFT_Idxs(Trajm2,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om
nufftStruct2 = nufft_init(BART2Fes_NUFT_Idxs(Trajm2m,Sz128), Sz128, [6 6], Sz128*2); % , [0 0] st.om

block1=nufft_adj(wi, nufftStruct);
block2=nufft_adj(wi, nufftStruct2);

fftkern=blocksToFFTKern(block1,block2);

fftkern=single(real(fftkern));


v=I1;
% v=v+v';
GGv=nufft_adj(nufft(v,nufftStruct), nufftStruct);

Padv=v;
Padv(N1*2,N1*2)=0;
FPadv=fft2(Padv);
HFPadv=fftkern.*FPadv;
FHFPadv=ifft2(HFPadv);
CFHFPadv=FHFPadv(1:N1,1:N2);

rDiff=abs(GGv(:)-CFHFPadv(:))./abs(GGv(:));
[max(rDiff) mean(rDiff)]*100
rDiff2=abs(GGv(:)-CFHFPadv(:))/mean(abs(GGv(:)));
[max(rDiff2) mean(rDiff2)]*100
%%
AI=FesNUFTOp*I1;
AHAI=FesNUFTOp'*AI;

PaddedI=padarray(I1,size(I1),0,'post');
% FI=fftn(I1,size(I1)*2);
FI=fft2(PaddedI);
KFI=FI.*fftkern1;
IFKFI=ifft2(KFI);
TI=IFKFI(1:N1,1:N2);

%%
N=128;
BTT=[-12.2,5.6];
TT=BART2Fes_NUFT_Idxs(BTT.',Sz128);
% nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
% nufft_args = {N, [6 6], 2*N, N/2, 'minmax:kb'};
Jd=[5 5];
st1=nufft_init(TT, Sz128, Jd, Sz128*2, 'minmax:kb');
TablsN=2^16;
st2=nufft_init(TT, Sz128, Jd, Sz128*2, 'table', TablsN, 'minmax:kb');
B=reshape(full(st1.p(abs(st1.p)>0)),[5 5]);
% v1=flipud( (st2.h{1}(   + (1:TablsN:Jd(1)*TablsN)))).';
% v1=flipud( (st2.h{1}(  floor(( (BTT(1)-floor(BTT(1)))*2+0.5)*TablsN)+ (1:TablsN:Jd(1)*TablsN)))).';
v1=flipud( (st2.h{1}(  mod(floor(( (BTT(1)-floor(BTT(1)))*2+0.5)*TablsN)+6*TablsN-1,TablsN)+1+ (1:TablsN:Jd(1)*TablsN)))).';
[(v1./norm(v1)); B(:,3).'./norm(B(:,3))]
v2=flipud( (st2.h{2}(  mod(floor(( (BTT(2)-floor(BTT(2)))*2+0.5)*TablsN)+6*TablsN-1,TablsN)+1+ (1:TablsN:Jd(2)*TablsN)))).';
% v2=conj(st2.h{2}((BTT(2)-floor(BTT(2))+0.5)*TablsN+ (1:TablsN:Jd(2)*TablsN))).';
% (v2./norm(v2))
([(v2./norm(v2)); B(3,:)./norm(B(3,:))])

Best=v1.'*v2;
abs(v1.'*v2)
abs(B)

angle(v1.'*v2)
angle(B)

[grmss(B) grmss(Best) grmss(Best-B)]

Pfull=reshape(full(st1.p),[256 256]);
[I1 I2]=find(abs(Pfull)>0);
[min(I1) min(I2) ]

Pest=zeros(256);
StartIs=mod(round(BTT*2)-2 +256,256)+1;
for i=1:Jd(1)
    for j=1:Jd(2)
        i1=mod(StartIs(1)+i-2+256*6,256)+1;
        j1=mod(StartIs(2)+j-2+256*6,256)+1;
        Pest(i1,j1)=Best(i,j);
    end
end

[grmss(Pfull) grmss(Pest) grmss(Pfull-Pest)]
