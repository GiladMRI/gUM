BaseHP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/meas_MID867_gBP_ASL_SMS_Spi_TI_VD1_ST15_FID51481/';
BaseHP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/meas_MID885_gBP_ASL_SMS_Spi_TI_VD1_ST13_low_FID51499/';
Mg=load([BaseHP 'Mg.mat']);
Mg=Mg.Mg;
B0Q2=load([BaseHP 'B0Q2.mat']);
B0Q2=B0Q2.B0Q2;

Bins=-500:1:500;
BinsC=(Bins(1:end-1)+Bins(2:end))/2;
MgFlat=Mg(:);
[N,IA,IB]=histcounts(B0Q2(:),Bins);
NW=N*0;
for i=1:numel(N)
    NW(i)=sum(MgFlat(IB==i));
end

NN=N./sum(N);
NWN=NW./sum(NW);
figure;plot(BinsC,NN,'k');hold on;plot(BinsC,NWN,'r');

BinsA=0:500;
BinsAC=(BinsA(1:end-1)+BinsA(2:end))/2;
MgFlat=Mg(:);
[NA,IAA,IBA]=histcounts(abs(B0Q2(:)),BinsA);
NWA=NA*0;
for i=1:numel(NA)
    NWA(i)=sum(MgFlat(IBA==i));
end

NNA=NA./sum(NA);
NWNA=NWA./sum(NWA);
figure;plot(BinsAC,NNA,'k');hold on;plot(BinsAC,NWNA,'r');

cNNA=cumsum(NNA);
cNWNA=cumsum(NWNA);
figure;plot(cNNA,'k');hold on;plot(cNWNA,'r');
AlmostAllHz=find(cNWNA>0.98,1);
NeedEveryxms=1000/AlmostAllHz/6;
TSPerms=1/NeedEveryxms;
%%
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';
nSlices=96;
nEchos=50;
TimePerEcho_ms=0.93;

Full50Echo=zeros(192,192,nSlices,nEchos);
for SliIdx=1:nSlices
    disp(SliIdx);
    % SliIdx=16;
    load([BaseP 'Sli' num2str(SliIdx) '.mat']); % CurSliS
    load([BaseP 'Sli' num2str(SliIdx) '_Sens.mat']); % SensS: [192×192×12 single]
    
    CurSli192=padarray(padarray(CurSliS(:,1:184,:,:),[4,8,0,0],'pre'),[17,0,0,0],'post');
    
    CurSliCombined=sum(CurSli192.*conj(permute43(SensS)),4);
    
    Full50Echo(:,:,SliIdx,:)=CurSliCombined;
end
Full50EchoS=single(Full50Echo);
save([BaseP 'Full3DAllEchos.mat'],'Full50EchoS');
%%
SensSFull=single(zeros(192,192,nSlices,12));
for SliIdx=1:nSlices
    disp(SliIdx);
    clear SensS
    load([BaseP 'Sli' num2str(SliIdx) '_Sens.mat']); % SensS: [192×192×12 single]
    SensSFull(:,:,SliIdx,:)=SensS;
end
save([BaseP 'SensSFull.mat'],'SensSFull');
%% get B0s
Sz=gsize(Full50EchoS,1:2);
HankelTemporalLen=2;
WhichIdxs=1:nEchos;
[HankelMat, HankelizingMat, DeHankelizingMat]=ghankel(numel(WhichIdxs),HankelTemporalLen);
HankelizingMatP=permute(HankelizingMat,[3:4, 1:2]);
DeHankelizingMatP = permute(DeHankelizingMat,[3:4, 1:2]);

H_for=@(x) reshape(sum(x.*HankelizingMatP,3),[Sz size(HankelMat)]);
H_inv=@(x) squeeze(sum(reshape(x,[Sz numel(HankelMat)]).*DeHankelizingMatP,3));
%%
% SliIdx=40;
for SliIdx=1:nSlices
    disp(SliIdx);
    [ U_LLR, s_LLR, V_LLR ] = batch_svd(H_for(squeeze(Full50EchoS(:,:,SliIdx,:))));
    
    R1=V_LLR(:,:,2,1)./V_LLR(:,:,1,1); % R1 is simply the decay
    CurB0Map=-angle(R1)*1e3/2/pi;
    CurT2S=-1./log(abs(R1));
    B0M(:,:,SliIdx)=CurB0Map;
    R1M(:,:,SliIdx)=R1;
    T2SM(:,:,SliIdx)=CurT2S;
end
B0M=B0M/TimePerEcho_ms;
T2SM=T2SM/TimePerEcho_ms;
save([BaseP 'Full3D_B0_T2S_R1.mat'],'B0M','T2SM','R1M');
fgmontage(B0M(:,:,5:3:end-4),[-500 500])
%%
DBx=B0M([2:end end],:,:)-B0M([1 1:end-1],:,:)./permute([1; ones(size(B0M,1)-2,1); 1],[1 2 3]);
DBy=B0M(:,[2:end end],:)-B0M(:,[1 1:end-1],:)./permute([1; ones(size(B0M,2)-2,1); 1],[3 1 2]);
DBz=B0M(:,:,[2:end end])-B0M(:,:,[1 1:end-1])./permute([1; ones(size(B0M,3)-2,1); 1],[2 3 1]);
DB=cat(4,DBx,DBy,DBz);
DBa=grmss(DB,4)*sqrt(3);
%%
fgmontage(DBx(:,:,8:5:end-5),[-50 50])
fgmontage(DBy(:,:,8:5:end-5),[-50 50])
fgmontage(DBz(:,:,8:5:end-5),[-50 50])
fgmontage(DBa(:,:,8:5:end-5),[0 50])

fgmontage(T2SM(:,:,8:5:end-5),[0 100])
% TimePerEcho_ms ?
%% Simple numerical simulation of one 1D vox
nSpins=900;
HzAcrossVox=150;
SpinsB0Hz=linspace(-HzAcrossVox/2,HzAcrossVox/2,nSpins);
nTimePointsSim=1000;
AcqLengthSim_ms=50;
tSim_ms=linspace(0,AcqLengthSim_ms,nTimePointsSim).';
SpinM=exp(1i*2*pi.*SpinsB0Hz.*tSim_ms/1000);
SumM=sum(SpinM,2)./nSpins;
figure;plot(tSim_ms,abs(SumM),'--');xlabel('ms');title([num2str(HzAcrossVox) 'Hz across voxel, 1D']);
hold on;
plot(tSim_ms,real(SumM));xlabel('ms');
Analytical=1./(pi*HzAcrossVox*tSim_ms/1000);
StartIdx=find(Analytical<abs(SumM)*1.01,1);
plot(tSim_ms(StartIdx:end),Analytical(StartIdx:end),'k');
AnalyticalS=sin(pi*HzAcrossVox*tSim_ms/1000)./(pi*HzAcrossVox*tSim_ms/1000);
plot(tSim_ms(2:end),AnalyticalS(2:end),'m--','LineWidth',3);
% grmss(imag(SumM))
%% Energy retained in small slabs
d=3-1;
AbsSum=abs(Full50EchoS(:,:,1:end-d,:)+Full50EchoS(:,:,d+1:end,:));
SumAbs=abs(Full50EchoS(:,:,1:end-d,:))+abs(Full50EchoS(:,:,d+1:end,:));
Retained=AbsSum./SumAbs;
Retained2=AbsSum./SumAbs;
fgmontage(Retained(:,:,5:25:end,5:7:end),[0 1])
%% figure time-segments linear approximation of phase evolution 
figure;
subplot(1,2,1);
StartPhi=0;
EndPhi=StartPhi+2*pi/6;
plot(exp(1i*linspace(StartPhi,EndPhi,100)),'r');hold on;
plot(exp(1i*[StartPhi EndPhi]),'b');
axis square;axis([-1 1 -1 1]*1.1);
legend({'Phase evolution','linear interpolation'});
subplot(1,2,2);
StartPhi=0;
EndPhi=StartPhi+2*pi/2;
plot(exp(1i*linspace(StartPhi,EndPhi,100)),'r');hold on;
plot(exp(1i*[StartPhi EndPhi]),'b');
axis square;axis([-1 1 -1 1]*1.1);