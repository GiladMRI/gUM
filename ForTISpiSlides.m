ScanP='/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/';
BaseFN='meas_MID01088_FID09952_gSpi2d_T12_d110_Dw11_40rep_VD123';
% BaseFN='meas_MID01106_FID09970_gSpi2d_T12_d110_Dw11_40rep_VD123_Thin';

TrajType=str2num(BaseFN(32:33));
ResType=floor((TrajType-10)/2)+1; % 2,1.5,1
TimingType=mod(TrajType-10,2)+1; % 5ms, 10ms
load('GAll5ms10ms.mat');
GTrajaCBase=GAll(:,TimingType,ResType);
if(TimingType==1)
    nInnerShots=10;
else
    nInnerShots=5;
end
ResStrs={'xxx','1.9mm','1.3mm'};
InnerShotLengthStrs={'5ms','10ms'};

figure; plot([real(GTrajaCBase) imag(GTrajaCBase)],'LineWidth',2); removeTicks; setXaxis([-100 numel(GTrajaCBase)+100])
%%
GradReduceFac=1.23;
GTrajaC=GTrajaCBase/GradReduceFac;

gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_us=10;
GradDwellTime_ms=GradDwellTime_us/1000;

g=GTrajaC;
k=cumsum([0; g])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
s=diff(g)/GradDwellTime_ms;

FOV_mm=220;
FOVx=220;

kK=k*FOV_mm/1000/2/pi;

Kmax=ceil(max(abs(kK)));
NTrg=Kmax*2;
res_mm=FOVx/(max(abs(kK))*2);

figure;plot(kK);axis([-80 80 -80 80]);

figure;plot3(1:numel(kK),real(kK),imag(kK));
%%
NK=1000;
kKs=interp1(1:numel(kK),kK,linspace(1,numel(kK),NK));
figure;
for i=1:NK-1
%     plot3(i,real(kKs(i)),imag(kKs(i)),'.','Color',hsv2rgb([i/NK,1,1]));
%     plot3([i i+1],real(kKs([i i+1])),imag(kKs([i i+1])),'Color',hsv2rgb([i/NK,1,1]));
    plot(real(kKs([i i+1])),imag(kKs([i i+1])),'Color',hsv2rgb([i/NK,1,1]));
    hold on;
end
axis([-80 80 -80 80]);
%%
figure;plot(abs(kK));axis([-300 numel(kK)+300 -10 Kmax+10]);
%%
figure;
n=4;
Len=2000;
StIdxs=ceil(linspace(1,numel(kK)-Len-1,n));
for i=1:n
    StIdx=StIdxs(i);
    CurIdxs=StIdx+(0:Len);
    plot(real(kK(CurIdxs))+(i-1)*150,imag(kK(CurIdxs)));hold on;
end
setYaxis([-80 80]);
%%
mainP='/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/meas_MID01106_FID09970_gSpi2d_T12_d110_Dw11_40rep_VD123_Thin';
load([mainP filesep 'PerSliceRec2.mat'],'Rec_CompgB0_MS','CompsPS','THLRMultiShotS');
nSlices=10;
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);
THLRMultiShotS=THLRMultiShotS(ROrd);
nTS=size(THLRMultiShotS{1},7);
TrgSz=gsize(THLRMultiShotS{1},1:2);
[~,~,~,H_AllTS]=ghankel(nTS,2,TrgSz);

AcqDwellTime_us=1.1;
nTrajPartMed=45000;
for s=1:nSlices
    disp(s);
    tmp=squeeze(THLRMultiShotS{s});
    [ ~, s_MS, V_MS] = batch_svd(H_AllTS*tmp);
    R1ms=V_MS(:,:,2,1)./V_MS(:,:,1,1); % R1 is simply the decay
    InnerTSDiff_PerShot_ms=nTrajPartMed*AcqDwellTime_us/1e3/(nTS-1);
    UpdatedT2SMap_ms=-InnerTSDiff_PerShot_ms./log(abs(R1ms));
    UpdatedB0Map_ms=-(angle(R1ms)/(2*pi))/(InnerTSDiff_PerShot_ms/1e3); % in Hz

    UpdatedB0Map_msS(:,:,s)=UpdatedB0Map_ms;
    UpdatedT2SMap_msS(:,:,s)=UpdatedT2SMap_ms;
    s_MSS(:,:,:,s)=s_MS;
    THLRMultiShotSM(:,:,:,s)=squeeze(tmp);
end

for s=1:nSlices
    tmp=grmss(THLRMultiShotSM(:,:,:,s),3);
%     tmp(tmp==0)=min(tmp(:)>0);
    [~,BMS,BNMS]=CalcSlicesSNR(tmp,false,7);
    BMSX(:,:,s)=~BNMS & tmp>0;
%     BMSX(:,:,s)=imfillholesBySlices(getLargestComponent(BMS));
end

THLRMultiShotSMN=THLRMultiShotSM./grms(THLRMultiShotSM,[1 2 3]);

fgmontage(squeeze(THLRMultiShotSMN(:,:,1,1:9)).*BMSX(:,:,1:9));title('"PD", TH-LR Multi-shot');

fgmontage(UpdatedB0Map_msS(:,:,1:9).*BMSX(:,:,1:9),[-100 100]);title('B_0, TH-LR Multi-shot');

fgmontage(UpdatedT2SMap_msS(:,:,1:9).*BMSX(:,:,1:9),[00 100]);title('T_2^*, TH-LR Multi-shot');
%%
QQ=squeeze(THLRMultiShotSMN(:,:,1,1:9)).*BMSX(:,:,1:9);
% QQ=UpdatedT2SMap_msS(:,:,1:9).*BMSX(:,:,1:9);
% QQ=UpdatedB0Map_msS(:,:,1:9).*BMSX(:,:,1:9);
A=PartitionDim(QQ,3,3);
B=CombineDims(A,[3 2]);
C=CombineDims(B,[3 1]);

fgmontage(C);removeTicks;axis equal;title('"PD", TH-LR Multi-shot');
% fgmontage(C,[0 100]);removeTicks;axis equal;title('T_2^*, TH-LR Multi-shot');
%%
s=1;
XX=squeeze(sum(Rec_CompgB0_MS{s}.*CompsPS{s},6)).*BMSX(:,:,s);
YY=cat(4,grmss(XX(:,:,:,3:7),4),grmss(XX(:,:,:,8:13),4),grmss(XX(:,:,:,14:20),4));
% ZZ=grmss(PartitionDim(XX(:,:,:,3:end),4,6),4);
% fgmontage(ZZ(:,:,3:7:end,:));
% ZZ=grmss(PartitionDim(XX,4,5),4);
% fgmontage(ZZ(:,:,3:7:end,2:end));
ZZ=grmss(PartitionDim(XX,4,4),4);
fgmontage(ZZ(:,:,3:11:end,2:end));
%%
Rec_CompgB0_MS=Rec_CompgB0_MS(ROrd);
CompsPS=CompsPS(ROrd);
%%
s=9;
XX=squeeze(sum(Rec_CompgB0_MS{s}.*CompsPS{s},6)).*BMSX(:,:,s);

Rep=5;

WhichEchos=4:20;
[~,~,~,H_SomeTS]=ghankel(numel(WhichEchos),2,TrgSz);

% ZZ=grmss(PartitionDim(XX,4,5),4);
% tmp=squeeze(ZZ(:,:,Rep,WhichEchos));
% InnerTSDiff_PerShot_ss=InnerTSDiff_PerShot_ms*4;
InnerTSDiff_PerShot_ss=InnerTSDiff_PerShot_ms*1;
tmp=squeeze(XX(:,:,Rep,WhichEchos));
[ ~, s_SS, V_SS] = batch_svd(H_SomeTS*tmp);
R1ss=V_SS(:,:,2,1)./V_SS(:,:,1,1); % R1 is simply the decay
UpdatedT2SMap_ss=-InnerTSDiff_PerShot_ss./log(abs(R1ss));
UpdatedB0Map_ss=-(angle(R1ss)/(2*pi))/(InnerTSDiff_PerShot_ss/1e3); % in Hz
fgmontage(UpdatedT2SMap_ss,[0 100])

UpdatedT2SMap_ssM(:,:,Rep)=UpdatedT2SMap_ss;

UpdatedT2SMap_ssM=max(0,UpdatedT2SMap_ssM);
UpdatedT2SMap_ssM=min(100,UpdatedT2SMap_ssM);
