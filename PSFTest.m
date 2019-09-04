g=CurGrad;
FOV_mm=200;
gammaMHz=42.574; % MHz/T
TwoPiGammaMHz=gammaMHz*2*pi;
GradDwellTime_ms=10e-3;

k=cumsum([0;g;0])*GradDwellTime_ms*TwoPiGammaMHz; % mT/m*ms * 2*pi*MHz/T = rad/m
if(ShowFigs)
    gPlotTraj_radm(g,FOV_mm);MaximizeFig
end

kK=k*FOV_mm/1000/2/pi;
%%
SpBW=1e6;
Acq_dT_us=1e6/SpBW;
AcqPointsPerGrad=GradDwellTime_ms*1e9/SpBW;

nAcqPoints=numel(CurGrad)*AcqPointsPerGrad;

TrajQ=interp1(0:numel(kK)-1,kK,((0:nAcqPoints-1)+AcqPointsPerGrad)/AcqPointsPerGrad).';
sum(isnan(TrajQ(:)))
BARTTraj=[real(TrajQ) imag(TrajQ) imag(TrajQ)*0].';

%% Apply traj on data
Nres=192;
HNres=floor(Nres/2);
Delta=zeros(Nres);
Delta(HNres,HNres)=1;
SensP=permute43(SensS);
IWithSens=Delta.*SensP;
Sig=bart('nufft ',BARTTraj, IWithSens);
%% Recon, no B0
Rec1c=bart('pics -t ',BARTTraj,Sig,SensP);
%% Recon, with added B0
TimePointsus=(1:nAcqPoints)*Acq_dT_us;
B0Hz=30;
PhaseAdd=exp(1i*2*pi*B0Hz*TimePointsus/1e6);
ModSig=Sig.*PhaseAdd;
Rec1cx=bart('pics -t ',BARTTraj,ModSig,SensP);
%% Now on sim data, no B0
SliIdx=26;
EchoIdx=10;
RefI=Full50EchoS(:,:,SliIdx,EchoIdx);
IWithSens=RefI.*SensSFull(:,:,SliIdx,:);
Sig=bart('nufft ',BARTTraj, IWithSens);
%% Recon, no B0
TrgSz=[96 96];
SensPSz=imresizeBySlices(SensSFull(:,:,SliIdx,:),TrgSz);
RefISz=imresizeBySlices(RefI,TrgSz);
RecIc=bart('pics -t ',BARTTraj,Sig,SensPSz);
RecIcN=abs(RecIc).*grmss(RefI)./grmss(RecIc);
ssimX=ssim(RecIcN,abs(RefISz));
figure;imshowpair(RecIcN,abs(RefISz),'montage');title(ssimX);



