BaseP='/home/g/Documents/UM_Stuff/MLN_NN_MB/';
M0='031_GRE_ASPIRE_2-1_test.nii';
P0='033_GRE_ASPIRE_2-1_test_combined.nii';

M0I=loadniidata([BaseP M0]);
P0I=loadniidata([BaseP P0]);
C0=M0I.*exp(1i*2*pi*P0I/4095);
clear M0I P0I
C0=single(C0);
I=permute(C0(:,:,:,1:2),[1 2 5 4 3]);
clear C0
nSlices=size(I,5);

% sTwix = mapVBVD(['/media/g/301A1F881A1F49EC/gUM/meas_MID38_gBP_ep2d_bold_multiecho_Noise_FID33419.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
% scanFreq_Hz=sTwix.hdr.Config.ScanFrequency; % SystemFrequency
scanFreq_Hz=297205524;
TEs_us=[4600 9340 14080];

%%
% I=I(:,:,:,:,50:1:end-50);
% nSlices=size(I,5);
%%
WhichTwo=[1 2];
% WhichTwo=[2 3];

% M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));

M=permute(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:),[1 2 3 5 4]);

% Mag=abs(M);
Mag=permute(mean(abs(I(:,:,:,WhichTwo,:)),4),[1 2 3 5 4]);
%%
Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));

% Combined=SmoothBySlices(Combined,[7 7],2.5);
% Combined=permute(sum(Mag.*exp(1i*angle(M)),3));
% Combined=sum(Combined,5);
gammaMHz=42.5774806;
gammaHz=gammaMHz*1e6;

deltaTE_us=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us/1e6)));
B0_Hz=dAngle/(2*pi*deltaTE_us/1e6);

% FirstEcho=squeeze(I(:,:,:,1,:));
FirstEcho=permute(I(:,:,:,1,:),[1 2 4 3 5]);

Mg=grmss(FirstEcho,3);
disp('ok');
%%
% fgmontage(Mg);
%%
MgAll=permute43(grmss(I,3));
RAll=MgAll(:,:,:,2:end)./MgAll(:,:,:,1:end-1);
dTEs_ms=diff(TEs_us,1,2)/1000;
LRAll=-log(RAll)./permute(dTEs_ms(1,1:size(RAll,4)),[1 4 3 2]);
T2SAll=1./LRAll;
% 
[Out B1 BN1]=CalcSlicesSNR(MgAll(:,:,:,1),false,5);
B2=imfillholesBySlices(~BN1);
T2SAll=T2SAll(:,:,:,1).*B2;
T2SAll(T2SAll<0 | T2SAll>300) =300;
T2SAll(T2SAll==0)=20;
T2SAll(~isfinite(T2SAll))=20;
% 
disp('ok t2*');
%%
% fgmontage(T2SAll(:,:,1:7:end),[0 100]);
% fgmontage(permute(T2SAll(1:9:end,:,:),[2 3 1]),[0 100]);
% fgmontage(permute(Mg(1:9:end,:,:),[2 3 1]));
%%
WhichSli=1:nSlices;
% WhichSli=1:38;
% WhichSli=15:20;

% [unwrapped] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);

unwrapped = single(robustunwrap(floor(size(dAngle)/2), double(dAngle), double(Mg)));
% unwrappedB = single(robustunwrap(floor(size(dAngle)/3), double(dAngle), double(Mg)));
% B=abs(unwrapped)<abs(unwrappedB);
% unwrapped=unwrapped.*B+unwrappedB.*(1-B);

B0_HzU=unwrapped/(2*pi*deltaTE_us/1e6);
disp('ok unwrapped to B0_HzU');
%% Restrict to range
ARange=[-500 500];
dHz=1e6/deltaTE_us;
B=B0_HzU<ARange(1);
d1=floor((B0_HzU-ARange(1))/dHz);
B0_HzU=B0_HzU-B.*(d1*dHz);

B=B0_HzU>ARange(2);
d2=ceil((B0_HzU-ARange(2))/dHz);
B0_HzU=B0_HzU-B.*(d2*dHz);
%%
fgmontage(B0_HzU(:,:,1:7:end),[-500 500]);colorbar
%%
% se = strel('disk',3,8);
% % SE = strel('sphere',7)
% %%
% S=B0_HzU(:,:,21);
% % SB2=imerode(B2(:,:,21),se)>-15+0.9;
% SM=Mg(:,:,21);
% Sz=size(S);
% dx=diff(S,1,1);
% dy=diff(S,1,2);
% dx(Sz(1),Sz(2))=0;
% dy(Sz(1),Sz(2))=0;
% ddS=sqrt(dx.*dx+dy.*dy);
% 
% dx=diff(SM,1,1);
% dy=diff(SM,1,2);
% dx(Sz(1),Sz(2))=0;
% dy(Sz(1),Sz(2))=0;
% ddSM=sqrt(dx.*dx+dy.*dy);
% 
% SD=ddSM>20;
% DSD=imdilate(SD,se);
% ddA=dd>(dHz/2) & (~DSD);
% %%
% fgmontage(B0_HzU(:,:,1:7:end),[-500 500]);colorbar