% run tgv_recon_4ch_brain_radial.m before
load('After_tgv_recon_4ch_brain_radial.mat');
%%
figure;
plot(real(k(:,1:1:end)),imag(k(:,1:1:end)))


imgTGV2 = tgv2_l2_2D_pd(imgSens, sqrt(repmat(w,[1,1,nCh])).*rawdata, FT,repmat(w,[1,1,nCh]), 2*alpha, alpha, maxit, reduction, innerIter);
disp(['Time: ', num2str(toc/60), ' min']);
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03')

ESens=RunESPIRiTForSensMaps(imgSens);
SensP=permute(ESens,[1 2 4 3]);

Traj=CombineDims(k,[2 1])*size(imgSens,1);

Data=CombineDims(rawdata,[2 1]);
DataP=permute(Data,[3 1 4 2]);

clear BARTTraj
BARTTraj(1,:)=real(Traj).';
BARTTraj(2,:)=imag(Traj).';
BARTTraj(3,1)=0;


RecW=bart(['pics -S -m -R W:3:0:' num2str(0.0003) ' -t'],BARTTraj, DataP, SensP);
RecT=bart(['pics -S -m -R T:3:0:' num2str(0.0003) ' -t'],BARTTraj, DataP, SensP);
%%
TGVBase=permute([1;0;0],[8 7 6 5 4 3 2 1]);
writecfl('/tmp/TGVBase',TGVBase);

% Data=CombineDims(rawdata,[2 1]);

SensP=permute(imgSens,[1 2 4 3]);

% [u      du/dx   du/dy;
%  p1     dp1/dx  dp1/dy;
%  p2     dp2/dx  dp2/dy]
TGVTensor=complex(zeros(5,3,3));
TGVTensor(1,:,:)= [0 1 0; -1 0 0; 0 0 0];% du/dx-p1
TGVTensor(2,:,:)= [0 0 1; 0 0 0; -1 0 0];% du/dy-p2
TGVTensor(3,2,2)=1; % dp1/dx
TGVTensor(4,3,3)=1; % dp2/dy
TGVTensor(5,2,3)=1/2; % dp1/dy + dp2/dx
TGVTensor(5,3,2)=1/2;
TGVTensor(6,2,3)=1/2; % dp1/dy + dp2/dx
TGVTensor(6,3,2)=1/2;

TGVTensor1=TGVTensor(1:2,:,:);
TGVTensor2=TGVTensor(3:6,:,:);
TGVTensor1=permute(TGVTensor1,[9 8 7 6 5 4 10 2 1 15 14 13 12 11 3]);
TGVTensor2=permute(TGVTensor2,[9 8 7 6 5 4 10 2 1 15 14 13 12 11 3]);
writecfl('/tmp/TGVTensor1',TGVTensor1);
writecfl('/tmp/TGVTensor2',TGVTensor2);

SensP2=SensP;
SensP2(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

DataP2=DataP;
DataP2(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

writecfl('/tmp/RR',repmat(RecW,[1 1 1 1 1 1 1 3]));

setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TGV')

Rho=2;

tic
% RegW=0.0015;BaseAdd=.5;
RegW=0.003;BaseAdd=.2;Range=[0 0.018];
% RegW=0.004;BaseAdd=.1;Range=[0 0.015];
% RegW=0.0095;BaseAdd=0.05;
Data=CombineDims(  (abs(repmat(w+BaseAdd,[1,1,nCh])).^0.5   ).*rawdata,[2 1]);
DataP=permute(Data,[3 1 4 2]);

RecG2G=bart(['pics -S -m  -R G1:3:0:' num2str(RegW) ' -R G2:3:0:' num2str(RegW*Rho) ' -t'],BARTTraj, DataP2, SensP2);
t1=toc;
tic
RecG2T=bart(['pics -S -m -R T:3:0:' num2str(RegW) ' -t'],BARTTraj, DataP2, SensP2);
t2=toc;
[t1 t2]

RecG2GS=squeeze(RecG2G);
RecG2TS=squeeze(RecG2T);

RecT=RecG2TS(:,:,1);
RecG=RecG2GS(:,:,1);

Xs=30:140;Ys=30:140;
Xs=20:230;Ys=Xs;
% Range=[0 0.01];
figure;
ha = tight_subplot(1,3,[.0001 .0001],[.1 .1],[.01 .01]);
axes(ha(1));
gmontage(RecT(Ys,Xs),Range);title('BART TV');axis equal
axes(ha(2));
gmontage(RecG(Ys,Xs),Range);title('BART TGV');axis equal
xlabel(['lambda=' num2str(RegW) ' Rho=' num2str(Rho) ' BaseAdd ' num2str(BaseAdd)]);
axes(ha(3));
gmontage(imgTGV2(Ys,Xs),[0 0.7]);title('imgTGV2');axis equal
linkaxes(ha);

% fgmontage(RecG2GS(:,:,2:3))
%%
Dx=cat(1,zeros(1,256),diff(RecG2GS(:,:,1),1,1));
Dy=cat(2,zeros(256,1),diff(RecG2GS(:,:,1),1,2));
Du=cat(3,Dx,Dy);
p=RecG2GS(:,:,2:3);
Dump=Du-p;
All=cat(4,p,Du,Dump);
% fgmontage(All)
SAll=grmss(All,3);
fgmontage(SAll)

figure;
ha = tight_subplot(2,2,[0 0],[.1 .1],[.01 .01]);
axes(ha(1));
gmontage(SAll(:,:,1),[0 4e-3]);axis equal
axes(ha(2));
gmontage(RecG(:,:,1),Range);axis equal
axes(ha(3));
gmontage(SAll(:,:,2),[0 4e-3]);axis equal
axes(ha(4));
gmontage(SAll(:,:,3),[0 4e-3]);axis equal
linkaxes(ha);

%%
% I=double(imread('TGV.png'));
I=double(rgb2gray(imread('double_gradient.png')));
% I=double(imread('j_jiip-2016-0051_fig_0059.png'));

Acc=1;
% Msk=squeeze(bart(['poisson -Y 128 -Z 128 -y ' num2str(Acc) ' -z ' num2str(Acc) ' -v -e -C 8']));
Msk=I*0+1;
FI=fft2cg(I);

SensA=FI*0+1;
SensA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

DataA=FI.*Msk;
DataA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

MskA=Msk;
MskA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

RegWA=3;
RegWA2=RegWA*2;
RecTA=bart(['pics -S -m  -R T:3:0:' num2str(RegWA) ' -p'],MskA, DataA, SensA);
RecGA=bart(['pics -S -m  -R G1:3:0:' num2str(RegWA) ' -R G2:3:0:' num2str(RegWA2) ' -p'],MskA, DataA, SensA);

figure;
subplot(1,2,1);
gmontage(RecGA(:,:,1));title('BART TGV');axis equal
subplot(1,2,2);
gmontage(RecTA(:,:,1));title('BART TV');axis equal

% fgmontage(RecGA(:,:,2:3))
%%
% I=double(rgb2gray(imread('PenguinNoise.png')));
I=double(rgb2gray(imread('PenguinNoise.png')));
Acc=1;
% Msk=squeeze(bart(['poisson -Y 128 -Z 128 -y ' num2str(Acc) ' -z ' num2str(Acc) ' -v -e -C 8']));
Msk=I*0+1;
FI=fft2cg(I);

SensA=FI*0+1;
SensA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

DataA=FI.*Msk;
DataA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

MskA=Msk;
MskA(1,1,1,1,1,1,1,3,1,1,1,1,1,1)=0;

RegWA=.03;
RegWA2=RegWA*2;
RecTA=bart(['pics -S -m  -R T:3:0:' num2str(RegWA/1) ' -p'],MskA, DataA, SensA);
RecGA=bart(['pics -S -m  -R G1:3:0:' num2str(RegWA) ' -R G2:3:0:' num2str(RegWA2) ' -p'],MskA, DataA, SensA);

figure;
subplot(1,2,1);
gmontage(RecGA(:,:,1));title('BART TGV');axis equal
subplot(1,2,2);
gmontage(RecTA(:,:,1));title('BART TV');axis equal

% fgmontage(RecGA(:,:,2:3))
%%
fgmontage(RecT);
fgmontage(RecG);
%% Verify
fgmontage(RecG2TS(:,:,2:3))
fgmontage(RecG2GS(:,:,2:3))
%%
save('StatusForBART_TGV_Test.mat');
%%
DG=RecG(:,:,:,:,:,:,:,1)-RecG1(:,:,:,:,:,:,:,1);
grmss(DG)
grmss(RecG(:,:,:,:,:,:,:,1))
grmss(RecG1(:,:,:,:,:,:,:,1))