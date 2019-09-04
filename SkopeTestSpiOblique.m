%% Data location: with matching TWIX file
dataFolder = '/autofs/cluster/kawin/Gilad/SkopeTest-10-04-19/';
scanId = 10;
DatFN='meas_MID355_BP_fieldmap_3echoes_tighly_spaced_FID46885.dat';
%% read data
scan = AqSysData(dataFolder,scanId);
%% SKOPE clean the data and show

%read raw data of 1st dynamic
raw = scan.getData('raw', [], [], [], 1);
%read phase data of 1st dynamic
% ophase = scan.getData('phase', [], [], [], 1);

%read (higher order) k-space coefficients of 2nd dynamic
%these coefficients relate to the basis set in kbase_spha.m and kbase_coco.m
kcoco =         scan.getData('kcoco', [], [], [], 1);
%these coefficients relate to kbase.m (spherical harmonics)
kspha =         scan.getData('kspha', [], [], [], 1);

%smooth data - use BW of gradient system
kFiltered = filterPhaseData(kspha(:,:,1,1), 1e-6, 1e-5);

%plot (linear) 2d k-space
figure, plot(kFiltered(:,2),kFiltered(:,3)); xlabel('kx [rad/m]'), ylabel('ky [rad/m]'), axis equal
figure, plot(kFiltered(:,1:4)); xlabel('Samples'), ylabel('k_0 [rad], k_x_y_z [rad/m]'), legend('k_0','k_x','k_y','k_z')
%% Take just kx,ky,kz
Data3=kFiltered(:,2:4);
%% Taking the rotation matrix from the TWIX data
twix_obj = mapVBVD([dataFolder DatFN],'removeOS');
% twix_obj.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal
RotMat = transpose(Quat2RotMat(twix_obj.image.slicePos(4:7, 100)));

ThetaOrig=-acos(RotMat(2,1));
%% Showing the trajectory
RData3=(RotMat'*(Data3.')).';
RFPoints=(1:1024*10)*2.5;
StartPoint=200;
clear IData
for i=1:3
    IData(:,i)=interp1(1:size(RData3,1),RData3(:,i),StartPoint+(RFPoints));
end

figure;
% subplot(1,2,1);
plot(IData(:,1:2));
hold on;
plot((IData(:,3)-mean(IData(:,3)))*20,'k','LineWidth',2);
title('Original rotation matrix');
removeTicks;
setXaxis([1,size(IData,1)]);
xlabel(num2str([ThetaOrig 0 pi/2]*180/pi,'%.2f'));
%% Find the peaks of the CAIPI blips
CAIPI_Period=200;
[~,MxCI]=max(RData3(:,3));

Xs=(CAIPI_Period*10+mod(MxCI,CAIPI_Period*2)):(CAIPI_Period*2):20000;

figure;plot(Xs(1):Xs(end),RData3(Xs(1):Xs(end),3)); hold on;
plot(Xs,RData3(Xs,3),'*');
%% Find the rotation matrix that gives the best CAIPI blips (i.e. the blips not going crazy around)
RotMatf=@(x) grotx(x(1))*groty(x(2))*grotz(x(3));

RotMatf2y=@(x) grotx(x(1))*groty(x(2))*grotz(pi/2);
RotMatf2z=@(x) grotx(x(1))*groty(0)*grotz(x(2));

x0=[ThetaOrig,0,pi/2];
x0y=[ThetaOrig,0];
x0z=[ThetaOrig,pi/2];
getValsFromRotated=@(X) X(Xs,3) ;
CostFunc=@(x) std(getValsFromRotated((RotMatf(x)'*(Data3.')).'));

CostFuncy=@(x) std(getValsFromRotated((RotMatf2y(x)'*(Data3.')).'));
CostFuncz=@(x) std(getValsFromRotated((RotMatf2z(x)'*(Data3.')).'));

[BestX,BestF]=fminsearch(CostFunc,x0);

[BestXy,BestFy]=fminsearch(CostFuncy,x0y);
[BestXz,BestFz]=fminsearch(CostFuncz,x0z);

% R3=RotMatf(BestX);

R3=RotMatf([BestXy pi/2]);
disp('Optimized rotation');
%% Showing the corrected trajectory
RData3=(RotMatf([BestXy pi/2])'*(Data3.')).';
% our scan is 25ms?
RFPoints=(1:1024*10)*2.5;
StartPoint=200;
clear IDatax
for i=1:3
    IDatax(:,i)=interp1(1:size(RData3,1),RData3(:,i),StartPoint+(RFPoints));
end

figure;
% subplot(1,2,2);
plot(IDatax(:,1:2));
hold on;
plot((IDatax(:,3)-mean(IDatax(:,3)))*20,'k','LineWidth',2);
removeTicks;
title('Optimized rotation matrix');
setXaxis([1,size(IDatax,1)]);
xlabel(num2str([BestXy pi/2]*180/pi,'%.2f'));