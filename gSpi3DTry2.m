ScanP='/media/a/DATA/SK_17Aug18/';
BaseFN='meas_MID92_gBP_ep3d_bold_multiecho_ASL_FID27295';
RefFldMapP='/media/a/DATA/SK_17Aug18/meas_MID35_BP_fieldmap_5echosX_FID27238/';

ScanP='/media/a/DATA/GL_Spi3D/';
BaseFN='meas_MID108_gBP_ep3d_bold_multiecho_ASL_2s_FID27699';
RefFldMapP=[ScanP 'meas_MID104_BP_fieldmap_5echosX_FID27695' filesep];

ScanP='/media/a/DATA/27Aug18_TO/';
BaseFN='meas_MID466_gBP_ep3d_bold_multiecho_ASL_2s_U05_FID28404';
RefFldMapP=[ScanP 'meas_MID467_BP_fieldmap_5echosX_FID28405' filesep]; % sNormal: dTra: 1

ScanP='/media/a/DATA/CW_03Sep18_3D2shotkCAIPI/';
% BaseFN='meas_MID426_gBP_ep3d_bold_multiecho_ASLb_FID29204'; % This is 1 interleave, sNormal: dTra: 1
BaseFN='meas_MID430_gBP_ep3d_bold_multiecho_ASLb_FID29208'; % 2 shots, sNormal: dTra: 1
RefFldMapP=[ScanP 'meas_MID424_BP_fieldmap_5echosX_FID29202' filesep]; % sNormal: dTra: 1

FN=[ScanP BaseFN '.dat'];
%%
mkdir([ScanP BaseFN]);
%% Read raw
AData = mapVBVD(FN);

% mrprot = rdMeas(FN) ;
% [mrprot, mdh, fid] = rdMeas(FN);
[mrprot, mdh, fid] = rdMeas(FN,'meas.asc',false);
% 32 channels, 6 reps, 5 ADCs = 960. *44 partitions = 42240
% size(fid): 1024       46336 : 4096 more, 128 per channel. It seems that
% the scan stqarts with 128 ADCs. Maybe for noise
% ulScanCounter: 1448
% QQ=AData.image.unsorted();

% mrprot.sWiPMemBlock.alFree.'
% mrprot.sWiPMemBlock.adFree(6:14).' % Interleaves seems like 11
disp('read data');
%%
Flds=fieldnames(mdh);
for i=1:numel(Flds)
    NInF(i)=numel(unique(mdh.(Flds{i})));
    [U IA IB]=unique(mdh.(Flds{i}));
    N{i}=hist(IB,0:numel(U));
end
[gmat2cell(num2str((1:numel(Flds))'),1) Flds gmat2cell(num2str(NInF'),1)]
%%
[U IA IB]=unique(mdh.ushSamplesInScan);
ChID=mdh.ulChannelId;
RepID=mdh.ushRepetition;
SegID=mdh.ushSeg;
nReps=max(RepID)+1;
% fidX=fid(:,IB==2); % to remove the noise 4096 = 128 blocks *32 channels
for r=1:(nReps-1) % Forget rep 0 that contains also noise scan
    for seg=0:max(SegID)
        for ch=0:max(ChID)
            CurB= ChID==ch & RepID==r & SegID==seg;
            Data(:,seg+1,ch+1,r,:)=squeeze(fid(:,CurB));
        end
    end
end
% size(Data) :  1024           5 ADCs          32 ch      5 rep 44 Kz
% figure;plot(grmss(Data,1:4)) % Kz: seems ok
DataC=CombineDims(Data,[2 1]);
% figure;plot(grmss(DataC,2:4)); % the ADCs seem to be one after the other
disp('ok');
%%
% DataCZ=fft1cg(DataC,4);
% %%
% CurData=DataCZ(:,:,1,22);
%%
FOVx=AData.hdr.Meas.ReadFOV;
dFOV=FOVx/1000;
%
paramLongROSamples = AData.hdr.MeasYaps.sWiPMemBlock.alFree{17};
spBW =AData.hdr.MeasYaps.sWiPMemBlock.adFree{15};
AccR =AData.hdr.MeasYaps.sWiPMemBlock.adFree{7};
paramLongInterleaves =AData.hdr.MeasYaps.sWiPMemBlock.adFree{11};
VD =AData.hdr.MeasYaps.sWiPMemBlock.adFree{6};
paramLongSpGradAmp =AData.hdr.MeasYaps.sWiPMemBlock.adFree{13};
paramLongSpSlewRate =AData.hdr.MeasYaps.sWiPMemBlock.adFree{12};

[kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,paramLongROSamples,spBW,AccR,...
        paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;
EffMaxRes_mm=FOVx/EffMaxRes;
% %%
% [kTraj, BaseRes, GradBuf, MaxGrad]=VDSpiralMex([dFOV,7*1024,spBW,2.7,...
%         paramLongInterleaves,VD,paramLongSpGradAmp,paramLongSpSlewRate,0]);
% EffMaxRes=sqrt(sum(((kTraj(end,:))*FOVx/2/pi/1000).^2))*2;
% EffMaxRes_mm=192/EffMaxRes
%%
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

% for i=1:numel(AData.hdr.Phoenix.sSliceArray.asSlice)
%     SLoc(i,1)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dSag;
%     SLoc(i,2)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dCor;
%     SLoc(i,3)=AData.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dTra;
% end

asSlice=AData.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[AData.hdr.Phoenix.sSliceArray.asSlice{:}];
end

nSlices=numel(asSlice);
for s=1:nSlices
    try
        SlbLoc(1,s)=asSlice(s).sPosition.dSag;
    catch
        SlbLoc(1,s)=0;
    end
    try
        SlbLoc(2,s)=asSlice(s).sPosition.dCor;
    catch
        SlbLoc(2,s)=0;
    end
    try
        SlbLoc(3,s)=asSlice(s).sPosition.dTra;
    catch
        SlbLoc(3,s)=0;
    end
end

RotMat = transpose(Quat2RotMat(AData.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;
disp('ok 1');
%%
MaxK=max(BARTTrajMS(:));
nTraj=size(BARTTrajMS,2);
Acc=ceil(MaxK*2).^2/nTraj;
figure;subplot(2,2,1);
plot(BARTTrajMS(1,:),BARTTrajMS(2,:),'.')
setXaxis([-1.1 1.1]*ceil(MaxK));
setYaxis([-1.1 1.1]*ceil(MaxK));
title(['MaxK=' num2str(MaxK) ' #Traj=' num2str(nTraj) ' Acc=' num2str(Acc)]);
subplot(2,2,2);
plot(GradBuf*MaxGrad*1000);title(['Grad, max=' num2str(MaxGrad*1000,'%.2f') 'mT/m'])
SlewBuf=diff(GradBuf*MaxGrad*1000,[],1);
subplot(2,2,4);
plot(SlewBuf*100);MaxSlew=max(max(abs(SlewBuf(20:end,:))));
title(['Slew, max~=' num2str(MaxSlew*100,'%.2f') 'mT/m/s'])
%%
FirstEcho=load([RefFldMapP 'FirstEcho.mat']);
FirstEcho=FirstEcho.FirstEcho;

FirstEcho=gflip(FirstEcho,1:2);
Mg=grmss(FirstEcho,3);

% Mg=circshift(Mg,-1,3);

SensB=load([RefFldMapP 'Sens.mat']);
SensB=SensB.SensB;
SnsSzB=gsize(SensB,1:2);

% SensB=circshift(SensB,-1,4);

SensX=SensB(:,:,:,6,1);
SensX=gflip(SensX,1:2);
disp('ok s');
%%
% SensXF=superkron(SensB(:,:,:,:,1),ones(1,1,1,4));
% SensXF=gflip(SensXF,1:2);
%%
% Txt=getLines('AddPhase_MID108.txt');
% Txt=getLines('AddPhase_466.txt');
Txt=getLines('AddPhase_430.txt');
Txt=Txt(1:2:end)';
clear AddPhase
for i=1:numel(Txt)
    AddPhase(i)=str2num(Txt{i});
end
AddPhase=AddPhase(3:end); % remove 2 ASL stuff in the beginning
nPars=size(DataC,4);
AddPhase=reshape(AddPhase,nPars,[]);
DataCA=DataC.*exp(-1i*permute(AddPhase(:,2:end)*pi/180,[3 4 2 1]));

try
    dx=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dSag/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
catch
    disp('No x shift!');
    dx=0;
end
% dy=-15;
% dy=AData.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition.dCor/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;

BARTTrajAct=BARTTrajMS(:,1:nTraj,:);

BARTTrajAct(1,:,:)=BARTTrajMS(1,:,:);
BARTTrajAct(2,:,:)=BARTTrajMS(2,:,:);

kx=BARTTrajAct(1,:,:)*2*pi;
ky=BARTTrajAct(2,:,:)*2*pi;

dy=RotatedLocs(1,1)/AData.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV;
% modx=double(exp(1i*(dx*kx+dy*ky))');
modx=permute(double(conj(exp(1i*(dx*kx+dy*ky)))),[2 1 3]);


DataCZ=ifft1cg(DataCA,4);
% CurData=DataCZ(:,:,5,8)
% CurData=DataCZ(:,:,5,32);

disp('ok c')
%%
s=22;
% for s=1:44
% SensX=SensB(:,:,:,s,1);
SensX=SensXF(:,:,:,s,1);
% SensX=SensB(:,:,:,5,1);
% SensX=SensB(:,:,:,17,1);
SensX=gflip(SensX,1:2);

% for r=1:5
r=9;
WhichI=2;
CurData=DataCZ(:,:,r,s);
% CurDataX=CurData(3:(nTraj+2),:).*modx;
CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));

%
% nukData=CurDataX(:,11).';
nukData=CurDataX.';
nukDataP=permute(nukData,[3 2 4 1]);

SensP=permute((SensX),[1 2 5 3 4]);

% SensP=ones(SnsSzB);

RecIfTVs=@(x) bart(['pics -S -m -R T:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP(:,:,:,:,1));
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP(:,:,:,:,1));


Lambda=1e-5;
% Rec=RecIfTVs(Lambda);
Rec=RecIfWs(Lambda);
% RecR(:,:,r)=Rec;
RecS(:,:,s)=Rec;
% end

fgmontage(RecS);
%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03')

nScc=5;

s=22;
r=9;
WhichI=2;
CurData=DataCZ(:,:,r,s);

CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));

%
% nukData=CurDataX(:,11).';
nukData=CurDataX.';
nukDataP=permute(nukData,[3 2 4 1]);

sccmtx = calcSCCMtx(permute(nukData,[3 2 1]));
SCCDATA = CC(permute(nukData,[3 2 1]),sccmtx(:,1:nScc));

% c=14;
for c=1:nScc
    disp(c);
%     nukDataP=permute(nukData,[3 2 4 1]);
nukDataP=SCCDATA(:,c).';

RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, ones(128));
Lambda=1e-5;
Rec=RecIfWs(Lambda);
RecC(:,:,c)=Rec;
end


%%
% s=22;
r=9;
WhichI=2;

for s=1:nPars
    CurData=DataCZ(:,:,r,s);
    
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
    
    nukData=CurDataX.';    
    sccmtx(:,:,s) = calcSCCMtx(permute(nukData,[3 2 1]));
    SCCDATA(:,:,s) = CC(permute(nukData,[3 2 1]),sccmtx(:,1:nScc,s));
end
disp('ok k');
% %%
% for s=1:nPars
%     for c=1:nScc
%         disp([s c]);
%         %     nukDataP=permute(nukData,[3 2 4 1]);
%         nukDataP=SCCDATA(:,c,s).';
%         
%         RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, ones(128));
%         Lambda=1e-5;
%         Rec=RecIfWs(Lambda);
%         RecSC(:,:,s,c)=Rec;
%     end
% end
% disp('ok sc');
%%
SensF=gflip(SensB(:,:,:,:,1),[1:2 4]);
%%
for s=1:nPars
    SensFCC(:,:,:,s) = CC(SensF(:,:,:,s),sccmtx(:,1:nScc,s));
end
disp('ok FCC');
%%
r=9;
WhichI=2;
for s=1:nPars
    disp(s);
    nukDataP=(permute(SCCDATA(:,:,s),[3 1 4 2]));
    
    SensP=permute(SensFCC(:,:,:,s),[1 2 4 3]);
    % SensP=gflip(SensP,1:2);
    
    RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP);
    Lambda=1e-5;
    Rec=RecIfWs(Lambda);
    RecSCC(:,:,s)=Rec;
end
disp('ok RecSCC');
%%
% r=9;
% WhichI=2;

for r=18:19
    WhichI=1+mod(r,2);
    for s=1:nPars
        disp([r s]);
        CurData=DataCZ(:,:,r,s);
        CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
        
        nukData=CurDataX.';
        SCCDATAR(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtx(:,1:nScc,s));
        
        disp(s);
        nukDataP=(permute(SCCDATAR(:,:,r,s),[3 1 4 2]));
        
        SensP=permute(SensFCC(:,:,:,s),[1 2 4 3]);
        % SensP=gflip(SensP,1:2);
        
        RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP);
        Lambda=1e-5;
        Rec=RecIfWs(Lambda);
        RecSCCR(:,:,r,s)=Rec;
    end
end
disp('ok RecSCCR');
%% Both interleaves
BARTTrajActBoth=[BARTTrajAct(:,:,1) BARTTrajAct(:,:,2)];
s=22;
rs=[18 19];
for s=1:nPars
    disp(s);
    CurDataA=DataCZ(:,:,rs(1),s);
    WhichIA=1+mod(rs(1),2);
    CurDataB=DataCZ(:,:,rs(2),s);
    WhichIB=1+mod(rs(2),2);
    CurDataXA=CurDataA(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIA));
    CurDataXB=CurDataB(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIB));
    
    nukData=[CurDataXA; CurDataXB].';
    SCCDATARBoth(:,:,s) = CC(permute(nukData,[3 2 1]),sccmtx(:,1:nScc,s));
    
    nukDataP=(permute(SCCDATARBoth(:,:,s),[3 1 4 2]));
    
    SensP=permute(SensFCC(:,:,:,s),[1 2 4 3]);
    
    RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajActBoth, nukDataP, SensP);
    Lambda=1e-5;
    Rec=RecIfWs(Lambda);
    RecSCCBoth(:,:,s)=Rec;
end
disp('ok RecSCCBoth');
%% As MB
DataCA2=DataCA(:,:,:,1:2:end);
DataCZ2=ifft1cg(DataCA2,4);

DataCA2b=DataCA(:,:,:,2:2:end);
DataCZ2b=ifft1cg(DataCA2b,4);
disp('ok m');
%%
SensFS=circshift(SensF,nPars/4,4);
% RecSCCS=circshift(RecSCC,nPars/4,3);
r=9;
WhichI=2;

half_nPars=nPars/2;

s=20;

for s=1:half_nPars
SliIs=[s s+half_nPars];
disp([r s]);
CurData=DataCZ2(:,:,r,s);
CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
        
nukData=CurDataX.';

sccmtxBoth(:,:,s) = calcSCCMtx(permute(nukData,[3 2 1]));
    
SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
        
SensFCCBoth(:,:,:,s) = CC(SensFS(:,:,:,s),sccmtxBoth(:,1:nScc,s));
SensFCCBoth(:,:,:,s+half_nPars) = CC(SensFS(:,:,:,s+half_nPars),sccmtxBoth(:,1:nScc,s));

SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);

nukDataP=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP);
Lambda=1e-5;
RecA=RecIfWs(Lambda);
RecAS(:,:,:,s)=squeeze(RecA);
end
%%
for s=1:half_nPars
    SliIs=[s s+half_nPars];
disp([r s]);

CurData=DataCZ2b(:,:,r,s);
CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
        
nukData=CurDataX.';

% sccmtxBoth(:,:,s) = calcSCCMtx(permute(nukData,[3 2 1]));
    
SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
        
% SensFCCBoth(:,:,:,s) = CC(SensFS(:,:,:,s),sccmtxBoth(:,1:nScc,s));
% SensFCCBoth(:,:,:,s+half_nPars) = CC(SensFS(:,:,:,s+half_nPars),sccmtxBoth(:,1:nScc,s));

SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);

nukDataP=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP);
Lambda=1e-5;
RecB=RecIfWs(Lambda);

RecBS(:,:,:,s)=squeeze(RecB);
end
%%
RecASX=circshift(CombineDims(RecAS,[3 4]),-nPars/4,3);
RecBSX=circshift(CombineDims(RecBS,[3 4]),-nPars/4,3);

W=abs(RecASX)+abs(RecBSX);
dAng=angle(RecASX)-angle(RecBSX);
WwithdAng=W.*exp(1i*dAng);

dAngExp=linspaceWithHalfStep(-pi,pi,nPars)-2*pi/(nPars*2);

dAngExp=linspace(-pi,pi,nPars+1);
dAngExp=dAngExp(1:end-1);

% dAngExp=linspacec(-pi,pi*(nPars-1)/nPars,nPars);
% dAngExp=linspace(-pi,pi,nPars);
figure;plot(angle(gsums(WwithdAng,1:2)));hold on;plot(dAngExp,'r');

dAngExpS=circshift(dAngExp,nPars/4);
%%
for s=1:half_nPars
SliIs=[s s+half_nPars];
disp([r s]);

CurData=DataCZ2b(:,:,r,s);
CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
        
nukData=CurDataX.';

% sccmtxBoth(:,:,s) = calcSCCMtx(permute(nukData,[3 2 1]));
    
SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
        
% SensFCCBoth(:,:,:,s) = CC(SensFS(:,:,:,s),sccmtxBoth(:,1:nScc,s));
% SensFCCBoth(:,:,:,s+half_nPars) = CC(SensFS(:,:,:,s+half_nPars),sccmtxBoth(:,1:nScc,s));

SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);

SensP=SensP.*exp(-1i*permute(dAngExpS(SliIs),[1 3 4 5 2]));

nukDataP=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
RecIfWs=@(x) bart(['pics -S -m -R W:7:0:' num2str(x) ' -t'],BARTTrajAct(:,:,WhichI), nukDataP, SensP);
Lambda=1e-5;
RecBb=RecIfWs(Lambda);

RecBbS(:,:,:,s)=squeeze(RecBb);
end
%%
RecBbSX=circshift(CombineDims(RecBbS,[3 4]),-nPars/4,3);

Wb=abs(RecASX)+abs(RecBbSX);
dAngb=angle(RecASX)-angle(RecBbSX);
WwithdAngb=Wb.*exp(1i*dAngb);

figure;plot(angle(gsums(WwithdAngb,1:2)));

%%
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03')

s=6;

for s=1:half_nPars
    SliIs=[s s+half_nPars];
    disp([r s]);
    CurData=DataCZ2(:,:,r,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
    
    nukData=CurDataX.';
    
    SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    
    SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);
    
    nukDataPA=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
    CurData=DataCZ2b(:,:,r,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
    
    nukData=CurDataX.';
    
    SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    
    nukDataPB=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
    
    setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TS')
    
    nukDataP=cat(2,nukDataPA,nukDataPB);
    
    nBands=2;
    nShots=2;
    
    nTS=1;
    TSBFA=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFB=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFAB=cat(2,TSBFA,TSBFB);
    TSBFAB=repmat(TSBFAB,[1 1 1 1 nBands]);
    
    TSBFAm=repmat(TSBFA,[1 1 1 1 nBands 1 1 1 1]);
    TSBFBm=TSBFAm.*exp(-1i*permute(dAngExpS(SliIs),[1 3 4 5 2]));
    TSBFm=cat(2,TSBFAm,TSBFBm);
    writecfl('/tmp/TSB',TSBFm);
    
    SensW=repmat(SensP(:,:,:,:,:),      [1 1 1 1 1 nTS]);
    
    TrajW=repmat(BARTTrajAct(:,:,WhichI),[1 nShots 1 nScc nBands nTS]);
    
    RecTS=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],TrajW, nukDataP, SensW);
    RecTSS(:,:,:,s)=RecTS;
end
RecTSX=circshift(CombineDims(RecTSS,[3 4]),-nPars/4,3);
%% Now same from 2 interleaves from whatever repetitions
for r=19
for s=1:half_nPars
    SliIs=[s s+half_nPars];
    disp([r s]);
    
    WhichI=1+mod(r,2);
    
    CurData=DataCZ2(:,:,r,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
    
    
    nukData=CurDataX.';
    
    SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    
    nukDataPA=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
    
    CurData=DataCZ2b(:,:,r,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichI));
    
    nukData=CurDataX.';
    
    SCCDATAR2(:,:,r,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    
    nukDataPB=(permute(SCCDATAR2(:,:,r,s),[3 1 4 2]));
        
    nukDataP=cat(2,nukDataPA,nukDataPB);
    
    nBands=2;
    nShots=2;
    
    nTS=1;
    TSBFA=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFB=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFAB=cat(2,TSBFA,TSBFB);
    TSBFAB=repmat(TSBFAB,[1 1 1 1 nBands]);
    
    TSBFAm=repmat(TSBFA,[1 1 1 1 nBands 1 1 1 1]);
    TSBFBm=TSBFAm.*exp(-1i*permute(dAngExpS(SliIs),[1 3 4 5 2]));
    TSBFm=cat(2,TSBFAm,TSBFBm);
    writecfl('/tmp/TSB',TSBFm);
    
    SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);
    SensW=repmat(SensP, [1 1 1 1 1 nTS]);
    TrajW=repmat(BARTTrajAct(:,:,WhichI),[1 nShots 1 nScc nBands nTS]);
    
    RecTS=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],TrajW, nukDataP, SensW);
    RecTSSR(:,:,r,:,s)=RecTS;
end
end
RecTSSRX=circshift(CombineDims(RecTSSR,[4 5]),-nPars/4,4);
% ShowAbsAngle(squeeze(RecTS))
%% Now mix both repetitions
s=4;
rA=19;
rB=19;
for s=1:half_nPars
SliIs=[s s+half_nPars];
    disp([r s]);
    
%     WhichIA=1+mod(rA,2);
%     WhichIB=1+mod(rB,2);
    
    WhichIA=1;
    WhichIB=3-WhichIA;
    
    CurData=DataCZ2(:,:,rA,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIA));
    nukData=CurDataX.';
    SCCDATAR2(:,:,rA,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    nukDataPA=(permute(SCCDATAR2(:,:,rA,s),[3 1 4 2]));
    
    CurData=DataCZ2b(:,:,rB,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIB));
    nukData=CurDataX.';
    SCCDATAR2(:,:,rB,s) = CC(permute(nukData,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    nukDataPB=(permute(SCCDATAR2(:,:,rB,s),[3 1 4 2]));
        
    nukDataP=cat(2,nukDataPA,nukDataPB);
    
    nBands=2;
    
    nTS=1;
    TSBFA=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFB=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFAB=cat(2,TSBFA,TSBFB);
    TSBFAB=repmat(TSBFAB,[1 1 1 1 nBands]);
    
    TSBFAm=repmat(TSBFA,[1 1 1 1 nBands 1 1 1 1]);
    TSBFBm=TSBFAm.*exp(-1i*permute(dAngExpS(SliIs),[1 3 4 5 2]));
    TSBFm=cat(2,TSBFAm,TSBFBm);
    writecfl('/tmp/TSB',TSBFm);
    
    SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);
    SensW=repmat(SensP, [1 1 1 1 1 nTS]);
    
    TrajW=cat(2,BARTTrajAct(:,:,WhichIA),BARTTrajAct(:,:,WhichIB));
    TrajW=repmat(TrajW,[1 1 1 nScc nBands nTS]);
    
    RecTS=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],TrajW, nukDataP, SensW);
    RecTSAB(:,:,:,s)=RecTS;
end
RecTSABX=circshift(CombineDims(RecTSAB,[3 4]),-nPars/4,3);

%     ShowAbsAngle(RecTS);
%% Now try B0
B0S=load([RefFldMapP 'B0S.mat']);
B0S=B0S.B0S;
B0S=gflip(B0S,1:3);
disp('ok');

Sz2=SnsSzB;
Mgc=gflip(imresizeBySlices(Mg,Sz2),1:3);
Mskc=Mgc>7e-5;
MskcE=imdilate(imfillholesBySlices(Mskc),strel('disk',5,8));
%%
TimeInMs2=(0:nTraj-1)*2.5/1e3;

[U_TimeInMs2, IA_TimeInMs2, IB_TimeInMs2]=unique(TimeInMs2);
nU_TimeInMs2=numel(U_TimeInMs2);

%% Fessler time segmentation
clear ErrTS
TS_Thresh=1e-5;
nTS=11;
TSCF=zeros([Sz2 nTS nPars]);
TSBF=zeros([nTraj nTS nPars]);

s=1;
for s=1:nPars
    disp(s);
    B0M2=B0S(:,:,s);
    
    AllB0C=exp(1i*2*pi*RepDotMult(B0M2,gpermute(TimeInMs2(IA_TimeInMs2)/1000,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    E=reshape(AllB0C,prod(Sz2),nU_TimeInMs2);
    MgcN=Mgc(:,:,s)./grmss(Mgc(:,:,s));
    WE=Col(MgcN);
    % WE=Col(MgcN)*0+1;
    
    WeightedE=WE.*E;
    
    FesTimePoints=linspace(0,TimeInMs2(end)/1000,nTS);
    TSC=exp(1i*2*pi*RepDotMult(B0M2,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);
    TSC2=reshape(TSC,prod(Sz2),nTS);
    WTSC2=WE.*TSC2;
    TSB=(WeightedE.')/(WTSC2.');% W in both sides
    
    ErrTS(s)=grmss(WeightedE-WTSC2*(TSB.')); %include W
    
    TSCF(:,:,:,s)=TSC;
    TSBF(:,:,s)=TSB;
    disp([datestr(now) ' nTS ' num2str(nTS) ' err=' num2str(ErrTS(s))]);
end
%%
save([ScanP BaseFN filesep 'TSData.mat'],'TSCF','TSBF');
%% Now with TS:
setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03TS')
% s=1;
for s=1:half_nPars
    disp(s);
    SliIsx=mod([s+nPars*3/4 s+nPars/4]-1,nPars)+1;
    rA=19;
    rB=19;

    WhichIA=1;
    WhichIB=3-WhichIA;
    
    SliIs=[s s+half_nPars];
        
    CurData=DataCZ2(:,:,rA,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIA));
    nukDataA=CurDataX.';
    
    CurData=DataCZ2b(:,:,rB,s);
    CurDataX=CurData(3:(nTraj+2),:).*squeeze(modx(:,:,WhichIB));
    nukDataB=CurDataX.';
    
    
    
    sccmtxBoth(:,:,s) = calcSCCMtx(permute(nukDataA,[3 2 1]));
    
    SensFCCBoth(:,:,:,s) = CC(SensFS(:,:,:,s),sccmtxBoth(:,1:nScc,s));
    SensFCCBoth(:,:,:,s+half_nPars) = CC(SensFS(:,:,:,s+half_nPars),sccmtxBoth(:,1:nScc,s));

    SCCDATAR2(:,:,rA,s) = CC(permute(nukDataA,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    nukDataPA=(permute(SCCDATAR2(:,:,rA,s),[3 1 4 2]));
    
    
    SCCDATAR2(:,:,rB,s) = CC(permute(nukDataB,[3 2 1]),sccmtxBoth(:,1:nScc,s));
    nukDataPB=(permute(SCCDATAR2(:,:,rB,s),[3 1 4 2]));
        
    nukDataP=cat(2,nukDataPA,nukDataPB);
    
    nBands=2;
    
%     nTS=11;
%     TSBFA=ones([1 nTraj 1 nScc 1 nTS])/nTS;
%     TSBFB=ones([1 nTraj 1 nScc 1 nTS])/nTS;
    TSBFA=permute(TSBF(:,:,SliIsx),[4 1 6 5 3 2]);
    TSBFB=permute(TSBF(:,:,SliIsx),[4 1 6 5 3 2]);
%     TSBFB=permute(TSBF(:,:,SliIsx(2)),[4 1 3 2])/nTS;
%     TSBFAB=cat(2,TSBFA,TSBFB);
%     TSBFAB=repmat(TSBFAB,[1 1 1 1 nBands]);
    
    TSBFAm=repmat(TSBFA,[1 1 1 nScc 1 1 1 1 1]);
    TSBFBm=TSBFAm.*exp(-1i*permute(dAngExpS(SliIs),[1 3 4 5 2]));
    TSBFm=cat(2,TSBFAm,TSBFBm);
    TSBFm=TSBFm/nTS;
    writecfl('/tmp/TSB',TSBFm);
    
    SensP=permute(SensFCCBoth(:,:,:,SliIs),[1 2 5 3 4]);
    SensW=SensP.*permute(TSCF(:,:,:,SliIsx),[1 2 5 6 4 3])*nTS;
    
    TrajW=cat(2,BARTTrajAct(:,:,WhichIA),BARTTrajAct(:,:,WhichIB));
    TrajW=repmat(TrajW,[1 1 1 nScc nBands nTS]);
    
    disp(datestr(now));
    RecTS=bart(['pics -S -m -R W:7:0:' num2str(1e-7) ' -t'],TrajW, nukDataP, SensW);
    disp(datestr(now));
    
    RecTSAB_B0SE(:,:,:,s)=RecTS;
end
RecTSAB_B0SEX=rot90(circshift(CombineDims(RecTSAB_B0SE,[3 4]),-nPars/4,3));
% save('kCAIPI_On_CW.mat')
% fgmontage(RecTS);
