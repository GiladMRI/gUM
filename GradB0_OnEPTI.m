%%
Ch_Dim=4;
TS_Dim=3;

Ch_Flag=2^(Ch_Dim-1);
TS_Flag=2^(TS_Dim-1);

ScriptFN=[BaseSP 'MPBDScriptCartMPx.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% Copy with sens/CS mask, FFT and sample mask
ITS_Cmnds={['fmac 0 0'],'fftc 3','fmac 1 0'};
% Open over PD and then like ITS
MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
WriteLinopToFile(ScriptFN,MP_Cmnds);
disp('Wrote LinopScript');
%% Sim
S100um=load('/autofs/cluster/kawin/Gilad/Slices100um_Acq.mat');
%%
NBig=3200;
SzBig=[NBig NBig];
NAcq=120;
SzAcq=[NAcq NAcq];

C=S100um.Ax(:,:,5);
C(1)=C(1)+1i*1e-5;
C=imresize(C,SzBig);

[X, Y]=ndgrid(linspace(-1,1,NBig),linspace(-1,1,NBig));
CenterP=[-0.4,0];
% 2 is 200mm. We want one sided range of 2cm, i.e. of 10%, so sigma=~0.2.
% That'll be 10Hz/mm.
B0Sigma=[0.15 0.3];
RFromCSqr=(((X-CenterP(1))/B0Sigma(1)).^2+((Y-CenterP(2))./B0Sigma(2)).^2);
B0Max=200;
B0_Hz=exp(-RFromCSqr)*B0Max;
T2S_ms=B0_Hz*0+50;

nEchos=70;
EchoSpacing_us=690;
% TAcq_ms=25; % 10Hz/mm : pi/2 per mm, pi at 2mm
TAcq_ms=EchoSpacing_us*(nEchos-1)/1e3;
TimeBetweenEchos_ms=TAcq_ms/(nEchos-1);
EchoTimes_ms=TimeBetweenEchos_ms*(0:nEchos-1);
EchoTimes_ms3=permute32(EchoTimes_ms);
disp('OK loaded values for sim');
%%
% BaseRes_mm
DB0=cat(3,symD(B0_Hz,1),symD(B0_Hz,2)); %./perm42(BaseRes_mm);
figure;subplot(1,2,1);gmontage(DB0(:,:,1));colorbar
subplot(1,2,2);gmontage(DB0(:,:,2));colorbar
%%
nChannels=7;

if(~exist('Ma7TSensCollection','var'))
    Ma7TSensCollection=load('/autofs/cluster/kawin/Gilad/CCSensMaps.mat');
end
Sens=rot90(perm43(imresizeBySlices(perm31(Ma7TSensCollection.SensCC(1:nChannels,:,:,15)),SzBig)));
SensMsk=grmss(Sens,3:4)>0.05;
disp('Loaded sens');
%% random poisson
TotalOSRatio=1.5;
AccPerEcho=nEchos/TotalOSRatio;
SigCenterSize=8;

sAccPerEcho=sqrt(AccPerEcho);

clear PoisMsk
for i=1:nEchos
    disp(i);
    PoisMsk(:,:,i)=squeeze(bart(['poisson -Y ' num2str(NAcq) ' -Z ' num2str(NAcq) ' -y ' num2str(sAccPerEcho) ' -z ' num2str(sAccPerEcho) ' -C ' num2str(SigCenterSize) ' -s ' num2str(rand*100000)]));
end
Msk=PoisMsk;
disp('ok Msk');
%% EPI mask
LineOrder=round(linspace(1,nEchos,NAcq)).';
clear EPIMsk
for i=1:nEchos
    EPIMsk(:,:,i)=repmat(LineOrder==i,[1 NAcq]);
end
Msk=EPIMsk;
%% random lines
LineOrder=randi(30,1,NAcq).';
clear RandLineMsk
for i=1:nEchos
    RandLineMsk(:,:,i)=repmat(LineOrder==i,[1 NAcq]);
end
Msk=RandLineMsk;
%%
ImWithEchos=C.*exp(1i*2*pi*B0_Hz.*EchoTimes_ms3/1000).*exp(-EchoTimes_ms3./T2S_ms);
ImWithEchosWSens=ImWithEchos.*Sens;
FullData=fft2cg(ImWithEchosWSens);
FullDataAcq=crop(FullData,[SzAcq nEchos nChannels]);
disp('Calculated data');
%%
MaskedData=FullDataAcq.*Msk;
disp('Calculated sig');
%%
ScriptFN=[BaseSP 'MPBDScriptCartMPx.txt'];
% SensCSMap is 0
% mask_sampleP is 1
% fmac with sens/CS mask, FFT and sample mask
ITS_Cmnds={['fmac 0 0'],'fftc 3','fmac 1 0'};
% Open over PD and then like ITS
MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
WriteLinopToFile(ScriptFN,MP_Cmnds);
disp('Wrote LinopScript');
%%
% OSTestScriptFN=[BaseSP 'OS.txt'];
% WriteLinopToFile(OSTestScriptFN,OS_Cmnds);
%
% ITest=RecC{NRec};
% ImSzTest=FillOnesTo16(size(ITest));
% OS=bart(['linopScript -d 5 ' OSTestScriptFN],ImSzTest,ITest);
%%
% OS_Cmnds={'fftc 3','dblszc 3','ifftc 3'};
% OSScriptFN=[BaseSP 'MPBDScriptCartMPxOS.txt'];
% OS_MP_Cmnds=[OS_Cmnds,{'fmac 2 0'}, ITS_Cmnds];
% WriteLinopToFile(OSScriptFN,OS_MP_Cmnds);
% disp('Wrote OS LinopScript');
%% Now rec
NRecs=[100 200 400 800];
% NRecs=[400];
for i=1:numel(NRecs)
    NRec=NRecs(i);
    % NRec=800;
    SzRec=[NRec NRec];
    NHalfDiff=(NRec-NAcq)/2;
    
    SensRec=imresizeBySlices(Sens,SzRec);
    SensMskRec=grmss(Sensr,3:4)>0.05;
    
    B0_HzRec=imresizeBySlices(B0_Hz,SzRec);
    T2S_msRec=imresizeBySlices(T2S_ms,SzRec);
    T2S_msRec=max(4,min(T2S_msRec,300));
    
    TSC_Rec=exp(1i*2*pi*B0_HzRec.*EchoTimes_ms3/1000).*exp(-EchoTimes_ms3./T2S_msRec);
%     TSC_Rec=TSC_Rec/(NRec/100);
    
    PaddedData=padarray(MaskedData,[NHalfDiff NHalfDiff],'both');
    MskRec=padarray(Msk,[NHalfDiff NHalfDiff],'both');
    %
    ImSz16Rec=FillOnesTo16(SzRec);
    % 25 sec for 100, 50 for 200, 140 for 400
    disp(datestr(now));
    RecC{NRec}=bart(['picsS -w 1 -R W:3:0:0.001 -m ' ScriptFN],ImSz16Rec,PaddedData,SensRec,MskRec,TSC_Rec);
end
%%
NRecs=[100 200 400 800];
% NRecs=[100 200 400];
% NRecs=800;
for i=1:numel(NRecs)
    NRec=NRecs(i);

    nDbls=log2(NRec/NAcq);
    OS_Cmnds=[{'fftc 3'},repmat({'dblszc 3'},[1 nDbls]),{'ifftc 3'}];
    OSScriptFN=[BaseSP 'MPBDScriptCartMPxOS.txt'];
    OS_MP_Cmnds=[OS_Cmnds,{'fmac 2 0'}, ITS_Cmnds];
    WriteLinopToFile(OSScriptFN,OS_MP_Cmnds);
    disp('Wrote OS LinopScript');
    
    ImSz16Acq=FillOnesTo16([NAcq NAcq]);
    
    SzRec=[NRec NRec];
    NHalfDiff=(NRec-NAcq)/2;
    
    SensRec=imresizeBySlices(Sens,SzRec);
    SensMskRec=grmss(Sensr,3:4)>0.05;
    
    B0_HzRec=imresizeBySlices(B0_Hz,SzRec);
    T2S_msRec=imresizeBySlices(T2S_ms,SzRec);
    T2S_msRec=max(4,min(T2S_msRec,300));
    
    TSC_Rec=exp(1i*2*pi*B0_HzRec.*EchoTimes_ms3/1000).*exp(-EchoTimes_ms3./T2S_msRec);
%     TSC_Rec=TSC_Rec/(NRec.^2/10000);
    
    PaddedData=padarray(MaskedData,[NHalfDiff NHalfDiff],'both');
    MskRec=padarray(Msk,[NHalfDiff NHalfDiff],'both');
    %
    % 3 sec for 100, 20 for 200, 60 for 400, 255 for 800
    disp(datestr(now));
    RecC_OS{NRec}=bart(['picsS -w 1 -R W:3:0:0.001 -m ' OSScriptFN],ImSz16Acq,PaddedData,SensRec,MskRec,TSC_Rec);
end
%%
RecC_OSM=cat(3,RecC_OS{100},RecC_OS{200},RecC_OS{400},RecC_OS{800});
RecC_OSM=RecC_OSM./grms(RecC_OSM,1:2);
fgmontagex(RecC_OSM);

RecCM=cat(3,RecC{100},RecC{200},RecC{400},RecC{800});
RecCM=RecCM./grms(RecCM,1:2);
fgmontagex(RecCM);
figure;
subplot(2,2,1);
gmontage(RecC{100});
subplot(2,2,2);
gmontage(RecC{200});
subplot(2,2,3);
gmontage(RecC{400});
subplot(2,2,4);
gmontage(RecC{800});
%%
% save('GradB0_OnSim','RecC_OS','RecC');
%%
load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/meas_MID04675_FID21637_gre_te4_9/sTwixX.mat'); % sTwixX
% sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}
% dPhaseFOV: 213.1250
%     dReadoutFOV: 220
%     sTwixX.hdr.MeasYaps.sKSpace.lBaseResolution 128
%     sTwixX.hdr.MeasYaps.sKSpace.dPhaseResolution 1
%     sTwixX.hdr.MeasYaps.sKSpace.lPhaseEncodingLines 124
ROFOV=sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
RORes=sTwixX.hdr.MeasYaps.sKSpace.lBaseResolution;
SSRes=sTwixX.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
BaseRes_mm=[ROFOV/RORes ROFOV/RORes SSRes];
QQ=load('/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/meas_MID04675_FID21637_gre_te4_9/B0T2S.mat');
%%
DB0_Hzmm=cat(4,symD(QQ.B0M_Hz,1),symD(QQ.B0M_Hz,2),symD(QQ.B0M_Hz,3))./perm42(BaseRes_mm);
MskA=squeeze(QQ.s_vals(:,:,1,:))>3e-4;

Dy=symD(QQ.B0M_Hz,2);
Dz=symD(QQ.B0M_Hz,3);



fgmontagex(D1.*MskA)
fgmontagex(D2.*MskA)