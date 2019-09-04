BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/SPIRAL_ASL/S01/RawData/';
FN='meas_MID23_BP_fieldmap_9echos_2mm_Full_FID33404';

BaseP='/media/deni/bigdrive/SPIRAL_ASL/May_7_2019_Scan_SKOPE_TWIX/DENI_INVIVO_TWIX/';
FN='meas_MID766_BP_fieldmap_9echos_2mm_Full_FID51380';

BaseP='/autofs/cluster/kawin/Gilad/Skope_7May19/CRAZY_TRAJECTORIES_TWIX/';
FN='meas_MID853_BP_fieldmap_9echos_2mm_Full_FID51467';
FN='meas_MID881_BP_fieldmap_9echos_2mm_Full_low_FID51495';

% BaseP='/autofs/cluster/kawin/Gilad/';
% FN='meas_MID764_BP_fieldmap_5echosX_FID51378';

BaseP='/autofs/cluster/kawin/Gilad/PhantomBay4/';
FN='meas_MID00526_FID03646_gre_te7_40';

BaseP='/autofs/cluster/kawin/Gilad/Bay3/';
FN='meas_MID209_MEMPRAGE_4e_p2_1mm_iso_FID65136';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX=sTwix;
Data=sTwixX.image();

mkdir([BaseP FN]);
disp('Read data');
%%
nSlices=numel(sTwixX.hdr.Phoenix.sSliceArray.asSlice);

asSlice=sTwixX.hdr.Phoenix.sSliceArray.asSlice;
if(iscell(asSlice(1)))
    asSlice=[sTwixX.hdr.Phoenix.sSliceArray.asSlice{:}];
end

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

SlabSz=[asSlice.dReadoutFOV asSlice.dPhaseFOV asSlice.dThickness].';

RotMat = transpose(Quat2RotMat(sTwixX.image.slicePos(4:7, 100)));
RotatedLocs=RotMat.'*SlbLoc;
ZLocs=RotatedLocs(3,:);
[~, Ord]=sort(ZLocs);
% Ord=[1:2:nSlices 2:2:nSlices];
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);
disp('ok');
%%
save([BaseP FN filesep 'Locs.mat'],'RotatedLocs');
disp(['Saved ' BaseP FN filesep 'Locs.mat']);
%%
RelativeDisplacement=SlbLoc./SlabSz;
D=gpermute(Data,[1 3 4 2 8 5 6 7]);

% Shift1=permute(exp(-1i*linspace(-pi,pi,size(D,1))*(size(D,1)-1)*RelativeDisplacement(1)),[2 3 1]);
% Shift2=permute(exp(-1i*linspace(-pi,pi,size(D,2))*(size(D,2)-1)*RelativeDisplacement(2)),[1 2 3]);
% Shift3=permute(exp(-1i*linspace(-pi,pi,size(D,3))*(size(D,3)-1)*RelativeDisplacement(3)),[3 1 2]);
% 
% D=D.*Shift1.*Shift2.*Shift3;
% PD=padLeft(D,24,2);
% D=D(:,:,:,:,:,ROrd);
% PD=RepDotMult(PD,permute( mod(1:nSlices,2)*2-1,[1 6 3 4 5 2]));
I=squeeze(ifft3cg(D));
I=circshift(I,-1,3);
I=circshift(I,-1,2);
% I=circshift(I,-1,1);
I=flip(I,3);
MgAll=grmss(I,4);
% fgmontage(MgAll(:,:,20,1));title('MgAll');MaximizeFig
%%
fgmontage(MgAll(:,:,:,1));title('MgAll')

TEs_us=[sTwixX.hdr.Phoenix.alTE{:}];
disp('ok');
%%
% AA=loadniidata('/autofs/cluster/kawin/Gilad/Bay3/Dicoms/M.nii');
% fgmontage(AA(:,:,:,1));title('AA')
%%
% save([BaseP FN filesep 'IandsTwix.mat'],'I','sTwix');
% disp(['Saved ' BaseP FN filesep 'IandsTwix.mat']);
%%
% PD=D(:,1:48,1,:,:,:);
% PD=PD.*permute(hanning(48),[2 1]);
% PD=PD.*hanning(96);
% PD=padBoth(PD,24,2);
% PD=PD(:,:,:,:,:,ROrd);
% I=squeeze(fft3cg(PD));

MgAll=grmss(I,4);
RAll=MgAll(:,:,:,2:end)./MgAll(:,:,:,1:end-1);
dTEs_ms=diff(TEs_us,1,2)/1000;
LRAll=-log(RAll)./permute(dTEs_ms(1,1:size(RAll,4)),[1 4 3 2]);
T2SAll=1./LRAll;
% 
[Out B1 BN1]=CalcSlicesSNR(MgAll(:,:,:,1),false,5);
B2=imfillholesBySlices(~BN1);
T2SAll=T2SAll(:,:,:,1).*B2;
T2SAll(T2SAll<0 | T2SAll>300) =0;
T2SAll(T2SAll==0)=20;
% 
disp('ok t2*');
%%
% RAll=MgAll(:,:,:,3:end)./MgAll(:,:,:,1:end-2);
% dTEs_ms=(TEs_us(3:end)-TEs_us(1:end-2))/1000;
% LRAll=-log(RAll)./permute(dTEs_ms(1,1:size(RAll,4)),[1 4 3 2]);
% T2SAll=1./LRAll;
% % 
% [Out B1 BN1]=CalcSlicesSNR(MgAll(:,:,:,1),false,5);
% B2=imfillholesBySlices(~BN1);
% T2SAll=T2SAll(:,:,:,1).*B2;
% T2SAll(T2SAll<0 | T2SAll>300) =0;
% T2SAll(T2SAll==0)=20;
% % 
% disp('ok t2*');

%%
fgmontage(T2SAll,[0 100]);
gprint(get(gcf,'Number'),[BaseP FN filesep 'T2S'],[]) 
close(gcf);

save([BaseP FN filesep 'T2S.mat'],'T2SAll');
disp(['Saved ' BaseP FN filesep 'T2S.mat']);

% Mn=min(1./LRAll(:,:,:,1:2),[],4);
% Mx=max(1./LRAll(:,:,:,1:2),[],4);
% 
% B=Mn;
% B(Mx-Mn>20)=Mx(Mx-Mn>20);
% B(~B2)=0;
%%
I=permute(I,[1 2 4 5 3]);
%%
WhichTwo=[1 2];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
% Mag=abs(M);
Mag=squeeze(mean(abs(I(:,:,:,WhichTwo,:)),4));

Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));
% Combined=sum(Combined,5);
gammaMHz=42.5774806;
gammaHz=gammaMHz*1e6;

deltaTE_us=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));
% scanFreq_Hz=sTwix{1}.hdr.Config.ScanFrequency; % SystemFrequency
scanFreq_Hz=sTwixX.hdr.Dicom.lFrequency;

dAngle=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us/1e6)));
B0_Hz=dAngle/(2*pi*deltaTE_us/1e6);

FirstEcho=squeeze(I(:,:,:,1,:));
Mg=grmss(FirstEcho,3);
disp('ok');
%%
fgmontage(Mg);
gprint(get(gcf,'Number'),[BaseP FN filesep 'Mg'],[]) 
close(gcf);

save([BaseP FN filesep 'FirstEcho.mat'],'FirstEcho');
disp(['Saved ' BaseP FN filesep 'FirstEcho.mat']);
% load([BaseP FN filesep 'FirstEcho.mat'])
%%
WhichSli=1:nSlices;
% WhichSli=1:38;
% WhichSli=15:20;
if size(dAngle,3)>1
    [unwrapped] = cusackUnwrap(dAngle(:,:,:), double(Mg(:,:,:))/3000);
else
    [unwrapped] = cusackUnwrap(repmat(dAngle(:,:,WhichSli),[1 1 2]), repmat(Mg(:,:,WhichSli),[1 1 2])/3000);
    unwrapped=unwrapped(:,:,1);
end
%%
% dAngleToUse=dAngle;
% dAngleToUse(Mg<5e-5)=0;
% [unwrapped] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
%% Run this to work on smooth interecho angle : seems more robust
% SmCombined=SmoothBySlices(Combined,[5 5],1.5);
% dAngleSm=double(angle(SmCombined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us/1e6)));
% 
% [unwrapped] = cusackUnwrap(dAngleSm(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
%% Apply Cusack in blocks
% Starts=floor(linspace(1,nSlices,6));
% Ends=[Starts(2:end-1)-1 nSlices];
% Starts=Starts(1:end-1);
% Starts=8;
% Ends=10;
% WhichSli=[19 10 9  ];
% for i=1:numel(Starts)
%     WhichSli=Starts(i):Ends(i);
%     [unwrapped(:,:,WhichSli)] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
% end
%
% fgmontage(unwrapped)
%%
% unwrapped=unwrappeda;
% unwrapped(:,:,9)=unwrapped9(:,:,9);
% unwrapped(:,:,8)=unwrapped8(:,:,8);
% unwrapped(:,:,10)=unwrapped10(:,:,10);
% unwrapped(:,:,11)=unwrapped11(:,:,11);
% fgmontage(unwrapped)
%%
B0_HzU=unwrapped/(2*pi*deltaTE_us/1e6);
disp('ok unwrapped to B0_HzU');
%%
fgmontage(B0_HzU,[-500 500]);colorbar

gprint(get(gcf,'Number'),[BaseP FN filesep 'B0_HzU'],[]) 
close(gcf);

save([BaseP FN filesep 'B0_HzU.mat'],'B0_HzU');
disp(['Saved ' BaseP FN filesep 'B0_HzU.mat']);
% load([BaseP FN filesep 'B0_HzU.mat'])
%%
% SnsSzB=[128 128];
% SnsSzB=[96 96];
SnsSzB=gsize(B0_HzU,1:2);
B0Q=imresizeBySlices(B0_HzU,SnsSzB);
Mgc=imresizeBySlices(Mg,SnsSzB);
% Mskc=Mgc>7e-5;

Mskc=imresizeBySlices(B2,SnsSzB)>0;
B0M2=-B0Q;
B0M2(~Mskc)=0;
% SymMskC=abs(B0M2-gflip(B0M2,1))>230;
% SymMskC(1:30,:,:)=true;
% B0M2(SymMskC & B0M2>150)=-20;
% MskcE=imdilate(Mskc,strel('disk',5,8));
% B0M2(~MskcE)=0;
%%
fgmontage(B0M2,[-500 500]);colorbar

gprint(get(gcf,'Number'),[BaseP FN filesep 'B0M2'],[]) 
close(gcf);

save([BaseP FN filesep 'B0M2.mat'],'B0M2');
disp(['Saved ' BaseP FN filesep 'B0M2.mat']);
%%
nSlices=size(B0M2,3);
B0MForSmoothing=B0M2;
B0MForSmoothing(~Mskc)=0;
B0S=zeros([SnsSzB nSlices]);
for s=1:nSlices % 5 sec per slice
    disp(['Sli #' num2str(s) ' ' datestr(now)]);
    B0S(:,:,s)=FitB0ToRemoveGaps(B0MForSmoothing(:,:,s),Mgc(:,:,s),5);
end
disp('Finished B0S');
%%
fgmontage(B0S,[-500 500]);colorbar

gprint(get(gcf,'Number'),[BaseP FN filesep 'B0S'],[]) 
close(gcf);

save([BaseP FN filesep 'B0S.mat'],'B0S');
disp(['Saved ' BaseP FN filesep 'B0S.mat']);
%%
% SnsSzB=[128 128];
% for SliI=1:nSlices
%     disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
%     SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);
% end
%%
for SliI=1:nSlices
    disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
    SensBMM{SliI}=RunESPIRiTForSensMapsMultiMap(FirstEcho(:,:,:,SliI),0,SnsSzB);
%     SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);1

end
SensB=permute(cat(5,SensBMM{:}),[1 2 3 5 4]);
%%
save([BaseP FN filesep 'Sens.mat'],'SensB');
disp(['Saved ' BaseP FN filesep 'Sens.mat']);
%%
for SliI=1:nSlices
    figure;
    subplot(2,2,1);gmontage(abs(SensB(:,:,:,SliI,1)),[0 0.7]);
    subplot(2,2,2);gmontage(angle(SensB(:,:,:,SliI,1)),[-pi pi]);
    subplot(2,2,3);gmontage(abs(SensB(:,:,:,SliI,2)),[0 0.7]);
    subplot(2,2,4);gmontage(angle(SensB(:,:,:,SliI,2)),[-pi pi]);
%     ShowAbsAngle(SensB(:,:,:,SliI,:))
    YLbl=['Sli' num2str(SliI,'%02d')];
    ylabel(YLbl);
    gprint(get(gcf,'Number'),[BaseP FN filesep 'Sens_' YLbl],[]) 
    close(gcf);
end
disp('printed sens images');