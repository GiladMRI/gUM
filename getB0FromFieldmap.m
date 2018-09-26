BaseP='/media/a/DATA/2018_01_25/';
FN='meas_MID131_dt_fieldmap_iso4mm_Ice_FID275';
FN='meas_MID131_dt_fieldmap_iso4mm_Ice_FID275';

BaseP='/media/a/DATA1/13May18/Phantom/';
FN='meas_MID377_BP_fieldmap_4echos_FID17766';

BaseP='/media/a/DATA1/13May18/Phantom/';
FN='meas_MID377_BP_fieldmap_4echos_FID17766';

BaseP='/media/a/DATA/13May18/Me/';
FN='meas_MID399_BP_fieldmap_4echos_FID17788';

BaseP='/media/a/DATA/14May18/Ben/';
FN='meas_MID123_BP_fieldmap_5echosX_FID17958';

BaseP='/media/a/DATA/180628_AK/';
FN='meas_MID265_BP_fieldmap_5echosX_FID22460';

BaseP='/media/a/DATA/11Jul18/RL/';
FN='meas_MID141_BP_fieldmap_5echosX_FID23838';

% BaseP='/media/a/DATA/PhantomCAIPI/';
% FN='meas_MID330_BP_fieldmap_5echosX_FID24027';

BaseP='/media/a/DATA/ASLSubjData/S01/';
FN='meas_MID399_BP_fieldmap_4echos_FID17788';

BaseP='/media/a/DATA/ASLSubjData/S02/';
FN='meas_MID123_BP_fieldmap_5echosX_FID17958';

% BaseP='/media/a/DATA/ASLSubjData/S03/';
% FN='meas_MID265_BP_fieldmap_5echosX_FID22460';
% 
% BaseP='/media/a/DATA/ASLSubjData/S04/';
% FN='meas_MID141_BP_fieldmap_5echosX_FID23838';

% BaseP='/media/a/DATA/FC/';
% FN='meas_MID197_BP_fieldmap_5echosX_FID24583';
% FN='meas_MID156_BP_fieldmap_5echosX_oblique_FID24549';

BaseP='/home/deni/DI_AK_ObliqueMBOutInOut/';
FN='meas_MID167_BP_fieldmap_5echosX_FID25537';

BaseP='/home/deni/GL/';
FN='meas_MID399_BP_fieldmap_4echos_FID17788';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnBP_9Aug18/';
FN='meas_MID484_BP_fieldmap_5echosX_FID26052';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnBP_9Aug18/DeeperRegion/';
FN='meas_MID507_BP_fieldmap_5echosX_FID26075';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/DI_AK_ObliqueMBOutInOut/';
FN='meas_MID167_BP_fieldmap_5echosX_FID25537';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnSK_17Aug18/';
FN='meas_MID35_BP_fieldmap_5echosX_FID27238';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/27Aug18_TO/';
FN='meas_MID446_BP_fieldmap_5echosX_FID28384';

BaseP='/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/Sep19_Deni/';
FN='meas_MID537_BP_fieldmap_5echosX_FID30388';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg','doAverage');
Data=sTwix.image();

mkdir([BaseP FN]);
%%
nSlices=numel(sTwix.hdr.Phoenix.sSliceArray.asSlice);
% Ord=[1:2:nSlices 2:2:nSlices];
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);
disp('ok');
%%
D=gpermute(Data,[1 3 4 2 8 5 6 7]);
PD=padLeft(D,24,2);
PD=PD(:,:,:,:,:,ROrd);
% PD=RepDotMult(PD,permute( mod(1:nSlices,2)*2-1,[1 6 3 4 5 2]));
I=squeeze(fft3cg(PD));

TEs_us=[sTwix.hdr.Phoenix.alTE{:}];
disp('ok');
%%
save([BaseP FN filesep 'IandsTwix.mat'],'I','sTwix');
disp(['Saved ' BaseP FN filesep 'IandsTwix.mat']);
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
scanFreq_Hz=sTwix.hdr.Config.ScanFrequency; % SystemFrequency

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
[unwrapped] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
%%
% dAngleToUse=dAngle;
% dAngleToUse(Mg<5e-5)=0;
% [unwrapped] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
%% Run this to work on smooth interecho angle : seems more robust
SmCombined=SmoothBySlices(Combined,[5 5],1.5);
dAngleSm=double(angle(SmCombined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us/1e6)));

[unwrapped] = cusackUnwrap(dAngleSm(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
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
fgmontage(unwrapped)
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
SnsSzB=[128 128];
B0Q=imresizeBySlices(B0_HzU,SnsSzB);
Mgc=imresizeBySlices(Mg,SnsSzB);
Mskc=Mgc>7e-5;
B0M2=-B0Q;
% B0M2(~Mskc)=0;
SymMskC=abs(B0M2-gflip(B0M2,1))>230;
SymMskC(1:30,:,:)=true;
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
%     SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);
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