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

BaseP='/autofs/cluster/kawin/Gilad/Bay4Prisma5ms10ms/';
FN='meas_MID01026_FID07879_gre_te7_40';

BaseP='/autofs/cluster/kawin/Gilad/Bay4Kawin5ms10ms/';
FN='meas_MID01104_FID09968_gre_te7_40';

BaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms/';
FN='meas_MID03482_FID20451_gre_te4_9';

BaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
FN='meas_MID04675_FID21637_gre_te4_9';

BaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68xms/';
FN='meas_MID00654_FID30513_gre_4echo_32_22';
FN2='meas_MID00655_FID30514_gre_4echo_43_22';

BaseP='/autofs/cluster/kawin/Gilad/EPTI_and_spi68msx_on_CL/';
FN='meas_MID00856_FID32091_gre_4echo_32_22';
FN2='meas_MID00857_FID32092_gre_4echo_43_22';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX=sTwix{end};
Data=sTwixX.image();

sTwix2 = mapVBVD([BaseP FN2 '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX2=sTwix2{end};
Data2=sTwixX2.image();

mkdir([BaseP FN]);
system(['chmod +777 -R ' BaseP FN]);
disp([BaseP FN ' Created']);
%%
save([BaseP FN filesep 'sTwixX.mat'],'sTwixX');
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
D=gpermute(Data,[1 3 4 2 8 5 6 7]);
D2=gpermute(Data2,[1 3 4 2 8 5 6 7]);
D=cat(5,D,D2);
D=D(:,:,:,:,[1 5 2 6 3 7 4 8],:);
clear D2 Data2
% PD=padLeft(D,24,2);
D=D(:,:,:,:,:,ROrd);
% PD=RepDotMult(PD,permute( mod(1:nSlices,2)*2-1,[1 6 3 4 5 2]));
I=squeeze(ifft2cg(D));

TEs_us=[sTwixX.hdr.Phoenix.alTE{:}];
TEs_us2=[sTwixX2.hdr.Phoenix.alTE{:}];
nEchos=size(I,4)/2;
TEs_us=cat(2,TEs_us(1:nEchos),TEs_us2(1:nEchos));
TEs_us=TEs_us([1 5 2 6 3 7 4 8]);

nEchos=size(I,4);
TEs_us=TEs_us(1:nEchos);
TEs_ms=TEs_us/1000;
dTEs_ms=diff(TEs_ms,1,2);

FirstEcho=squeeze(I(:,:,:,1,:));

SnsSzB=gsize(I,1:2);

disp('ok');
%%
for SliI=1:nSlices
    disp([num2str(SliI) ' ' datestr(now)]); % 45 sec per slice!
    SensBMM{SliI}=RunESPIRiTForSensMapsMultiMap(FirstEcho(:,:,:,SliI),0,SnsSzB);
%     SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);1
end
SensB=permute(cat(5,SensBMM{:}),[1 2 3 5 4]);

SensMsk=grmss(SensB(:,:,:,:,1),3)>0.01;
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
%%
% SensB [X Y Ch Sli Maps]
% I [X Y Ch Echo Slices]
Combined=gpermute(sum(I.*gpermute(conj(SensB(:,:,:,:,1)),[5 4]),3),[5 3]);
%%
WhichEchosToUse=1:nEchos;
for i=1:nSlices
    disp(i);
    [PDBase(:,:,i), UpdatedB0Map_Hz(:,:,i), UpdatedT2SMap_ms(:,:,i), s_vals(:,:,:,i), Fitted0(:,:,:,i), PDBase0(:,:,i)]=FitToModel_MPBD1CSf(Combined(:,:,i,:),WhichEchosToUse,dTEs_ms(1),TEs_ms(1));
end
%%
fgmontagex(abs(UpdatedT2SMap_ms(:,:,3:8:end)),[0 100]);title('GRE-ME T_2^*');
fgmontagex(UpdatedB0Map_Hz(:,:,3:8:end),[-300 300]);title(['GRE-ME B_0, TEs: ' num2str(TEs_ms,'%.1f ')]);
ShowAbsAngle(PDBase0(:,:,3:8:end),1,'Size',[2 2])
fgmontagex(s_vals(:,:,:,3:8:end));title('GRE-ME SV maps');
%%
% PDBase0x=min(PDBase0,6*grmss(PDBase0));
PDBase0x=min(PDBase0,6*median(abs(PDBase0(:))));
[Out B1 BN1]=CalcSlicesSNR(abs(PDBase0x(:,:,:)),false,5);
% [Out B1 BN1]=CalcSlicesSNR(abs(PDBase0x(:,:,:)),false,2);
B2=(~BN1).*SensMsk;
B2D=imdilate(B2,strel('disk',3,8));
% B2D=imdilate(B2,strel('disk',5,8));
B3=imfillholesBySlices( B2D );
for i=1:nSlices
    B4(:,:,i)=getLargestComponent(B3(:,:,i));
end
% B4=B3;
B4=B4.*SensMsk;
%%
dAngle=UpdatedB0Map_Hz*2*pi*dTEs_ms(1)/1000;
[unwrapped] = cusackUnwrap(dAngle, grmss(Combined,4));
unwrapped=unwrapped.*B4;
B0M_Hz=unwrapped*1000/2/pi/dTEs_ms(1);
fgmontagex(B0M_Hz(:,:,1:16),[-300 300]);title(['GRE-ME B_0 unwrapped, TEs: ' num2str(TEs_ms,'%.1f ')]);
fgmontagex(UpdatedB0Map_Hz(:,:,1:16),[-300 300]);title(['GRE-ME B_0, TEs: ' num2str(TEs_ms,'%.1f ')]);
%%
save([BaseP FN filesep 'B0T2S.mat'],'SensMsk','B1','BN1','B2','B2D','B3','B4',...
    'B0M_Hz','UpdatedB0Map_Hz','UpdatedT2SMap_ms','s_vals','PDBase0','TEs_ms');
disp('Saved');