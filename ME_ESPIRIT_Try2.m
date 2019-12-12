BaseP='/autofs/space/daisy_001/users/data/Gilad/gep_CL/';
FN='meas_MID01928_FID43869_gre_4echo_24_26_G2';
FN2='meas_MID01929_FID43870_gre_4echo_37_26_G2';

sTwix = mapVBVD([BaseP FN '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX=sTwix{end};

sTwix2 = mapVBVD([BaseP FN2 '.dat'],'removeOS','ignoreSeg','doAverage','rampSampRegrid');
sTwixX2=sTwix2{end};
%%
nSlices=numel(sTwixX.hdr.Phoenix.sSliceArray.asSlice);
Ord=[2:2:nSlices 1:2:nSlices];
[~,ROrd]=sort(Ord);

load([BaseP FN filesep 'resAll.mat'],'resAll');
resAll1=resAll;
load([BaseP FN2 filesep 'resAll.mat'],'resAll');
%%
D=cat(5,resAll1,resAll);
D=D(:,:,:,:,[1 5 2 6 3 7 4 8],:);
D=permute(D,[1 2 3 6 5 4]);
D=D(:,:,:,:,:,ROrd);
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

load([BaseP FN filesep 'Sens.mat'],'SensB');

Combined=gpermute(sum(I.*gpermute(conj(SensB(:,:,:,:,1)),[5 4]),3),[5 3]);
%%
WhichEchosToUse=1:nEchos;
[PDBase, UpdatedB0Map_Hz, UpdatedT2SMap_ms, s_vals, Fitted0, PDBase0]=...
    FitToModel_MPBD1CSf(perm43(Combined),WhichEchosToUse,dTEs_ms(1),TEs_ms(1));
%%
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;
%%
FCombined=fft2cg(Combined); % [x y Slices Echos]
Sz=gsize(FCombined,1:2);
%%
% save('ME_ESPIRIT_Try2.mat');
% load('ME_ESPIRIT_Try2.mat'); % reasy here
%%
% CalibRegion=[30 30];
CalibRegion=Sz;
SliI=12;
HankelTemporalLen=2;
CurFCombined=squeeze(FCombined(:,:,SliI,:));
[~, ~, ~,H]=ghankel(size(CurFCombined,3),HankelTemporalLen,gsize(CurFCombined,1:2));
HCurSli=H*squeeze(Combined(:,:,SliI,:));
% HCurFCombined=H*CurFCombined;
HCurFCombined=fft2cg(HCurSli);
HCurFCombinedP=perm43(HCurFCombined);
% kSize=[6,6];
% kSize=[9,9];
kSize=[12,12];
% kSize=[16,16];
% kSize=[24,24];
kSize=[36,36];
calibME = crop(HCurFCombinedP,CalibRegion(1),CalibRegion(2),HankelTemporalLen,size(CurFCombined,3)-HankelTemporalLen+1);

[M,W] = ME_ESPIRIT(calibME,kSize,eigThresh_1);

maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,HankelTemporalLen]);
mapsB=imresize(maps,Sz);
% tmp=mapsB(:,:,2)./mapsB(:,:,1);
% fgmontagex(angle(tmp))

% M1=squeeze(M(:,:,:,2)./M(:,:,:,1));
M2=squeeze(M(:,:,2,:)./M(:,:,1,:));

% M1a=squeeze(M(:,:,:,2).*M(:,:,:,1));
% M2a=squeeze(M(:,:,2,:).*M(:,:,1,:));

tmp1=imresize(M2(:,:,2),Sz);
tmp=gflip(angle(conj(tmp1)),[]);
% fgmontagex(tmp)
% fgmontagex(tmp,[-pi pi]);title(num2str([CalibRegion kSize],' %d'));

MxHz=500/dTEs_ms(1);
tmp_Hz=tmp/pi*MxHz;

fgmontagex(tmp_Hz,[-MxHz MxHz]);title(num2str([CalibRegion kSize],' %d'));%colorbar

%% Refs
fgmontagex(UpdatedB0Map_Hz(:,:,SliI))
fgmontagex(grmss(s_vals(:,:,:,SliI),3))
%%
CalibRegion=Sz;
HankelTemporalLen=2;
[~, ~, ~,H]=ghankel(size(CurFCombined,3),HankelTemporalLen,gsize(CurFCombined,1:2));

kSzs=[12 16 24 36];

MxHz=500/dTEs_ms(1);
dispstat('','init') 
AllB0_EME=zeros([Sz size(FCombined,3) numel(kSzs)]);
disp(['Starting ' datestr(now)]);
for k=1:numel(kSzs)
    kSize=[1 1]*kSzs(k);
    % SliI=12;
    for SliI=1:size(FCombined,3)
        dispstat(num2str([k SliI]),'timestamp');
        CurFCombined=squeeze(FCombined(:,:,SliI,:));
        HCurSli=H*squeeze(Combined(:,:,SliI,:));
        HCurFCombined=fft2cg(HCurSli);
        HCurFCombinedP=perm43(HCurFCombined);
        
        calibME = crop(HCurFCombinedP,CalibRegion(1),CalibRegion(2),HankelTemporalLen,size(CurFCombined,3)-HankelTemporalLen+1);
        
        M = ME_ESPIRIT(calibME,kSize,eigThresh_1);
        tmp=angle(M(:,:,1,2)./M(:,:,2,2));
        AllB0_EME(:,:,SliI,k)=tmp;
    end
end
disp('Finished');
AllB0_EME=AllB0_EME/pi*MxHz;
save([BaseP FN filesep 'AllB0_EME.mat'],'AllB0_EME');
disp('saved');
%%

tmp_Hz=tmp/pi*MxHz;
fgmontagex(tmp_Hz,[-MxHz MxHz]);title(num2str([CalibRegion kSize],' %d'));%colorbar
%%

%% BART alternative: find the smootherst B that solves. -R D:A:B:C	L2 finite differences
size(HCurSli)
Ops='fmac 0 0';
ScriptFN='/autofs/space/daisy_002/users/Gilad/gUM/fmac00.txt';
WriteLinopToFile(ScriptFN,Ops);

ImSz16=FillOnesTo16(Sz);
Lambdas=10.^(-5:.5:3);
for i=1:numel(Lambdas)
    RecsL(:,:,i)=bart(['picsS -m -S -R D:3:3:' num2str(Lambdas(i)) ' ' ScriptFN],ImSz16,HCurSli(:,:,:,2),HCurSli(:,:,:,1));
end