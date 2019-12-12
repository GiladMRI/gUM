FigsP='/autofs/cluster/kawin/Gilad/Figs/';
mkdir(FigsP);
system(['chmod +777 -R ' FigsP]);
disp([FigsP ' Created']);
%%
FigP=[FigsP 'dB0Base' filesep]; 
mkdir(FigP);
system(['chmod +777 -R ' FigP]);
disp([FigP ' Created']);
%%
% save([FigP 'B0MapToUse_Hz.mat'],'B0MapToUse_Hz');
load([FigP 'B0MapToUse_Hz.mat'],'B0MapToUse_Hz');


Msk=imfillholesBySlices(RefSvals1r>0.0001);
%%
Thickness_mm=3;
DistBetweenSlices_mm=3;
TrgThickness_mm=2;
%
dB0dx=symD(B0MapToUse_Hz,1);
dB0dy=symD(B0MapToUse_Hz,2);
dB0dz=symD(B0MapToUse_Hz,3)*TrgThickness_mm/DistBetweenSlices_mm;
%%
SliIToShow=[2 6 15];
CurB0=B0MapToUse_Hz(:,:,SliIToShow);
CurdB0=cat(4,dB0dx(:,:,SliIToShow),dB0dy(:,:,SliIToShow),dB0dz(:,:,SliIToShow));
CurxB0=cat(4,RefSvals1r(:,:,SliIToShow)*8e4,CurB0,CurdB0*5);
MskToShow=cat(3,RefSvals1r(:,:,SliIToShow)>0.00015);
MskToShowx=imfillholesBySlices(MskToShow);
fgmontagex(CurxB0.*MskToShow,[-300 300]);
%%
cmap=[[0 0 0];jet(256)];
cmap=[[0 0 0]+0.5;jet(256)];
fgmontagex((CurxB0+300).*MskToShowx,[0 600]);colormap(cmap)
%%
fgmontage((1:256).');colormap(jet(256))
%%
fgmontagex((CurxB0(:,:,:,1)+300).*MskToShow,[300 600],'Size',[3 1]);colormap(gray)

% RGBMg=repmat(CurxB0(,[1 1 3]);
%%
EBaseP =    '/autofs/cluster/kawin/Gilad/EPTI_and_spi68ms_on_CL/';
load([EBaseP 'AllRecs4.mat'],'AllRecs4');
load([EBaseP 'RecSubspaceTDPEROS.mat'],'RecSubspaceTDPEROS');
%%
QQa=abs(perm43(squeeze(abs(AllRecs4(:,:,55,:,[1 2 4]))-abs(AllRecs4(:,:,55,:,3)))));
QQa=abs(perm43(squeeze(abs(AllRecs4(:,:,55,:,[2:4]))-abs(AllRecs4(:,:,55,:,1:3)))));

% SliX=8;
SliX=6;
SliX=2;
SliX=16;
fgmontagex(AllRecs4(:,:,55,SliX,:));%caxis(caxis/2)
figure;imagesc(QQa(:,:,:,SliX));removeTicks;daspect([1 1 1]);
fgmontagex(SRefB0r(:,:,SliX));
%%
M=cat(4,QQa(:,:,:,SliX),abs(repmat(perm54(AllRecs4(:,:,55,SliX,:))/10,[1 1 3])));
M=cat(4,M,abs(repmat(perm54(AllRecs4(:,:,25,SliX,[1 2 4]))/10,[1 1 3])));
% M=CombineDims(M,[4 2]);
M=CombineDims(CombineDims(PartitionDim(M(:,:,:,[2 3 4 5 1 6 7 8]),4,2),[5 1]),[4 2]);
figure;imagesc(M);removeTicks;daspect([1 1 1]);
%%
M=cat(4,QQa(:,:,:,SliX),abs(repmat(perm54(AllRecs4(:,:,55,SliX,:))/10,[1 1 3])));
M=CombineDims(M(:,:,:,[1 2 5]),[4 2]);
figure;imagesc(M);removeTicks;daspect([1 1 1]);
%%
