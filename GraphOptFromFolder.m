BaseTFRes='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';
Prefix='RegridTry1C2_TS';
Prefix='RegridTry3C2_TSB';
% Prefix='RegridTry1C2_TS2';
DB=dir([BaseTFRes Prefix '*']);
[SDB,DOrd]=sort([DB.datenum]);
LastDir=DB(DOrd(end)).name;
%
% BaseGP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/RegridTry1C2_TS__2018-06-27_15-30-21_train/';
BaseGP=[BaseTFRes LastDir filesep];
D=dir([BaseGP '*.png']);
Datenums=[D.datenum];
DNames={D.name};
[SDatenums, Ord]=sort(Datenums);
SDNames=DNames(Ord);
SDNames=SDNames(4:end);
SDatenums=SDatenums(4:end);

X=cat(1,SDNames{:});
X=textscan(X.','batch%06d_out_%f.png');
BatchN=X{1};
ScrN=X{2};
MinN=(SDatenums-SDatenums(1))*86400/60;

figure;
subplot(1,2,1);
plot(BatchN,ScrN,'.-');xlabel('batch #');ylabel('Scr');
% subplot(1,3,2);
% plot(MinN,ScrN,'.-');xlabel('Minutes');ylabel('Scr');
ax1 = gca; % current axes
ax=axis();
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','origin',...
    'Color','none');

line(MinN,ScrN*0,'Parent',ax2,'Color','k')
setYaxis(ax(3:4));
xlabel('Minutes')
subplot(1,2,2);
ParamsFNX=[BaseGP 'ParamsUsed.txt'];
lines = getLines(ParamsFNX);
plot([0 1],[0 1],'.')
for i=1:numel(lines)
    text(0.1,1-0.03*i,lines{i},'Interpreter','none');
end
title(LastDir,'Interpreter','None')
%%
I=imread([BaseTFRes LastDir filesep SDNames{end}]);
fgmontage(I(1:128,128*5+(1:128),1))
NNRes=double(I(1:128,128*5+(1:128),1))/(255*RealDataFac);
%%
clear AllRes
AllRes(:,:,:,1)=cat(3,Rec,RecMM(:,:,1));
AllRes(:,:,:,2)=cat(3,resL1ESPIRiT(:,:,1), resL1ESPIRiTCC(:,:,1));
AllRes(:,:,3,1)=NNRes;
AllRes(:,:,3,2)=PerfMap;
AllRes=gflip(permute(AllRes,[2 1 3 4]),1);
fgmontage(AllRes)
%%
NNResR=rot90(NNRes);
PerfMapR=rot90(PerfMap);
imB = NNResR;                       % Background image
imF = PerfMapR;   % Foreground image
% function [hF,hB] = imoverlay(B,F,climF,climB,cmap,alpha,haxes)
% [hf,hb] = imoverlay(imB,imF,[40,180],[0,0.6],'jet',0.6);
[hf,hb] = imoverlayx(imB,imF,[0 7e-3],[0.0012 0.009],'parula',1) %,'jet',PerfMapR>0.002);
% colormap('parula'); % figure colormap still applies
%%
PerfMapR=rot90(PerfMap);

PerfMapR(NNResR<0.001)=0;
ClimF=[0.0012 0.009];
A=(min(max(PerfMapR,ClimF(1)),ClimF(2))-ClimF(1))/(ClimF(2)-ClimF(1));
RGB=ind2rgb(round(A*255)+1,parula(256));
Msk=PerfMapR>ClimF(1);

ClimB=[0 7e-3];
B=(min(max(NNResR,ClimB(1)),ClimB(2))-ClimB(1))/(ClimB(2)-ClimB(1));
X=RGB.*Msk+B.*(1-Msk);
% figure;imshow(X)

AllRes3=repmat(AllRes,[1 1 1 1 3]);
AllRes3(:,:,3,2,:)=X*ClimB(2);

AllResI=CombineDims(AllRes3,[3 2]);
AllResI=CombineDims(AllResI,[3 1]);
figure;imshow(abs(AllResI)/ClimB(2))

FS=11;
Clr=[1 1 1];
text(10,10,'BART 1 map','Color',Clr,'FontSize',FS)
text(128*1+0,10,'BART 1 map','Color',Clr,'FontSize',FS)
text(10,128*1+10,'ESPIRIT 2 maps, 15TS','Color',Clr,'FontSize',FS)
text(128*1+0,128*1+10,'ESPIRIT 2 maps, 15TS, CC->13','Color',Clr,'FontSize',FS)

text(128*2+0,10,'Linear net, 7TS','Color',Clr,'FontSize',FS)
text(128*2+0,128*1+10,'Perfusion, 7TS','Color',Clr,'FontSize',FS)

% text(100,100,'BART')
%%
% AllResX=repmat(AllRes,[1 
% %%
% figure;
% TheGrayscaleImage=NNRes;
% TheColorImage=PerfMap;
% YourTransparencyData=PerfMap>0.002;
% imagesc(TheGrayscaleImage);
% colormap(gray(256));
% hold on
% imagesc(TheColorImage, 'AlphaData', YourTransparencyData)
%%
fgmontage(cat(3,Rec,RecMM(:,:,1)));
xlabel('Left - Single map, ADMM. Right - ESPIRIT 2 maps');

title(['BART No B0 W=' num2str(Lambda)]);
YLbl=['Sli' num2str(SliI,'%02d')];
ylabel(YLbl);

% fgmontage(im_res,[0 7e-3]);
% XFMStrFull=['[' XFMStr ',' num2str(filterSize) ',' num2str(wavScale) ',' num2str(xfmWeight) ']'];
% %     XFMStr ' Size=' num2str(filterSize) ' Scale=' num2str(wavScale) ' W=' num2str(xfmWeight)
% XLbl=['L' num2str(param.pNorm) ',TVW=' num2str(param.TVWeight) ',' XFMStrFull ];
% xlabel(XLbl)
% YLbl=['Sli' num2str(SliI,'%02d')];
% ylabel(YLbl);

fgmontage(resL1ESPIRiT(:,:,1),[0 7e-3])
title('resL1ESPIRiT 2 maps, B0');
ylabel(YLbl);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);

fgmontage(cat(3,resL1ESPIRiT(:,:,1), resL1ESPIRiTCC(:,:,1)),[0 7e-3])
title(['resL1ESPIRiT 2 maps, B0. Right - with CC -> ' num2str(ncc)]);
ylabel(YLbl);
xlabel(['Daubechies_TI lam ' num2str(lam) ' splitWeight ' num2str(splitWeight)]);
