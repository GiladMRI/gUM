%% SCC
% 1ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-33-46_train/';
%%
D=dir([VP 'Tra*.mat']);
D=D([D.bytes]>1000);
Q=load([VP D(end).name]);
disp(D(end).name);

G_LossV=Q.G_LossV;
var_list=Q.var_list;
Q=rmfield(Q,{'G_LossV','var_list'});

Q=CombineRIFlds(Q);

Flds=fieldnames(Q);
SFlds=sort(Flds);

for i=1:numel(SFlds)
    disp([PadStringWithBlanks(SFlds{i},65) num2str(size(Q.(SFlds{i})),'% 9d         ')]);
end
%%
MapsForTF=load('/media/a/H1/maps128x128x8.mat');
X=squeeze(Q.gene_GEN_L004_C2D_weightC);
OutCh=permute(MultMatTensor(X.',permute(MapsForTF.maps,[3 1 2])),[2 3 1]);
%%
ShowAbsAngle(MapsForTF.maps,[],'Size',[2 4])
subplot(1,2,1);title('Maps of 8 channel used');
%%
ShowAbsAngle(OutCh)
subplot(1,2,1);title('SCC to 1ch');
%%
% 2ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-37-15_train/';
ShowAbsAngle(OutCh)
subplot(1,2,1);title('SCC to 2ch');
%%
% 3ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-41-18_train/';
ShowAbsAngle(OutCh,[],'Size',[1 3])
subplot(1,2,1);title('SCC to 3ch');
%%
% 5ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-28-30_train/';
ShowAbsAngle(OutCh,[],'Size',[1 5])
subplot(1,2,1);title('SCC to 5ch');
%% GCC
% 2ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC2__2018-06-29_16-49-45_train/';
% 3ch
VP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC3__2018-06-29_17-00-15_train/';
%%
X=squeeze(Q.gene_GEN_L004_einsum_weightC);
% Y=permute(X,[1 4 2 3]);
Y=permute(X,[4 1 2 3]);
Z=Y.*MapsForTF.maps;
F=squeeze(sum(Z,3));
ShowAbsAngle(F,[],'Size',[1 size(F,3)]);subplot(1,2,1);title(['GCC to ' num2str(size(F,3)) ' channels']);
%%
BaseBGP='/media/a/DATA/DataForFigCC/';
D=dir(BaseBGP);
D=D([D.isdir]);
D=D(3:end);

PC={'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-33-46_train/','SCC-1';    
'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-37-15_train/', 'SCC-2';
'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-41-18_train/','SCC-3';
'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_SCC__2018-06-29_16-28-30_train/','SCC-5';
'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC2__2018-06-29_16-49-45_train/','GCC-2';
'/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/SMASH_GCC3__2018-06-29_17-00-15_train/','GCC-3'};
%%
clear ScrNC BatchNC MinNC
% for i=1:numel(D)
%     [ScrNC{i},BatchNC{i},MinNC{i}]=GraphOptFromFolderf([BaseBGP D(i).name filesep]);
% end
for i=1:size(PC,1)
    [ScrNC{i},BatchNC{i},MinNC{i},LastFN{i}]=GraphOptFromFolderf(PC{i,1});
end
%%
figure;
Clr={'r','r','r','r','b','b'};
Mrkr={'-','--','--','-','-','--'};
LW=[1,1,2,2,1,1]*2;
for i=1:6
%     plot(BatchNC{i},ScrNC{i},'.-');hold on;
    plot(MinNC{i},ScrNC{i},[Clr{i} '*' Mrkr{i}],'LineWidth',LW(i));hold on;
end
% legend({D.name})
legend(PC{:,2})
xlabel('Minutes');
ylabel('Score');
axis([-0.1636    3.0576   -0.0067    0.0886]);
%%
clear X
for i=1:6
    X(:,:,:,i)=imread([PC{i,1} LastFN{i}]);
end
%%
fgmontage(X(1:128,128*5+(1:128),1,:),[],'Size',[1 6])
%%
fgmontage(X(1:128,128*5+(1:128),1,:),[0 20],'Size',[1 6])
for i=1:6
    text(i*320-210,160,PC{i,2},'Color','black','FontSize',14);
end