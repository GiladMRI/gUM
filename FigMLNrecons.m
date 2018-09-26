%%
OnPhantomP='RegridTry3C2_TS_AK_S8P__2018-07-05_18-20-51_train/';
TS='RegridTry3C2_TS_AK_S8__2018-07-04_17-55-39_train/';
SharedTS='RegridTry3C2_TSB_AK_S8__2018-07-04_13-47-38_train/';

BaseBGP='/media/a/f38a5baa-d293-4a00-9f21-ea97f318f647/home/a/TF/srez/';

Paths={OnPhantomP,'OnPhantom'; SharedTS, 'Shared segments Ws.'; TS,'7 Independent segments Ws.' ;
    'RegridTry3C2_1TS_AK_S8__2018-07-05_21-12-11_train/','Single segment'
    'RegridTry3C2_15TS_AK_S8__2018-07-06_01-03-19_train/','Ind., 15 seg.'};
%%
clear ScrNC BatchNC MinNC LastFN
for i=1:size(Paths,1)
    [ScrNC{i},BatchNC{i},MinNC{i},LastFN{i}]=GraphOptFromFolderf([BaseBGP Paths{i,1}]);
end
%%
clear X
for i=1:size(Paths,1)
    X(:,:,:,i)=imread([BaseBGP Paths{i} LastFN{i}]);
end
%%
fgmontage(rot90(X(1:128,128*5+(1:128),1,:)),[0 120],'Size',[1 size(Paths,1)])
XLim=get(gca,'XLim');
NNN=floor(XLim(2))/size(Paths,1);
for i=1:size(Paths,1)
    text((i-1)*NNN+10,20,Paths{i,2},'Color','white','FontSize',12);
end
%%

fgmontage(rot90(X(52:86,128*5+(5:45),1,:)),[0 150],'Size',[1 size(Paths,1)])
XLim=get(gca,'XLim');
NNN=floor(XLim(2))/size(Paths,1);
for i=1:size(Paths,1)
    text((i-1)*NNN+10,20,Paths{i,2},'Color','white','FontSize',12);
end
%%
