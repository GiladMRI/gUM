%ImResRX=ImResR(:,:,SliIs,:);
ImResRX=ImResR;
ImResRX=RecTSAllX;
%% Perfusion
mImResRX=mean(ImResRX,4);
GMsk=abs(mImResRX)>2e-3;
%%
Even=ImResRX(:,:,:,2:2:end);
Odd=ImResRX(:,:,:,1:2:end);
Diff=Odd-Even;
mDiff=mean(Diff,4);
sDiff=std(Diff,1,4);
tDiff=mDiff./sDiff;
fgmontage(tDiff);
%%
Timepoints_sec=(1:nReps)*AData.hdr.Meas.alTR(1)/1e6;
StimOn=mod(Timepoints_sec,60)>40;
%%
Timepoints_sec2=(1:(nReps/2))*AData.hdr.Meas.alTR(1)*2/1e6;
StimOn2=mod(Timepoints_sec2-3,60)>40;
%% ttest
Sli=13;
x=50;
y=50;
tic
for Sli=6 % :24
    disp(Sli);
    for x=1:size(Diff,1)
        for y=1:size(Diff,2)
            [h(x,y,Sli),p(x,y,Sli)] = ttest2(squeeze(Diff(x,y,Sli,StimOn2)),squeeze(Diff(x,y,Sli,~StimOn2)),'Vartype','unequal');
        end
    end
end
toc
pMsked=-log(p).*GMsk;

TMsk=pMsked>1.5;
TMsk(:,:,setdiff(1:24,[6:9]))=0;
TMsk(:,60:end,:)=0;
fgmontage(TMsk)

TwoD=Reshape4d22d(Diff,TMsk);
figure;plot(mean(TwoD,1))
hold on
plot(StimOn2*5e-5)
%% BOLD
ImResRX_On=ImResRX(:,:,:,StimOn);
ImResRX_Off=ImResRX(:,:,:,~StimOn);

mImResRX_On=mean(ImResRX_On,4);
mImResRX_Off=mean(ImResRX_Off,4);

rmImResRX=(mImResRX_Off-mImResRX_On)./std(ImResRX_Off,[],4);
%%
TMsk=rmImResRX>0.8;
TwoD=Reshape4d22d(ImResRX,TMsk);
figure;plot(mean(abs(TwoD),1)-0.16)
hold on
plot(circshift(StimOn,-6,2)*5e-2)

%% ASL
Labeled=mod(1:nReps,2)==1;

ImResRX_OnL=ImResRX(:,:,:,StimOn & Labeled);
ImResRX_OnNL=ImResRX(:,:,:,StimOn & ~Labeled);

ImResRX_OffL=ImResRX(:,:,:,~StimOn & Labeled);
ImResRX_OffNL=ImResRX(:,:,:,~StimOn & ~Labeled);

mImResRX_OnL=mean(ImResRX_OnL,4);
mImResRX_OnNL=mean(ImResRX_OnNL,4);
mImResRX_OffL=mean(ImResRX_OffL,4);
mImResRX_OffNL=mean(ImResRX_OffNL,4);

mSigPerfOn=mImResRX_OnL-mImResRX_OnNL;
mSigPerfOff=mImResRX_OffL-mImResRX_OffNL;

mSigDiffPerf=mSigPerfOn-mSigPerfOff;

fgmontage(mSigDiffPerf);
%%
TMsk=mDiff*0;
TMsk(30:end-40,17:57,13)=mDiff(30:end-40,17:57,13)>1e-5;

TwoD=Reshape4d22d(Diff,TMsk);


figure;plot(mean(TwoD,1))
hold on
plot(StimOn*5e-5)