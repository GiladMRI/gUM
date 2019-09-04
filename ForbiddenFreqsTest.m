Grads=load('GAll.mat');
CurGrad=Grads.GAll(:,2,1);

GradDwellTime=10;
nGradSteps=numel(CurGrad);
ForbiddenEchospacings=[790 920]; % Prisma: 0.79-0.92ms echospacing forbidden
GradSteps=ForbiddenEchospacings/GradDwellTime;
nForbiddenEchospacings=numel(ForbiddenEchospacings);
GForbidden=zeros(nGradSteps,nForbiddenEchospacings);
for i=1:nForbiddenEchospacings % PE kind of blips
    GForbidden(7:GradSteps(i):end,i)=OneBlockmTm*3;
    GForbidden(6:GradSteps(i):end,i)=OneBlockmTm*2;
    GForbidden(8:GradSteps(i):end,i)=OneBlockmTm*2;
    GForbidden(5:GradSteps(i):end,i)=OneBlockmTm;
    GForbidden(9:GradSteps(i):end,i)=OneBlockmTm;
end

GForbidden=zeros(nGradSteps,nForbiddenEchospacings);
RiseSteps=3;
for i=1:nForbiddenEchospacings % RO kind of blips
    nEffectiveSteps=GradSteps(i)-RiseSteps*2;
    tmp=[ones(RiseSteps,1); ones(nEffectiveSteps,1)*TotalGradAmpmTm/nEffectiveSteps ; ones(RiseSteps,1)];
    tmp=[tmp; -tmp];
    tmp=repmat(tmp,[1000,1]);
    GForbidden(:,i)=tmp(1:nGradSteps);
end
FGForbidden=fft1cg(GForbidden,1); 
FCurGrad=fft1cg([real(CurGrad) imag(CurGrad)],1);
AllI=[real(CurGrad) imag(CurGrad) GForbidden];
AllF=fft1cg(AllI,1);
SAllI=SmoothBySlices(AllI,[40 1],1);
SAllF=fft1cg(SAllI,1);
SAllF(nGradSteps/2+(-10:10),:)=0;

figure;plot(abs(SAllF),'LineWidth',2);hold on;
DD=nGradSteps*GradDwellTime/ForbiddenEchospacings(1);
setXaxis(floor(nGradSteps/2 + [0 1]*2.3*DD));
plot(nGradSteps/2+nGradSteps*GradDwellTime./ForbiddenEchospacings,[50 50],'r','LineWidth',4);
plot(nGradSteps/2+nGradSteps*GradDwellTime/2./ForbiddenEchospacings,[50 50],'Color',[0.5 0 0],'LineWidth',4);
title('Forbidden frequencies analysis');
removeTicks;
%%
All=[abs(FCurGrad) abs(FGForbidden)*100];
SAll=SmoothBySlices(All,[40 1],15);
figure;plot(abs(FCurGrad));hold on;plot(abs(FGForbidden)*100);

figure;bar(1:nGradSteps,[abs(FCurGrad) abs(FGForbidden)*100]);
%% Standard blip: PE
EFOV_mm=200; % 200 mm
Res_mm=0.5; % 2 mm
GAcc=2; % GRAPPA 2
nBlips=EFOV_mm/Res_mm/GAcc; % 50 blips
gamma_MhzT=42.5774806; % Mhz/T
% TotalBlip % Creating pi every res_mm
% AccumulatedPhasePerGradStepOneVoxel=2*pi*GradDwellTime*gamma_MhzT*GradAmpmTm/1000*Res_mm/1000;
TotalGradAmpmTm=pi/(2*pi*GradDwellTime*gamma_MhzT/1000*Res_mm/1000);
TotalGradAmpmTmOneBlip=TotalGradAmpmTm/nBlips;
OneBlockmTm=TotalGradAmpmTmOneBlip/(5+3+1);
%% Standard blip: RO - same but a single blip
