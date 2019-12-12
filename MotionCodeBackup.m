%% get center-mass
[OutA,OutB]=system(['fslstats -t Movingx.nii -C']);
%% Check applyting flirt transform
% Tmean=loadniidata([MovFN '_mcf_mean.nii.gz']);
% Tmean=Tmean(:,:,3);
% 
% [~,TmeanCMStr]=system(['fslstats ' MovFN '_mcf_mean.nii.gz -C']);
% TmeanCM=str2num(TmeanCMStr);
% 
% TmeanRot=loadniidata([MovFN '_mcf_mean_to_' RefFN '.nii.gz']);
% TmeanRot=TmeanRot(:,:,3);
% 
% RefSpace=imref2d(Sz);
% tXb = mean(RefSpace.XWorldLimits);
% tYb = mean(RefSpace.YWorldLimits);
% tTranslationToCenterAtOriginb = [1 0 0; 0 1 0; -tXb -tYb,1];
% tTranslationBackToOriginalCenterb = [1 0 0; 0 1 0; tXb tYb,1];
% 
% tXc = -1;
% tYc = -120;
% 
% tTranslationToCenterAtOriginc = gTransMat2D([-tXc -tYc]);
% tTranslationBackToOriginalCenterc = gTransMat2D([tXc tYc]);
% 
% TMatx=TMat(1:3,1:3).';
% % TMatx(3,1:2)=TMat([1 2],4);
% TMatx(3,1)=-TMat(1,4);
% TMatx(3,2)=-TMat(2,4);
% TMatx(3,1)=TMat(1,4);
% 
% % TMatx(3,1)=-0.5; % rightward
% % TMatx(3,2)=-0.5; % downward
% % TMatx(3,1)=0; % rightward
% % TMatx(3,2)=0; % downward
% TMatx(3,3)=1;
% TMaty=tTranslationToCenterAtOrigin*TMatx*tTranslationBackToOriginalCenter;
% % TMaty=tTranslationToCenterAtOriginb*TMatx*tTranslationBackToOriginalCenterb;
% % TMaty=tTranslationToCenterAtOriginc*TMatx*tTranslationBackToOriginalCenterc;
% tformx = affine2d(TMaty);
% 
% WPos=[680   553   560   420];
% movingRegistered = imwarp(Tmean,tformx,'OutputView',RefSpace);
% fgmontagex(movingRegistered);set(gcf,'Position',WPos)
% 
% fgmontagex(Tmean);set(gcf,'Position',WPos)
% fgmontagex(TmeanRot);set(gcf,'Position',WPos)
%% Apply mcflirt using some matrix
% tmpMat=eye(4);
% tmpMat(1:3,1:3)=grotzd(13);
% tmpMat(1,4)=15; % positive is upward in FLIRT
% tmpMat(2,4)=10; % positive is rightward in FLIRT
% tmpMat(4,4)=1;
% Tf=ApplyTransformByFLIRT(ForRef,tmpMat);
% Tm=ApplyTransformAsFLIRT(ForRef,tmpMat);
% 
% tmpMatM=gTMatFLIRTtoM(tmpMat);
% Tm3=ApplyTransformAsFLIRT(ForRef,tmpMatM);
% 
% 
% % tmpMatx=tTranslationToCenterAtOriginbc*tmpMatM*tTranslationBackToOriginalCenterbc;
% tmpMatx=tTranslationBackToOriginalCenterbc*tmpMatM*tTranslationToCenterAtOriginbc;
% Tmc=ApplyTransformFromCenter(ForRef,tmpMatx);
% 
% 
% % % tmpMat=TMat;
% % tmpMat=TMats(:,:,1)*TMat;
% % putLines('tmp.mat',gmat2cell(num2str(tmpMat,'%f '),1));
% % % Cmd=['flirt -in ' RefFN '.nii.gz -ref ' RefFN '.nii -applyxfm -init tmp.mat -out FlirtOut.nii'];
% % Cmd=['flirt -in ' MovFN '.nii.gz -ref ' RefFN '.nii -applyxfm -init tmp.mat -out FlirtOut.nii'];
% % system(Cmd);
% % FlirtOut=loadniidata('FlirtOut.nii.gz');
% % 
% % fgmontagex(FlirtOut(:,:,3,1));set(gcf,'Position',WPos);title('FLIRT move');
% % 
% % tmpMat=MultTensorMat1(TMats,TMat);
% % 
% % tmpMatM=perm21(tmpMat(1:3,1:3,:));
% % tmpMatM(3,2,:)=-tmpMat(1,4,:);
% % tmpMatM(3,1,:)=tmpMat(2,4,:);
% % tmpMatM(3,3,:)=1;
% % % tX = -RefCM(1);
% % % tY = -RefCM(2);
% % tX = -1;
% % tY = -120;
% % 
% % tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
% % tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
% % 
% % % tmpMatM=tTranslationToCenterAtOriginb*tmpMatM*tTranslationBackToOriginalCenterb;
% % % tmpMatM=tTranslationToCenterAtOrigin*tmpMatM*tTranslationBackToOriginalCenter;
% % % tmpMatM=tTranslationBackToOriginalCenter*tmpMatM*tTranslationToCenterAtOrigin;
% % tmpMatM=MultMatTensor(tTranslationBackToOriginalCenter,MultTensorMat1(tmpMatM,tTranslationToCenterAtOrigin));
% % clear tform movingRegistered
% % for i=1:size(tmpMatM,3)
% %     tform(i) = affine2d(tmpMatM(:,:,i));
% % end
% % for i=1:size(tmpMatM,3)
% %     movingRegistered(:,:,i) = imwarp(ToMove(:,:,1,i),tform(i),'OutputView',RefSpace);
% % %     movingRegistered(:,:,i) = imwarp(ForRef,tform(i),'OutputView',RefSpace);
% % end
% % fgmontagex(movingRegistered);set(gcf,'Position',WPos);title('MATLAB move');
%%
% fgmontagex(ForRef);set(gcf,'Position',WPos);title('Base');
% 
% [~,RefCMStr]=system(['fslstats ' RefFN '.nii.gz -C']);
% RefCM=str2num(RefCMStr);

%%
% In=squeeze(ToMove);
% 
% % fixed=ForRef;
% [optimizer, metric] = imregconfig('multimodal');
% % [optimizer, metric] = imregconfig('monomodal');
% optimizer = registration.optimizer.OnePlusOneEvolutionary
% optimizer.InitialRadius = 0.009;
% optimizer.Epsilon = 1.5e-4;
% optimizer.GrowthFactor = 1.01;
% optimizer.MaximumIterations = 300;
% 
% dispstat('','init')
% for Idx=1:size(In,3)
%     dispstat(['Shot ' num2str(Idx)],'timestamp');
%     moving=In(:,:,Idx);
%     tform = imregtform(moving, Moved(:,:,Idx), 'rigid', optimizer, metric);
%     movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
%     TMats(:,:,Idx)=tform.T;
%     Out(:,:,Idx)=movingRegistered;
% end
% disp('ok');
%%
Rec_CompgB0_RSS_MXCT=grmss(cat(3,Rec_CompgB0_RSS_MXC{:}),4);
BaseI=Rec_CompgB0_RSS_MXCT(:,:,1);
TBase=eye(3);
TBase=gTransMat2D(-dy_vox,-dx_vox);
%%
% DxTtl_voxm=[dx_vox  0   dx_vox+7    0           dx_vox  0       dx_vox+7    0           dx_vox+7    dx_vox+7    ];
% DyTtl_voxm=[dy_vox  0   0           dy_vox+5    dy_vox  0       0           dy_vox+5    dy_vox+5    dy_vox+5    ];
% RotsDm    =[0       0   0           0           7.5     7.5     7.5         7.5         7.5         -15         ];

TrgI=Rec_CompgB0_RSS_MXCT(:,:,10);
% TMat=eye(3);
TMat=grotzd(-15); % positive is anti-clockwise
TMat(3,1)=dy_vox+5; % positive is rightward
TMat(3,2)=dx_vox+7; % positive is downward

% TMat=tTranslationToCenterAtOriginb*TBase*TMat*tTranslationBackToOriginalCenterb;
% 
% % tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
% % tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
% % tmpMatM=tTranslationToCenterAtOriginb*tmpMatM*tTranslationBackToOriginalCenterb;
% % tmpMatM=tTranslationToCenterAtOrigin*tmpMatM*tTranslationBackToOriginalCenter;
% % tmpMatM=tTranslationBackToOriginalCenter*tmpMatM*tTranslationToCenterAtOrigin;
% % tmpMatM=MultMatTensor(tTranslationBackToOriginalCenter,MultTensorMat1(tmpMatM,tTranslationToCenterAtOrigin));
% clear tform movingRegistered
% tform = affine2d(TMat);
% movingRegistered = imwarp(BaseI,tform,'OutputView',RefSpace);

movingRegistered = ApplyTransformFromCenter(BaseI,TBase*TMat);

fgmontagex(movingRegistered)



%%
%%
% TMat=grotz(-7.1569*pi/180);
TMat=grotz(7.5*pi/180);
% TMat=grotz(0*pi/180);
% TMat(3,1)=7; % positive is rightward
% TMat(3,2)=5; % positive is downward
RefSpace=imref2d(Sz);

tX = mean(RefSpace.XWorldLimits);
tY = mean(RefSpace.YWorldLimits);
tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];
tDown=[1 0 0; 0 1 0; 0 5 1];
tRight=[1 0 0; 0 1 0; 7 0 1];
tform1 = affine2d(TMat);
tform1 = affine2d(tTranslationToCenterAtOrigin*TMat*tTranslationBackToOriginalCenter*tRight);

% RefSpace.XWorldLimits=[-60 60];
% RefSpace.YWorldLimits=[-60 60];
% RefSpace.XIntrinsicLimits=[-60 60];
% RefSpace.YIntrinsicLimits=[-60 60];
% tform = imregtform(moving, fixed, 'rigid', optimizer, metric)
movingRegistered = imwarp(grmss(Rec_CompgB0_RSS_MX0,4),tform1,'OutputView',RefSpace);
fgmontagex(movingRegistered)