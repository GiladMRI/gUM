function Out=ApplyTransformAsFLIRT(Imgs,TMats)
if(size(TMats,1)==4)
    TMats=TMats([1 2 4],[1 2 4],:,:,:,:,:,:,:,:);
end
if(TMats(1,3)~=0 || TMats(2,3)~=0) % FLIRT format
    TMats=perm21(TMats(1:3,1:3,:,:,:,:,:,:,:,:));
    TMats(3,[2 1],:)=TMats(3,[1 2],:);
    TMats(3,2,:)=-TMats(3,2,:);
    TMats(3,3,:)=1;
end

nImgs=size(Imgs,3);
nTMats=size(TMats,3);
if(nImgs==1)
    Imgs=repmat(Imgs,[1 1 nTMats]);
end
if(nTMats==1)
    TMats=repmat(TMats,[1 1 nImgs]);
end
nImgs=size(Imgs,3);
nTMats=size(TMats,3);

Sz=gsize(Imgs,1:2);
RefSpace=imref2d(Sz);

tX = mean(RefSpace.XWorldLimits);
tY = mean(RefSpace.YWorldLimits);

tTranslationToCenterAtOrigin = [1 0 0; 0 1 0; -tX -tY,1];
tTranslationBackToOriginalCenter = [1 0 0; 0 1 0; tX tY,1];

Out=zeros(size(Imgs));
for i=1:nImgs
    tmpMatM=tTranslationToCenterAtOrigin*TMats(:,:,i)*tTranslationBackToOriginalCenter;
    tform = affine2d(tmpMatM);
    Out(:,:,i) = imwarp(Imgs(:,:,i),tform,'OutputView',RefSpace);
end
