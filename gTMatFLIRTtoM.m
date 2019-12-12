function TMats=gTMatFLIRTtoM(TMats)
TMats=TMats([1 2 4],[1 2 4],:,:,:,:,:,:,:,:);
TMats=perm21(TMats(1:3,1:3,:,:,:,:,:,:,:,:));
TMats(3,[2 1],:)=TMats(3,[1 2],:);
TMats(3,2,:)=-TMats(3,2,:);
TMats(3,3,:)=1;
