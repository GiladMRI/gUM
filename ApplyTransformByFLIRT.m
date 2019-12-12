function Out=ApplyTransformByFLIRT(Imgs,TMats)

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

Out=zeros(size(Imgs));
for i=1:nImgs
    MovFN='tmp';
    RefFN=MovFN;
    Raw2Nii(Imgs(:,:,i),[MovFN '.nii'],'float32');
    putLines('tmp.mat',gmat2cell(num2str(TMats(:,:,i),'%f '),1));
    % Cmd=['flirt -in ' RefFN '.nii.gz -ref ' RefFN '.nii -applyxfm -init tmp.mat -out FlirtOut.nii'];
    Cmd=['flirt -in ' MovFN '.nii.gz -ref ' RefFN '.nii -applyxfm -init tmp.mat -out FlirtOut.nii'];
    system(Cmd);
    Out(:,:,i)=loadniidata('FlirtOut.nii.gz');
end
