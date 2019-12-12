[optimizer, metric] = imregconfig('multimodal')

fixed = dicomread('knee1.dcm');
moving = dicomread('knee2.dcm');



fixed = MMM(:,:,3,1);
moving = MMMC(:,:,3,1);

[movingRegistered, RegS] = imregister(moving, fixed, 'rigid', optimizer, metric);


.bashrc
/autofs/cluster/pubsw/arch/CentOS6-x86_64/packages/ANTS/2.3.1/bin
ANTSPATH

mcflirt -in Movingx.nii -cost normcorr
flirt -in Movingx_mcf_mean.nii.gz -ref Ref.nii -cost normcorr -2D -out Movingx_mcf_mean_to_Ref.nii -omat Movingx_mcf_mean_to_Ref.mat

tform = imregtform(moving, fixed, 'rigid', optimizer, metric)
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
figure;imshowpair(fixed, movingRegistered,'Checkerboard')

acosd(tform.T(1,1))


figure;imshowpair(fixed, moving,'Checkerboard')

fixed=RefSV1x(:,:,Ord(SliIs));
for Idx=5:8
    disp(Idx);
moving=Rec_CompgB0_RSS_MXCMY(:,:,SliIs,Idx);
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

movingRegistered = imregister(moving, fixed, 'rigid', optimizer, metric);
MMMCC{Idx}=movingRegistered ;
end
disp('ok');
% figure;imshowpair(fixed, movingRegistered,'Scaling','joint')
figure;imshowpair(fixed, movingRegistered,'Checkerboard')


MMM=repmat(Rec_CompgB0_RSS_MXCMY(:,:,SliIs,:),[1 1 5 1]);
RRR=repmat(RefSV1x(:,:,Ord(SliIs)),[1 1 5 1]);


MMMC4=CombineDims(CombineDims(PartitionDim(cat(3,MMMCC{5:8}),3,2),[3 1]),[3 2]);
