QQ=load('/home/deni/TO_FromGilad/E1/RecTSAll.mat');

QQ=load('/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/Sep19_Deni/meas_MID543_gBP_ep2d_bold_multiecho_ASL_SMS_Spic_4min_FID30394/Sparse_AllRep.mat');
Sparse_AllRep=QQ.ImResR;
niftiwrite(abs(Sparse_AllRep),'DK_Sparse_E1.nii');



Combined=CombineDims(QQ.RecTSAll,[3 4]);
niftiwrite(abs(Combined),'TO_Bart_E1.nii');

QQ=load('/home/deni/TO_FromGilad/E2/RecTSAll.mat');
Combined=CombineDims(QQ.RecTSAll,[3 4]);
niftiwrite(abs(Combined),'TO_Bart_E2.nii');

% function PushLocIntoNifti
%%
FN='/home/deni/gUM/TO_Bart_E1.nii';
FN='/home/deni/gUM/DK_Sparse_E1.nii';
NII=load_untouch_nii(FN);
hdr=NII.hdr;

NRO=size(NII.img,1);
NPE=size(NII.img,2);

ROFov=AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
PixDimRO=ROFov/NRO;

PEFov=AData.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV;
PixDimPE=ROFov/NPE;

SliceThickness=AData.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;

TR=AData.hdr.Phoenix.alTR{1}/1e6;

hdr.dime.pixdim(1:4)=[PixDimRO PixDimPE SliceThickness TR];

CenterPos=mean(SlbLoc(:,nSlices/2+[0 1]),2);
hdr.hist.srow_x=[RotMat(1,:) CenterPos(1)];
hdr.hist.srow_y=[RotMat(2,:) CenterPos(2)];
hdr.hist.srow_z=[RotMat(3,:) CenterPos(3)];

CenterPosQ=mean(Qoffset(:,nSlices/2+[0 1]),2);
hdr.hist.qoffset_x=CenterPosQ(1);
hdr.hist.qoffset_y=CenterPosQ(2);
hdr.hist.qoffset_z=CenterPosQ(3);

hdr.hist.quatern_b=AData.image.slicePos(5,1);
hdr.hist.quatern_c=AData.image.slicePos(6,1);
hdr.hist.quatern_d=AData.image.slicePos(7,1);



FNOut=[FN(1:end-4) 'L.nii'];
save_untouch_nii(NII,FNOut);
disp(['Saved ' FNOut]);