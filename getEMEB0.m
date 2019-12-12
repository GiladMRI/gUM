function Out=getEMEB0(In,kSize,dT_ms,eigThresh_1)
if(nargin<4)
    eigThresh_1=0.02;
end
dispstat('','init');
Sz=gsize(In,1:2);
if(numel(kSize)<3)
    CalibRegion=Sz;
else
    CalibRegion=kSize(3:4);
end
kSize=kSize(1:2);
HankelTemporalLen=2;
nEchoes=size(In,3);
[~, ~, ~,H]=ghankel(nEchoes,HankelTemporalLen,Sz);

nSlices=size(In,4);
Out=zeros([Sz nSlices]);
for i=1:nSlices
    dispstat(['Slice ' num2str(i)],'timestamp')
    HCurSli=H*In(:,:,:,i);
    HCurFCombined=fft2cg(HCurSli);
    HCurFCombinedP=perm43(HCurFCombined);
    
    calibME = crop(HCurFCombinedP,CalibRegion(1),CalibRegion(2),HankelTemporalLen,nEchoes-HankelTemporalLen+1);
    
    M = ME_ESPIRIT(calibME,kSize,eigThresh_1);
    Out(:,:,i)=imresize(angle(M(:,:,1,2)./M(:,:,2,2)),Sz);
end
Out=Out*500/pi/dT_ms;