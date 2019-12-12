function [M,W]=ME_ESPIRIT(kCalib,ksize,eigThresh_1)
% kCalib : [x y channels contrasts]
% dat2Kernel
[sx,sy,nc,ne] = size(kCalib);
imSize = [sx,sy] ;

for i=1:ne
    tmpC{i} = im2row(kCalib(:,:,:,i),ksize);
end
tmpM=cat(1,tmpC{:});
[tsx,tsy,tsz] = size(tmpM);
A=reshape(tmpM,tsx,tsy*tsz);

[U,S,V] = svd(A,'econ');
    
kernel = reshape(V,ksize(1),ksize(2),nc,size(V,2));
S = diag(S);S = S(:);

idx = max(find(S >= S(1)*eigThresh_1));
[M,W] = kernelEig(kernel(:,:,:,1:idx),[sx,sy]);