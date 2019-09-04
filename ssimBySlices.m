function [OutA, OutB]=ssimBySlices(A,B)
for i=1:size(A,3)
    [OutA(i), OutB(:,:,i)]=ssim(abs(A(:,:,i)),abs(B(:,:,i)));
end