function [ alpha_thresh ] = TV_thresh_rotationInvariant16D(alpha, lambda)
% function [ alpha_thresh ] = TV_thresh(alpha, lambda)

TValpha=TVOP_MSlice16D*alpha;
TValphaAbs=grms(TValpha,17);
Direction=TValpha./max(TValphaAbs,eps);
% NewNorm=SoftThresh(TValphaAbs,lambda);
NewNorm=max(TValphaAbs-lambda,0);
alpha_thresh=TVOP_MSlice16D'*(NewNorm.*Direction);
% alpha_thresh = TVOP'*SoftThresh(TVOP*alpha, lambda);