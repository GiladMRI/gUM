function [ alpha_thresh ] = TV_thresh16D(alpha, lambda)
% function [ alpha_thresh ] = TV_thresh(alpha, lambda)

% TValpha=TVOP_MSlice16D*alpha;
% NewNorm=SoftThresh(TValphaAbs,lambda);
% NewNorm=max(TValphaAbs-lambda,0);
% alpha_thresh=TVOP_MSlice16D'*(NewNorm.*Direction);
alpha_thresh = TVOP_MSlice16D'*SoftThresh(TVOP_MSlice16D*alpha, lambda);