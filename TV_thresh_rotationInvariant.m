function [ alpha_thresh ] = TV_thresh_rotationInvariant(alpha, lambda)
% function [ alpha_thresh ] = TV_thresh(alpha, lambda)

TValpha=TVOP*alpha;
TValphaAbs=grms(TValpha,3);
Direction=TValpha./max(TValphaAbs,eps);
% NewNorm=SoftThresh(TValphaAbs,lambda);
NewNorm=max(TValphaAbs-lambda,0);
alpha_thresh=TVOP'*(NewNorm.*Direction);
% alpha_thresh = TVOP'*SoftThresh(TVOP*alpha, lambda);