function [ alpha_thresh ] = TV_thresh(alpha, lambda)
% function [ alpha_thresh ] = TV_thresh(alpha, lambda)

alpha_thresh = TVOP'*SoftThresh(TVOP*alpha, lambda);