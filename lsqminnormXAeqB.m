function Out=lsqminnormXAeqB(A,B)
% X = lsqminnorm(A,B) returns an array X that solves the linear equation AX = B an
% XA=B
% so A'X'=B'
Out=lsqminnorm(A',B')';