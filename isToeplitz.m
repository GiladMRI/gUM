function Out=isToeplitz(M)
Out=grmss(M-asToeplitz(M))/grmss(M);