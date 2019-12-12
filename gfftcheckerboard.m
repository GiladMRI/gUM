function Out=gfftcheckerboard(Sz)
Out=(checkerboard(1,Sz(1)/2,Sz(2)/2)<0.5)*2-1;