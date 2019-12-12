function Out=gTMat2D3to4(In)
Out([1 2 4],[1 2 4],:,:,:,:)=In;
Out(3,3,:)=1;