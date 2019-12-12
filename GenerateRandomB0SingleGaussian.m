function Out=GenerateRandomB0SingleGaussian(N,MaxB0,SigFac,MinSig)
if(nargin<2)
    MaxB0=200;
end
if(nargin<3)
    SigFac=1;
end
if(nargin<4)
    MinSig=0.05;
end
[X Y]=ndgrid(linspace(-0.5,0.5,N),linspace(-0.5,0.5,N));
Mu=rand(1,2)-0.5;
Sigx=rand*SigFac+MinSig;
Sigy=rand*SigFac+MinSig;
X=(X-Mu(1))/Sigx;
Y=(Y-Mu(2))/Sigy;
Amp=MaxB0*(rand*2-1);
Out=exp(- (X.^2+Y.^2)/2)*Amp;