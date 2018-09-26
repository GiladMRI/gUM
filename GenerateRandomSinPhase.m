function Out=GenerateRandomSinPhase(N,LFac,QFac,RandV)
if(numel(N)==1)
    N=[N N];
end
if(nargin<4)
    RandV=rand(11,1);
end
if(nargin<3)
    QFac=1;
end
if(nargin<2)
    LFac=5;
end
Linx=linspace(-pi,pi,N(1));
Liny=linspace(-pi,pi,N(2));

[X,Y]=ndgrid(Linx,Liny);

AL=(RandV(1)-0.5)*LFac*X+(RandV(2)-0.5)*LFac*Y+(RandV(3)-0.5)*QFac*(X.^2)+(RandV(4)-0.5)*QFac*(Y.^2);
BL=(RandV(5)-0.5)*LFac*X+(RandV(6)-0.5)*LFac*Y+(RandV(7)-0.5)*QFac*(X.^2)+(RandV(8)-0.5)*QFac*(Y.^2);

PX=(RandV(9)-0.5)*sin(AL)+(RandV(10)-0.5)*sin(BL);

DCPhase=RandV(11)*2*pi-pi;
PX=PX*pi+DCPhase;

Out=exp(1i*PX);