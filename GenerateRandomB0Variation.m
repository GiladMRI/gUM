function Out=GenerateRandomB0Variation(N,LFac,QFac,SFac,RandV)
if(numel(N)==1)
    N=[N N];
end
if(nargin<5)
    RandV=rand(13,1);
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

AL=(RandV(1)-0.5)*LFac*X+(RandV(2)-0.5)*LFac*Y+(RandV(3)-0.5)*QFac*(X.^2)+(RandV(4)-0.5)*QFac*(Y.^2)+(RandV(5)-0.5)*QFac*(X.*Y);
BL=(RandV(6)-0.5)*LFac*X+(RandV(7)-0.5)*LFac*Y+(RandV(8)-0.5)*QFac*(X.^2)+(RandV(9)-0.5)*QFac*(Y.^2)+(RandV(10)-0.5)*QFac*(X.*Y);

PX=(RandV(11)-0.5)*sin(AL)+(RandV(12)-0.5)*sin(BL);

% DCPhase=RandV(11)*2*pi-pi;
% PX=PX*2*pi*SFac+DCPhase;
PX=(PX+RandV(13)-0.5)*SFac;

Out=PX;