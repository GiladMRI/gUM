function Out=GenerateRandomBrokenPhase(N,LFac,QFac,nP)
if(numel(N)==1)
    N=[N N];
end
if(nargin<4)
    nP=2;
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
% 
% P=(rand-0.5)*5*X+(rand-0.5)*5*Y+(rand-0.5)*5*X.^2+(rand-0.5)*5*Y.^2;
% ShowAbsAngle(exp(1i*P))
% %%
% P1=((rand-0.5)*5*X+(rand-0.5)*5*Y)>(rand-0.5)*10;
% fgmontage(P1)
% %%
% LFac=5;
% QFac=1;
Base=(rand-0.5)*5*LFac*X+(rand-0.5)*5*LFac*Y+(rand-0.5)*5*QFac*X.^2+(rand-0.5)*5*QFac*Y.^2;
MskAll=false(size(Base));
PX=zeros(size(Base));
% nP=2;
for i=1:nP
    P=(rand-0.5)*LFac*X+(rand-0.5)*LFac*Y+(rand-0.5)*QFac*X.^2+(rand-0.5)*QFac*Y.^2;
    P1=((rand-0.5)*5*X+(rand-0.5)*5*Y)>(rand-0.5)*10;
    PX=PX+P.*P1;
    MskAll(P1)=true;
end
PX(~MskAll)=Base(~MskAll);
Out=exp(1i*PX);
% ShowAbsAngle(exp(1i*PX))