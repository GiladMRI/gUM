W=grmss(RecByPart(:,:,2:3),3);
W(:,70:end)=0;
Trg=exp(1i*angle(DD(:,:,2)));
x=[0.1 0.5 0.3];
x=[0 -[0.1 0.25]*5];
SW=SmoothBySlices(W.*Trg,[20 20],5);
Trgx=exp(1i.*angle(SW));

[X, Y]=ndgrid(linspace(-0.5,0.5,Sz(1)),linspace(-0.5,0.5,Sz(2)));

LinFieldFromParams=@(x) exp(1i.*(x(1)+ X*x(2)+Y*x(3)));

CostFunc=@(x) grmss(W.*Trgx-W.*LinFieldFromParams(x));

[BestX, BestVal]=fminsearch(CostFunc,x);

ShowAbsAngle(W.*LinFieldFromParams(BestX))

ShowAbsAngle(W.*Trg)

ShowAbsAngle(W.*(Trg./LinFieldFromParams(BestX)))


angle(gsum(W.*(Trg./LinFieldFromParams(BestX))))

fgmontagex(angle(Trg./LinFieldFromParams(BestX)),[-pi pi]/10);colormap jet


fgmontagex(angle(Trg),[-pi pi]/10);colormap jet
fgmontagex(angle(LinFieldFromParams(x)),[-pi pi]/10);colormap jet
fgmontagex(angle(LinFieldFromParams(BestX)),[-pi pi]/10);colormap jet