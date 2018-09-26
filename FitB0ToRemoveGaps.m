function X2=FitB0ToRemoveGaps(ToFitSli,WSli,DForFit)
% ToFit=B0M2;
W=abs(WSli);
[x,y]=meshgrid(1:size(ToFitSli,2),1:size(ToFitSli,1));
xF=x(:);
yF=y(:);
wF=W(:);
zF=ToFitSli(:);

% DForFit=5;
% 5: 1.14 + 4.19 sec
x=xF(1:DForFit:end);
y=yF(1:DForFit:end);
z=zF(1:DForFit:end);
w=wF(1:DForFit:end);
fo = fitoptions('Weights',w,'Method','LowessFit','Span',0.01);
tic
sf = fit([x, y],double(z),'lowess',fo);
% sf = fit([x, y],z,fo)
toc
% figure;plot(sf,[x,y],z)
tic
X=sf([xF(1:1:end),yF(1:1:end)]);
toc
X2=reshape(X,size(ToFitSli));

% fgmontage(ToFitSli,[-800 400])
% fgmontage(X2,[-800 400])