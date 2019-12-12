function Out=getB0_ExpPolyfit(In,dT_ms,pOrder);
% In=tmp(:,:,WhichTSToUset,:);
% dT_ms=dT_Med_ms;
% %%
dispstat('','init');

Sz=gsize(In,1:2);
HankelTemporalLen=2;
nEchoes=size(In,3);
[~, ~, ~,H]=ghankel(nEchoes,HankelTemporalLen,Sz);

nSlices=size(In,4);
Out=zeros([Sz nSlices]);

[X,Y]=ndgrid(linspace(-1,1,Sz(1)),linspace(-1,1,Sz(2)));
switch pOrder
    case 0
        RegsM=ones(Sz);
    case 1
        RegsM=cat(3,ones(Sz),X,Y);
    case 2
        RegsM=cat(3,ones(Sz),X,Y,X.^2,X.*Y,Y.^2);
    case 3
        RegsM=cat(3,ones(Sz),X,Y,X.^2,X.*Y,Y.^2,X.^3,(X.^2).*Y,(Y.^2).*X,Y.^3);
    otherwise
        error('wrong pOrder');
end 
x0=zeros(1,size(RegsM,3));
FldFunc=@(x) sum(RegsM.*perm32(x),3);

for i=1:nSlices
    dispstat(['Slice ' num2str(i)],'timestamp')
    HIn=H*In(:,:,:,i);
    
    CostFunc=@(x) grmss(HIn(:,:,:,2)-  HIn(:,:,:,1).*exp(-1i.*FldFunc(x)));
    
    [BestX, BestCost]=fminsearch(CostFunc,x0);
    Out(:,:,i)=FldFunc(BestX)*500/pi/dT_ms;
end