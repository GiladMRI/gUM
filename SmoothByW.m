function Out=SmoothByW(In,W,SRange,SSigs)
if(nargin<3)
    SRange=[21 21];
end
if(nargin<4)
    SSigs=[0.0001 0.1:0.1:29];
end 
M0B0=W.*In;
Sz=gsize(In,1:2);
nSlces=size(In,3);
Out=zeros([Sz nSlces]);
SM0B0=zeros([Sz numel(SSigs)]);
SM0=zeros([Sz numel(SSigs)]);
for s=1:nSlces
    disp(s);
    CurW=W(:,:,s);
    for i=1:numel(SSigs)
        SM0B0(:,:,i)=SmoothBySlices(M0B0(:,:,s),SRange,SSigs(i));
        SM0(:,:,i)=SmoothBySlices(CurW,SRange,SSigs(i));
    end
    mThresh=median(CurW(CurW>0));
    MskSM0=(SM0.*perm32(SSigs))>mThresh;
    MskSM0(:,:,end)=1;
    FirstTimeOverT=numel(SSigs)+1-sum(MskSM0,3);
    SB0=SM0B0./SM0;
    clear SB0x
    SB0x=zeros(gsize(SB0,1:2));
    for i=1:size(SB0,1)
        for j=1:size(SB0,2)
            SB0x(i,j)=SB0(i,j,FirstTimeOverT(i,j));
        end
    end
    SB0x(isnan(SB0x))=0;
    Out(:,:,s)=SB0x;
end