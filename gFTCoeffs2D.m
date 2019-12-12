function Out=gFTCoeffs2D(Sz2,Idxs2)

nPoints=size(Idxs2,1);
[U1, IA1, IB1]=unique(Idxs2(:,1));
[U2, IA2, IB2]=unique(Idxs2(:,2));

Coeffs1=gFTCoeffs1D(Sz2(1),U1);
Coeffs2=gFTCoeffs1D(Sz2(2),U2);

Coeffs1=Coeffs1(IB1,:);
Coeffs2=Coeffs2(IB2,:);

Out=NaN([Sz2 nPoints]);
for i=1:nPoints
    Out(:,:,i)=Coeffs1(i,:).'*Coeffs2(i,:);
end
Out=Out./sqrt(prod(Sz2));