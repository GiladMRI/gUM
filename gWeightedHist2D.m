function WH=gWeightedHist2D(In,Bins,W)
% BinCenters=(Bins(1:end-1)+Bins(2:end))/2;
[Ha,~,bin] = histcounts(In(:),Bins);
WH=Ha*0;
for i=1:numel(Ha)
    WH(i)=sum(W(bin==i));
end
% SWHT2s=max(WHT2s,sum(WHT2s)*0.03/numel(WHT2s));

% clear WDecays
% for i=2:numel(HT2s)-1
%     WDecays(i,:)=SWHT2s(i)*exp(-TimePointsMed_ms./T2sCenters(i-1));
% end