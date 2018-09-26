function [TSB,TSC,nTSX,CondGHG]=GetTSCoeffsByMinMaxHist(B0MapToUse,TimeInMs2,MaxBoHz,MaxCond)
if(nargin<3)
    MaxBoHz=500;
end
if(nargin<4)
    MaxCond=100;
end
%% now use histogram
HistEdges=linspaceWithHalfStep(-MaxBoHz,MaxBoHz,MaxBoHz*2);
BinCenters=(HistEdges(1:end-1)+HistEdges(2:end))/2;
B0Hist=histcounts(B0MapToUse(:),HistEdges);
B0HistC=B0Hist.';
%%
TimeIns=TimeInMs2/1000;

clear GHGC
for nTSX=4:20
    tau_s=(TimeInMs2(end)/1000)/(nTSX-1);
    for i=0:nTSX-1
        for j=0:nTSX-1
            dT=(i-j)*tau_s;
            GHGC{nTSX}(i+1,j+1)=exp(-1i*2*pi*BinCenters*dT)*B0HistC;
        end
    end
    CondGHG(nTSX)=cond(GHGC{nTSX});
    if(CondGHG(nTSX)>MaxCond)
        disp(['Reached cond ' num2str(CondGHG(nTSX),'%.2g') ' for #TS=' num2str(nTSX)]);
        break;
    end
end
nTSX=nTSX-1;
disp(['Using #TS=' num2str(nTSX) ' with cond ' num2str(CondGHG(nTSX),'%.2g')]);
GHG=GHGC{nTSX};
% IGHG=inv(GHG);

tau_s=(TimeInMs2(end)/1000)/(nTSX-1);
clear GHb
for i=0:nTSX-1
    GHb(:,i+1)=exp(1i*2*pi*(TimeIns.'-tau_s*i)*BinCenters)*B0HistC;
end
% TSB=IGHG*(Gsb.');
TSB=GHG\(GHb.');
TSB=TSB.';
FesTimePoints=linspace(TimeInMs2(1),TimeInMs2(end)/1000,nTSX);
TSC=exp(1i*2*pi*RepDotMult(B0MapToUse,gpermute(FesTimePoints,[3 2])));  % exp(1i*2*pi*(TimeInMs2/1000)*B0M2);

disp('GetTSCoeffsByMinMaxHist end')