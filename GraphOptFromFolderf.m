function [ScrN,BatchN,MinN,LastFN]=GraphOptFromFolderf(BaseGP,ShowFig)
if(nargin<2)
    ShowFig=false;
end
if(BaseGP(end)~=filesep)
    BaseGP=[BaseGP filesep];
end
D=dir([BaseGP '*.png']);
Datenums=[D.datenum];
DNames={D.name};
[SDatenums, Ord]=sort(Datenums);
SDNames=DNames(Ord);
SDNames=SDNames(4:end);
SDatenums=SDatenums(4:end);
LastFN=SDNames{end};

X=cat(1,SDNames{:});
X=textscan(X.','batch%06d_out_%f.png');
BatchN=X{1};
ScrN=X{2};
MinN=(SDatenums-SDatenums(1))*86400/60;
if(ShowFig)
    figure;plot(MinN(3:end),ScrN(3:end));
    ylabel('Score');
    xlabel('Minutes');
end