function [Out Idxs]=gmin(in,dims,opts,Msk)

if(nargin>3)
    in(~Msk)=Inf;
end

if(~exist('dims','var'))
    % Get the index and value of the minimum value
    [Out a]=max(in(:));
    % If we return more than 1 argument it means we want the index as well
    if(nargout>1)
        % Get the subscript index out of the linear one
        [Idxs{1:ndims(in)}]=ind2sub(size(in),a);
        Idxs=[Idxs{:}];
    end
    return;
end

Out=max(in,[],dims(1));

for i=2:length(dims)
    Out=min(Out,[],dims(i));
end

if(~exist('opts','var'))
    return;
end

if(ismember(1,opts))
    Out=squeeze(Out);
end