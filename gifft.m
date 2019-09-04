function X=gifft(X,dims)
for i=1:numel(dims)
    X=ifft(X,[],dims(i));
end