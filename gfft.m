function X=gfft(X,dims)
for i=1:numel(dims)
    X=fft(X,[],dims(i));
end