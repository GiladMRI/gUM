function Out=FES2BART_NUFT_Idxs(In,Sz)
for i=1:numel(Sz)
    Out(i,:)=(In(:,i).'*Sz(i)/(2*pi));
end
for i=(numel(Sz)+1):3
    Out(i,:)=Out(1,:)*0;
end