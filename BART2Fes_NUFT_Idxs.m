function Out=BART2Fes_NUFT_Idxs(In,Sz)
for i=1:numel(Sz)
    Out(:,i)=(In(i,:).'*2*pi/Sz(i));
end