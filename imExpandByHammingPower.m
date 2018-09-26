function Out=imExpandByHammingPower(In,TrgSz,Pow)
if(nargin<3)
    Pow=5;
end
if(numel(TrgSz)==1)
    TrgSz(2)=size(In,2);
end
Pad=(TrgSz-gsize(In,1:2))/2;
Wnd1=hamming(size(In,1)).^Pow;
Wnd2=(hamming(size(In,2)).^Pow).';
Out=In;
if(Pad(1)>0)
    Out=fft1cg(Out,1);
    Out=RepDotMult(Out,Wnd1);
end
if(Pad(2)>0)
    Out=fft1cg(Out,2);
    Out=RepDotMult(Out,Wnd2);
end

Out=padBoth(Out,Pad);

if(Pad(1)>0)
    Out=ifft1cg(Out,1);
end
if(Pad(2)>0)
    Out=ifft1cg(Out,2);
end