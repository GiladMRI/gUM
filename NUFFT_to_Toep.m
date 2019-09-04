function fftkernBigc=NUFFT_to_Toep(nufftStruct,wx)
Sz128=nufftStruct.Nd;
N1=Sz128(1);
N2=Sz128(2);

SN=nufftStruct.sn;
try
    P=nufftStruct.p.G;
catch
    P=nufftStruct.p;
end

TakeTopLeftBlock=@(x,Sz) x(1:Sz(1),1:Sz(2),:,:,:,:);

fft2osN=@(x) fft2(x,size(x,1)*2,size(x,2)*2)/sqrt(prod(gsize(x,1:2)*2));
ifft2osN=@(x) TakeTopLeftBlock(ifft2(x),[size(x,1)/2,size(x,2)/2])*sqrt(prod(gsize(x,1:2)));

forw=@(x) P*reshape(fft2osN(x.*SN),[],1);
back=@(x) ifft2osN(reshape(P'*x,[N1*2,N2*2])).*conj(SN);

TrgSz=[Sz128*2, gsize(wx,3:20)];
wx=reshape(wx,size(wx,1),[]);

v11=zeros(Sz128);
v12=zeros(Sz128);
v21=zeros(Sz128);
v22=zeros(Sz128);
v11(1,1)=1;
v12(1,end)=1;
v21(end,1)=1;
v22(end,end)=1;
f11=forw(v11);
f12=forw(v12);
f21=forw(v21);
f22=forw(v22);
fftkernBigc=zeros([Sz128*2, size(wx,2)]);
for i=1:size(wx,2)
    block11=back(f11.*wx(:,i));
    block12=back(f12.*wx(:,i));
    block21=back(f21.*wx(:,i));
    block22=back(f22.*wx(:,i));
    
    Big=zeros(Sz128*2);
    Big(1:N1,1:N2)=block22;
    Big(N1:end-1,N2:end-1)=block11;
    Big(1:N1,N2:end-1)=block21;
    Big(N1:end-1,1:N2)=block12;
    
    Bigc=circshift(Big,[-N1+1,-N2+1]);
    fftkernBigc(:,:,i) = fft2(Bigc);
end
fftkernBigc=reshape(fftkernBigc,TrgSz);
%%
% Sig=nufft(v11,nufftStruct);
% Sig=Sig.*wx;
% block11=nufft_adj(Sig, nufftStruct);
% Sig=nufft(v12,nufftStruct);
% Sig=Sig.*wx;
% block12=nufft_adj(Sig, nufftStruct);
% Sig=nufft(v21,nufftStruct);
% Sig=Sig.*wx;
% block21=nufft_adj(Sig, nufftStruct);
% Sig=nufft(v22,nufftStruct);
% Sig=Sig.*wx;
% block22=nufft_adj(Sig, nufftStruct);
%%
% SN=nufftStruct.sn;
% P=nufftStruct.p.G;
% forw=@(x) P*reshape(fft2(padarray(x.*SN,[N1,N2],'post')),[],1);
% back=@(x) subsref(ifft2(reshape(P'*x,[N1*2,N2*2])),struct('type','()','subs',{{1:N1,1:N2}})).*conj(SN);
% block11A=back(forw(v11).*wx*N1*2*N2*2);
% block12A=back(forw(v12).*wx*N1*2*N2*2);
% block21A=back(forw(v21).*wx*N1*2*N2*2);
% block22A=back(forw(v22).*wx*N1*2*N2*2);
% 
% block11A(end,end)
% conj(block22A(1,1))