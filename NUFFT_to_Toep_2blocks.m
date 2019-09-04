function fftkern=NUFFT_to_Toep_2blocks(nufftStruct,wi)

if(nargin<2)
    wi=ones(nufftStruct.M,1);
end

Sz=nufftStruct.Nd;
N1=Sz(1);
N2=Sz(2);

wi=double(wi);

Trajm2m=nufftStruct.om;
Trajm2m(:,1)=-Trajm2m(:,1);

nufftStruct2 = nufft_init(Trajm2m, Sz, nufftStruct.Jd, Sz*2, nufftStruct.n_shift); % , [0 0] st.om
z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
n1 = 0:(N1*2-1);
n2 = 0:(N2*2-1);
[n1, n2] = ndgrid(n1, n2);
m1 = mod(-n1, N1*2); % circular symmetry
m2 = mod(-n2, N2*2);

kernA=zeros([N1*2,N2*2,gsize(wi,2:20)]);
% fftkern=zeros([N1*2,N2*2,gsize(wi,2:20)]);

block1A=nufft_adj(wi(:,:), nufftStruct);
block2A=nufft_adj(wi(:,:), nufftStruct2);
    
for i=1:prod(gsize(wi,2:20))
%     block1=nufft_adj(wi(:,i), nufftStruct);
%     block2=nufft_adj(wi(:,i), nufftStruct2);
    block1=block1A(:,:,i);
    block2=block2A(:,:,i);
    
    kern = [
        [block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
        zeros(1,2*N2);
        [flipud(block2(2:N1,:)) z2 rot90(conj(block1(2:N1, 2:N2)),2)]
        ]; % [(2Nd)]
    tmp1 = kern(sub2ind([N1*2 N2*2], 1+n1, 1+n2));
    tmp2 = kern(sub2ind([N1*2 N2*2], 1+m1, 1+m2));
    tmp2 = conj(tmp2);
%     kern = (tmp1+tmp2)/2; % force it to be Hermitian
%     kern = kern/(prod(Sz)^2*4);
    kernA(:,:,i) = (tmp1+tmp2)/2; % force it to be Hermitian
    kernA(:,:,i) = kernA(:,:,i)/(prod(Sz)^2*4);
%     fftkern(:,:,i) = fftn(kern);
end
fftkern = fft2(kernA);
fftkern = real(fftkern);
% fftkern = single(fftkern);
% disp('fftnern ok');