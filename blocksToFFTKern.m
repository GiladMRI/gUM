function fftkern=blocksToFFTKern(block1,block2)

N1=size(block1,1);
N2=size(block1,2);
z1 = zeros(N1,1);
z2 = zeros(N1-1,1);
kern = [
    [block1 z1 conj(fliplr([block1(1,2:N2); block2(2:N1,2:N2)]))];
    zeros(1,2*N2);
    [flipud(block2(2:N1,:)) z2 fliplr(flipud(conj(block1(2:N1, 2:N2))))]
    ]; % [(2Nd)]
n1 = 0:(N1*2-1);
n2 = 0:(N2*2-1);
[n1, n2] = ndgrid(n1, n2);
tmp1 = kern(sub2ind([N1*2 N2*2], 1+n1, 1+n2));
m1 = mod(-n1, N1*2); % circular symmetry
m2 = mod(-n2, N2*2);
tmp2 = kern(sub2ind([N1*2 N2*2], 1+m1, 1+m2));
tmp2 = conj(tmp2);
kern = (tmp1+tmp2)/2; % force it to be Hermitian
fftkern = fftn(kern);
% fftkern = single(real(fftkern));
% disp('fftnern ok');