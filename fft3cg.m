function res = fft3cg(x)

% res = fft2c(x)
% 
% orthonormal forward 3D FFT
%
% (c) Michael Lustig 2005

% res = 1/sqrt(prod(gsize(x,1:2)))*fftshift(fft2(ifftshift(x)));
res = 1/sqrt(prod(gsize(x,1:3)))*gfftshift(gfft(gifftshift(x,1:3),1:3),1:3);

