function res = ifft3cg(x)

res = sqrt(prod(gsize(x,1:3)))*gfftshift(ifft(ifft2(gifftshift(x,1:3)),[],3),1:3);