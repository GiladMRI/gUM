% rmpath('/media/g/aa384808-2266-4edb-87a7-637bc772dc2f/SPENOnline/CS/ESPIRiT/nufft_files')
ig = image_geom('nx', 32, 'ny', 28, 'dx', 1, 'offsets', 'dsp');
ig.mask = ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
N = ig.dim;

[kspace, omega, wi] = mri_trajectory('spiral0', {}, N, ig.fov, {'voronoi'});
warn 'todo: crudely fixing dcf at edge of spiral'
wi(350:end) = wi(350);

% nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
nufft_args = {N, [6 6], 2*N, N/2, 'minmax:kb'};

% Gm = Gmri(kspace, ig.mask, 'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);
Gm = Gmri(kspace, ig.mask, 'fov', ig.fov, 'basis', {'dirac'}, 'nufft', nufft_args);

Tm = build_gram(Gm, 1);

v=rand(N);
Tmv=Tm*v;


x0 = v;
x0m=x0(ig.mask);
A=Gm*x0m;
A2=Gm.Gnufft*x0m;
x1 = ig.embed(Gm' * (Gm * x0(ig.mask)));
x2 = ig.embed(Tm * x0(ig.mask));
% im clf, im(stackup(x1,x2))

figure;plot(real(x1(:)),'r.');hold on;plot(real(x1(:)),'go');