nTrajX=3000;
nLoopsX=12;
R=linspace(0,1,nTrajX)*96;
% R=R+real(exp(1i*linspace(0,2*pi*nLoopsX*19.3125432,nTrajX))).*(96/12/3);
Phi=linspace(0,2*pi*nLoopsX,nTrajX);
C=R.*exp(1i*Phi);
CS=[0 cumsum(abs(diff(C)))];
CS=CS./CS(end);
R=R+real(exp(1i*(CS*2*pi*nLoopsX*19.46624))).*(96/12/3);
C=R.*exp(1i*Phi);
figure;
plot(real(C),imag(C));
axis equal
setXaxis([-100 100]);
setYaxis([-100 100]);
removeTicks