function C=gSpiral(N,nRounds,PowCoeff,phiBase)
% Simple VD Spiral
% Gilad

if(nargin<4)
    phiBase=0;
end
% N=500;
% nRounds=10;
% PowCoeff=2;

t=linspace(0,1,N);
r=(t.^PowCoeff);
phi=phiBase+2*pi*t*nRounds;

C=r.*exp(1i*phi);

% figure;plot(real(C),imag(C),'*-')

