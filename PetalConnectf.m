function Out=PetalConnectf(r1,r2,ArcPhiStart,dPoint)
% r1=Ks*KRadAlpha;
% r2=Ks;

% xy = [0 0; 0 r*0.9; 0 r; 7 12; 17 5].';

% xy = [0 0; 0 5; 0 10; 2 14; 7 12; 17 5].';
% P2=r2.*exp(1i*pi/5);
P2=r2.*exp(1i*ArcPhiStart);

% PerpVec=[1 -imag(P2)/real(P2)];
PerpVecC=P2*exp(1i*2*pi/4);
P3=P2+PerpVecC*0.5;

xy = [0 0; 0 r1; imag(P2)*0.18 r1+ (real(P2)-r1)*1.14; imag(P2) real(P2); imag(P3) real(P3)].';
npts=size(xy,2);
% plot(xy(1,:),xy(2,:),'ro','LineWidth',2);
% text(xy(1,:), xy(2,:),[repmat('  ',npts,1), num2str((1:npts)')])
% ax = gca;
% ax.XTick = [];
% ax.YTick = [];

Curve=cscvn(xy);
% hold on
% fnplt(cscvn(xy),'r')
% axis square
% axis equal 
% hold off
%
% Ixy=fnval(Curve,Curve.breaks(2):0.01:Curve.breaks(4));
Ixy=fnval(Curve,linspace(Curve.breaks(2),Curve.breaks(4),100));
IxyC=Ixy(1,:)+1i*Ixy(2,:);
L=sum(abs(diff(IxyC)));
Ixy=fnval(Curve,linspace(Curve.breaks(2),Curve.breaks(4),ceil(L/dPoint)));
% figure;
Out=Ixy(1,:)+1i*Ixy(2,:);
% plot(Ixy(1,:),Ixy(2,:),'b')



