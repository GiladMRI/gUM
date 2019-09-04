dPhi=111.25*2*pi/360;
dPhi=200*2*pi/360;

N=9;
for N=3:30
    C=exp(1i*dPhi*(0:N-1));
    Phis=angle(C);
    S1=sort(Phis+pi);
    S2=[S1 S1(1)+2*pi];
    dPhis=diff(S2);
    Scr(N)=max(dPhis);
end
[mean(Scr(3:30)) max(Scr(3:30))]
%%
Func=@(x) ScrGoldenAnglef(x);
[BestPhi BestScr]=fminsearch(Func,2);
[BestPhi*360/2/pi BestScr]