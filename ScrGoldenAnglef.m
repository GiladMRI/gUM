function Out=ScrGoldenAnglef(dPhi)

Ns=3:70;
for N=Ns
    C=exp(1i*dPhi*(0:N-1));
%     C=[C -C];
    Phis=angle(C);
    S1=sort(Phis+pi);
    S2=[S1 S1(1)+2*pi];
    dPhis=diff(S2);
    Scr(N)=max(dPhis);
end
% [mean(Scr(3:30)) max(Scr(3:30))]
Out=max(Scr(Ns));
% Out=mean(Scr(Ns));