t=linspace(0,1,10000);
S=sin(2*pi*t*11.25-pi/2);
Big=repmat(S,[100 1]);
Sm=imresize(Big,[100 100]);

fgmontage(Sm.');
%%
B0Here=gflip(B0RealEx.',2);
B0Tmp=imresize(B0Here,[100 100]);
B0TmpB=imresize(B0Here,[100 10000]);

i=50;

Fac=4;
NewLoc=t*size(B0TmpB,2)+B0TmpB(i,:)*Fac;
figure;
subplot(4,1,1);plot(t,NewLoc);
subplot(4,1,2);plot(t,B0TmpB(i,:)*Fac);
subplot(4,1,3);plot(NewLoc,S);setXaxis([0 10000]);
[N,edges,bin]=histcounts(NewLoc,linspace(0,10000,100));

for j=1:100
    Out(j)=sum(S(bin==j));
end
subplot(4,1,4);plot(Out);
%%
clear OutM
for i=1:100
    NewLoc=t*size(B0TmpB,2)+B0TmpB(i,:)*Fac;
    [N,edges,bin]=histcounts(NewLoc,linspace(0,10000,100));
    for j=1:100
        OutM(i,j)=sum(S(bin==j));
    end
end