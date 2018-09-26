function Out=gBARTnuFT(Traj,Data,Sz)

if(mod(Sz(1),2)==1)
    Traj(1,:)=Traj(1,:)-0.5;
end
if(mod(Sz(2),2)==1)
    Traj(2,:)=Traj(2,:)-0.5;
end

Out=bart('nufft ',Traj,Data);

if(mod(Sz(1),2)==1)
    Out=Out.*exp(1i* 2*( (Traj(1,:)+0.5)*pi/Sz(1)));
end
if(mod(Sz(2),2)==1)
    Out=Out.*exp(1i* 2*( (Traj(2,:)+0.5)*pi/Sz(2)));
end
