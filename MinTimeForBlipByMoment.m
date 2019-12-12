% MinTimeForBlipByMoment(MaxAmp,MaxSlew)
MaxAmp=35; % mT/m
MaxSlew=155; % mT/m/ms
gamma=42.5774806; % MHz/T
Separation=36; % mm 
DesiredMoment=2/(Separation*gamma);
% gamma*G*Separation*T=.5; % MHz/T*mT/m*mm*us
% Triangle:
% Area is T*H/2, or HalfT*H
% H is HalfT*MaxSlew
% Area is HalfT*^2*MaxSlew
HalfT=1000*sqrt(DesiredMoment/(MaxSlew));
HalfTRounded=ceil(HalfT/10)*10;
nSteps=HalfTRounded/10;
FHalfnSteps=floor(nSteps/2);
if(mod(nSteps,2)==0)
    nBlocks=FHalfnSteps*(FHalfnSteps+1);
else
    
    nBlocks=FHalfnSteps*(FHalfnSteps+1)+FHalfnSteps+1;
end
BlockMoment=DesiredMoment/nBlocks;
BlockH=BlockMoment/10;
for i=1:FHalfnSteps
    Grad(i)=BlockH*i;
    Grad(nSteps+1-i)=BlockH*i;
end
if(mod(nSteps,2)==1)
    Grad(FHalfnSteps+1)=BlockH*(FHalfnSteps+1);
end