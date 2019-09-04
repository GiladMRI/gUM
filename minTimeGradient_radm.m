function [k,g,s]=minTimeGradient_radm(InTraj,RV, g0, gfin, Gmax_mTm, Smax_Tms,DwellTimeGrad_ms)
% InTraj=kTrajQ;
if(size(InTraj,2)==2 || size(InTraj,2)==3)
    InTraj=InTraj.';
end
if(~isreal(InTraj))
    InTraj=[real(InTraj(:)), imag(InTraj(:))];
end
if(isempty(RV))
    RV=0;
end
if(isempty(DwellTimeGrad_ms))
    DwellTimeGrad_ms=10e-3;
end
if(isempty(g0))
    g0=0;
end
% if(isempty(gfin))
%     gfin=-1;
% end
if(~isreal(gfin))
    gfin=[real(gfin), imag(gfin)];
end

Gcm2mTm=10;
GCmms2Tms=10;
radm2cm=1/(2*pi*100);

Gmax_GCcm=Gmax_mTm/Gcm2mTm;
Smax_GCmms=Smax_Tms/GCmms2Tms;

Traj1cm=InTraj*radm2cm;
which minTimeGradient
if(strhas(which('minTimeGradient'),'mex'))
    error('minTimeGradient_radm : Dont use the mex version');
end
% [C,time_ms,g,s,k, phi, sta, stb] = minTimeGradient(Traj1cm,RV, g0, gfin, Gmax_GCcm, Smax_GCmms,DwellTimeGrad_ms);
[C,time_ms,g,s,k, phi, sta, stb] = minTimeGradient([Traj1cm Traj1cm(:,1)*0],RV, g0, gfin, Gmax_GCcm, Smax_GCmms,DwellTimeGrad_ms);

k=k/radm2cm;
g=g*Gcm2mTm;
s=s*GCmms2Tms;

k=k(:,1:2);
g=g(:,1:2);
s=s(:,1:2);

k=k(:,1)+1i*k(:,2);
g=g(:,1)+1i*g(:,2);
s=s(:,1)+1i*s(:,2);