function Cg=mintimegrad_radm(Nest,g0_mTm,gf_mTm,deltakMoment_radm,DwellTimeGrad_ms,Gmax_mTm, Smax_Tms)
if(~isreal(g0_mTm))
    g0_mTm=[real(g0_mTm) imag(g0_mTm)];
end
if(~isreal(gf_mTm))
    gf_mTm=[real(gf_mTm) imag(gf_mTm)];
end
if(~isreal(deltakMoment_radm))
    deltakMoment_radm=[real(deltakMoment_radm) imag(deltakMoment_radm)];
end

mTmToGcm=0.1;
TmsToGcms=100;
radmTo1cm=1/100/2/pi;
g0=g0_mTm*mTmToGcm;
gf=gf_mTm*mTmToGcm;
deltakMoment=deltakMoment_radm*radmTo1cm; % /cm, s/cm, s^2/cm
T		= DwellTimeGrad_ms/1000;  % Sampling time (s).
gmax=Gmax_mTm*mTmToGcm;
smax=Smax_Tms*TmsToGcms;
t0	=0;
type=3;

[g,v] = mintimegrad(Nest,g0,gf,deltakMoment,T,gmax,smax,t0,type);
Cg=(g(:,1)+1i*g(:,2))/mTmToGcm;