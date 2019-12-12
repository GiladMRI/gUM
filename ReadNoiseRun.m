BaseP='/media/a/DATA/NoiseRun/';
FNoise='meas_MID665_GRE_2mm_noise_ref0V_FID34346.dat';

sTwix = mapVBVD([BaseP FNoise],'removeOS','ignoreSeg');
Data=sTwix.image();

P=permute(Data,[1 3 4 2]);
R=reshape(P,[],32);
BadRows=all(R==0,2);
R=R(~BadRows,:);

RR=R'*R;

ShowAbsAngle(RR)
gprint(get(gcf,'Number'),[BaseP FNoise(1:end-4) '_Mat'],[])

figure;scatterhist(real(R(1:2:end,20)),real(R(1:2:end,21)))
gprint(get(gcf,'Number'),[BaseP FNoise(1:end-4) '_Ch20-21'],[])
save([BaseP FNoise(1:end-4) '_Noise.mat'],'R','RR');

%%
FN1='meas_MID661_GRE_2mm_FID34342.dat';
sTwix = mapVBVD([BaseP FN1],'removeOS','ignoreSeg');
Data1=sTwix.image();
Data1=squeeze(Data1);
I1=fft3cg(Data1);
figure(1);clf;gmontage(I1);
gprint(1,[BaseP FN1(1:end-4)],[])
%%
FN2='meas_MID667_GRE_1p25_iso_p2_FoV160_FID34348.dat';
sTwix = mapVBVD([BaseP FN2],'removeOS','ignoreSeg');
Data2=sTwix.image();
Data2=squeeze(Data2);
P2=permute(Data2,[1 3 4 2]);
I2=fft3cg(P2);

figure(1);clf;gmontage(AddSepLines(grmss(I2,4)))
gprint(1,[BaseP FN2(1:end-4)],[])
ShowAbsAngle(AddSepLines( squeeze(I2(:,:,60,:))))
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Ch32'],[])
%%
Ref=sTwix.refscan();
RefP=permute(Ref,[1 3 4 2]);
IRef=fft3cg(RefP);

figure(1);clf;gmontage(AddSepLines(grmss(IRef,4)),[0 15e-4])

figure(1);clf;gmontage(AddSepLines(squeeze(IRef(:,:,60,:))),[0 15e-4],'Size',[4 8])
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Ref_Ch32'],[])

%%
% [unfolded_data weights ] = doGrappaRecon_v2(Data2,sTwix.refscan )

%%
% [unfolded_data weights ] = doGrappaRecon_v2(varargin )
% function [unfolded_data weights ] = doGrappaRecon_v2(folded_data, ref_data, acc_fac, kernelsize )      if the kernel has not been fitted yet
% or 
% function [unfolded_data weights ] = doGrappaRecon_v2(folded_data, weights, acc_fac,  kernelsize )      if the kernel is already known
%
% this function performs a box standard GRAPPA reconstruction;  somewhat adapted form Mark Griswohld's OpenGrappa.m 
% INPUT: 
% folded_data    - k-space in format [COL LIN CHA .... SLC .... ]  SLC *must* be in DIM=10
% kernelsize     - size of grappa kernel [ x y ] (x odd and y even numbers)
% acc_fac         - acceleration factor
% ref_data       - ACS reference data as [COL LIN CHA SLC]
% weights        - grappa weights (OPTIONAL -- if they are not passed, they will be fitted) 
% OUTPUT:
% unfolded_data  - reconstructed data (LIN dimension will be acc_fac*LIN of the folded data)
% weights{slice} - Grappa weights for each slice for 
%%
FN2='meas_MID668_GRE_1p25_iso_p2_FoV160_ref0V_FID34349.dat';
sTwix = mapVBVD([BaseP FN2],'removeOS','ignoreSeg');
Data2=sTwix.image();
Data2=squeeze(Data2);
P2=permute(Data2,[1 3 4 2]);
R2=reshape(P2,[],32);

BadRows=all(R2==0,2);
R2=R2(~BadRows,:);

RR2=R2'*R2;
ShowAbsAngle(RR2)
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Mat'],[])

figure;scatterhist(real(R2(1:2:end,20)),real(R2(1:2:end,21)))
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Ch20-21'],[])
save([BaseP FN2(1:end-4) '_Noise.mat'],'R2','RR2');
%%
FN2='meas_MID446_GRE_2mm_noise_ref0V_FID44049.dat';
sTwix = mapVBVD([BaseP FN2],'removeOS','ignoreSeg');
Data2=sTwix.image();
Data2=squeeze(Data2);
P2=permute(Data2,[1 3 4 2]);
R3=reshape(P2,[],32);

BadRows=all(R3==0,2);
R3=R3(~BadRows,:);

for i=1:32
    HR(:,i)=histcounts(real(R3(:,i)),linspace(-1.5e-5,1.5e-5,100));
    HI(:,i)=histcounts(imag(R3(:,i)),linspace(-1.5e-5,1.5e-5,100));
end
%%
figure;
subplot(1,2,1);plot(HR);title('real')
subplot(1,2,2);plot(HI);title('imag')
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Hists'],[])
%%
figure;
Ns=32;
d=30;
dl=2e-5;
for i=1:Ns
    for j=(i+1):Ns
%         gsubplot(Ns,Ns,i,j);
%         scatter(real(R3(1:d:end,i)),real(R3(1:d:end,j)),'.')
%         plot(real(R3(1:d:end,i)),real(R3(1:d:end,j)),'.')
        plot(real(R3(1:d:end,i))+dl*i,real(R3(1:d:end,j))+dl*j,'b.')
        hold on
        
        plot(imag(R3(1:d:end,i))+dl*j,imag(R3(1:d:end,j))+dl*i,'b.')
        removeTicks
        
%         setXaxis([-1 1]*1e-5);
%         setYaxis([-1 1]*1e-5);
%         removeTicks;
        
%         gsubplot(Ns,Ns,j,i);
% %         scatter(imag(R3(1:d:end,i)),imag(R3(1:d:end,j)),'.')
%         plot(imag(R3(1:d:end,i)),imag(R3(1:d:end,j)),'.')
%         setXaxis([-1 1]*1e-5);
%         setYaxis([-1 1]*1e-5);
%         removeTicks;
    end
end
setXaxis([0 dl*Ns+2e-5]);
setYaxis([0 dl*Ns+2e-5]);
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_PlotMat'],[])
%%
RR3=R3'*R3;
ShowAbsAngle(RR3)
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Mat'],[])

figure;scatterhist(real(R3(1:2:end,20)),real(R3(1:2:end,21)))
gprint(get(gcf,'Number'),[BaseP FN2(1:end-4) '_Ch20-21'],[])
save([BaseP FN2(1:end-4) '_Noise.mat'],'R3','RR3');
