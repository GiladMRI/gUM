EchoSpacing_ms=0.93;
BaseP='/autofs/cluster/kawin/Gilad/50EchoData/';
load([BaseP 'Full3DAllEchos.mat']);

load('/autofs/cluster/kawin/Gilad/50EchoData/meas_prot_Fully_Sample.mat');
%%
nEchos=7;
Combined=Full50EchoS(:,:,5:10:end,1:nEchos);
nSlices=size(Combined,3);
dTEs_ms=parameters.iEffectiveEpiEchoSpacing/1000;
TEs_ms=2;
%%
WhichEchosToUse=1:nEchos;
for i=1:nSlices
    disp(i);
    [PDBase(:,:,i), UpdatedB0Map_Hz(:,:,i), UpdatedT2SMap_ms(:,:,i), s_vals(:,:,:,i), Fitted0(:,:,:,i), PDBase0(:,:,i)]=FitToModel_MPBD1CSf(Combined(:,:,i,:),WhichEchosToUse,dTEs_ms(1),TEs_ms(1));
end
%%
fgmontagex(abs(UpdatedT2SMap_ms(:,:,3:8:end)),[0 100]);title('GRE-ME T_2^*');
fgmontagex(UpdatedB0Map_Hz(:,:,3:8:end),[-300 300]);title(['GRE-ME B_0, TEs: ' num2str(TEs_ms,'%.1f ')]);
ShowAbsAngle(PDBase0(:,:,3:8:end),1,'Size',[2 2])
fgmontagex(s_vals(:,:,:,3:8:end));title('GRE-ME SV maps');
%%
PDBase0x=min(PDBase0,6*grmss(PDBase0));
[Out B1 BN1]=CalcSlicesSNR(abs(PDBase0x(:,:,:)),false,5);
B2=(~BN1).*SensMsk;
B2D=imdilate(B2,strel('disk',3,8));
B3=imfillholesBySlices( B2D );
for i=1:nSlices
    B4(:,:,i)=getLargestComponent(B3(:,:,i));
end
B4=B4.*SensMsk;
%%
dAngle=UpdatedB0Map_Hz*2*pi*dTEs_ms(1)/1000;
[unwrapped] = cusackUnwrap(dAngle, grmss(Combined,4));
unwrapped=unwrapped.*B4;
B0M_Hz=unwrapped*1000/2/pi/dTEs_ms(1);
fgmontagex(B0M_Hz(:,:,1:16),[-300 300]);title(['GRE-ME B_0 unwrapped, TEs: ' num2str(TEs_ms,'%.1f ')]);
fgmontagex(UpdatedB0Map_Hz(:,:,1:16),[-300 300]);title(['GRE-ME B_0, TEs: ' num2str(TEs_ms,'%.1f ')]);
