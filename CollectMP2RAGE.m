BP='/media/a/DATA/MP2RAGE_WithPhase/';
D=dir(BP);
D=D([D.isdir]);
D=D(3:end);
for i=1:numel(D)
    D2=dir([BP D(i).name filesep '*.nii']);
    % numel(D2)
    M=loadniidata([BP D(i).name filesep D2(1).name]);
    P=loadniidata([BP D(i).name filesep D2(2).name])*2*pi/4095;
    C(:,:,:,:,i)=M.*exp(1i*P);
end