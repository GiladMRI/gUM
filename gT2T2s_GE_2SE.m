function [S0,S02,T2,T2s,fitting_map_GRE,fitting_map_SE] = gT2T2s_GE_2SE(imgs,TEs_GRE,TEs_SE,te_SE,mask,threshold,np_SE1)

%  fit T2 and T2* 
% Inputs
% imgs:          input images at each echo time 
% TEs:            echo times 
% te2:            time of spin-echo 
% mask:          mask of voxels to fit 
% np_SE1:        time point index of the spin echo from the SE series
% Outputs
% So T2 T2s rmse : three parameter fit and rmse 
np_SE1=find(TEs_SE>te_SE,1)-1;

[nx,ny,~]=size(imgs);
np_GRE=size(TEs_GRE,1);
% np_SE1=(size(TEs_SE,1)+1)/2;
np_SE2=size(TEs_SE,1)-np_SE1;

A1 = [ones(np_GRE,1),zeros(np_GRE,1),-TEs_GRE(1:np_GRE),-TEs_GRE(1:np_GRE)];   %  ln(s0) ln (s02) R2' R2
A2 = [zeros(np_SE1,1), ones(np_SE1,1), TEs_SE(1:np_SE1)-te_SE, -TEs_SE(1:np_SE1)];   % ln(s0) ln (s02) R2' R2
A3 = [zeros(np_SE2,1), ones(np_SE2,1), te_SE-TEs_SE(np_SE1+1:end), -TEs_SE(np_SE1+1:end)];   % ln(s0) ln (s02) R2' R2
A = [A1;A2;A3];

S0=zeros(nx,ny);
S02=zeros(nx,ny);
T2=zeros(nx,ny);
T2s=zeros(nx,ny);
fitting_map_GRE=zeros(nx,ny,np_GRE);
fitting_map_SE=zeros(nx,ny,np_SE1+np_SE2);

for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=squeeze(imgs(ii,jj,:));
           maskk=abs(signal)>threshold;
           signal=signal(maskk);
           AA=A(repmat(maskk(:),[1 4]));
           AA=reshape(AA,[size(signal,1),4]);
           param = AA\log(signal);
%          param = lsqr(AA,log(signal),[],[],[],[],50*ones(4,1));
           %  result  
            S0(ii,jj) = exp(param(1));
            S02(ii,jj) = exp(param(2));
            T2(ii,jj) = 1/param(4);
            T2s(ii,jj) = 1/(param(3)+param(4));
            fitting_map_GRE(ii,jj,:)=S0(ii,jj)*exp(-TEs_GRE(:)*(param(3)+param(4)));
            fitting_map_SE(ii,jj,1:np_SE1)=S02(ii,jj)*exp(-te_SE*param(3)).*exp(-TEs_SE(1:np_SE1)*(param(4)-param(3)));
            fitting_map_SE(ii,jj,np_SE1+1:end)=S02(ii,jj)*exp(te_SE*param(3)).*exp(-TEs_SE(np_SE1+1:end)*(param(4)+param(3)));
        end
    end % end jj 
end % end ii

S02(isnan(S02)) = 0;
S0(isnan(S0)) = 0;
T2(isnan(T2)) = 0;
T2s(isnan(T2s)) = 0; 

