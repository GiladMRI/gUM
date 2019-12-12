load('splitProx-4Echo_demo.mat');
DX=RecX-GTX;

StepBRes=cat(4,GTX,RecA,RecX,DA*ErrShowFac,DAF*ErrShowFac,DX*ErrShowFac);
fgmontagex(StepBRes,[0 1e3]);title(num2str([grmss(DA) grmss(DAF) grmss(DX)],' %.2f '));

% ErrShowFac=5;
StepBRes=cat(4,RecA,RecX,DA*ErrShowFac,DX*ErrShowFac);
fgmontagex(mean(abs(StepBRes),3),[0 1e3],'Size',[2 2]);

DBase=GTX-gmean(GTX);
Denom=norm(DBase(:));
NRMSEA=norm(DA(:))./Denom;
NRMSEX=norm(DX(:))./Denom;
%%
load('splitProx-4Echo_demo_Pois.mat')

DA=RecA-GTX;
DX=RecX-GTX;

% StepBRes=cat(4,GTX,RecA,RecX,DA*ErrShowFac,DAF*ErrShowFac,DX*ErrShowFac);
% fgmontagex(StepBRes,[0 1e3]);title(num2str([grmss(DA) grmss(DAF) grmss(DX)],' %.2f '));

% ErrShowFac=5;
StepBRes=cat(4,RecA,RecX,DA*ErrShowFac,DX*ErrShowFac);
fgmontagex(mean(abs(StepBRes),3),[0 1e3],'Size',[2 2]);

DBase=GTX-gmean(GTX);
Denom=norm(DBase(:));
NRMSEA=norm(DA(:))./Denom;
NRMSEX=norm(DX(:))./Denom;
%%
for i=1:4
    MM(:,:,i)=bart(['poisson -Y 40 -Z 40 -C 4 -y 2 -z 2 -s ' num2str(rand*100000)]);
end
clear MMX
MMX(:,:,:,1)=MM(:,:,1).*permute32([1 0 0]);
MMX(:,:,:,2)=MM(:,:,2).*permute32([0 1 0]);
MMX(:,:,:,3)=MM(:,:,3).*permute32([0 0 1]);
MMX(:,:,:,4)=MM(:,:,4).*permute32([1 0.1 0.6]);
MMX(50,50,1,1)=0;
MMX=circshift(circshift(MMX,5,1),5,2);
figure;imagesc(CombineDims(MMX,[4 2]));daspect([1 1 1]);