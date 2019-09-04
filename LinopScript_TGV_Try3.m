setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/LSBart/bart-toUpload')
%% TGV denoising
% This is an example of using a multi-linop script for image reconstruction
% with tailored image-to-signal operator and tailored image to constrain
% operators.
%
% Following "Second Order Total Generalized Variation (TGV) for MRI", Knoll
% et al, MRM 2011.
%
% 2D TGV2 can be formulated as follows:
% The "image domain" will be extended (along a dimension denoted "TGV") to include 2 more variables, except
% the image estimate I, also u,and v, the variables controlling if to
% locally use first or second derivative. So finally they'll be similar to x,y derivative of I.
%
% TGV aims to optimize more or less the following function: $$|A*I-signal| + | \nabla(I)-\{u,v\}|_1+| \nabla\{u,v\}|_1 $$
%
% For simplicity this demonstrated image denoising.
% So, we'll have 3 linops:
% 
% 1. the "image-to-signal operator is just choosing I from {I,u,v}, we do that by
% multiplying with [1 0 0] and summing (fmac).78
%
% 2.  $$\nabla(I)-\{u,v\} $$
%
% 3.  $$\nabla(\{u,v\}) $$
%
N=128;
HalfN=N/2;
Signal=phantom(N);
Img=Signal;
Img(:,:,1,1,1,3)=0;
% TGV image-to-signal 
ChooseFirst=[1;0;0];
ChooseFirstP=permute(ChooseFirst,[6 5 4 3 2 1]);
TGV_Denoise_MainLinop_FN='TGV_Denoise_Main.txt';
TGV_Main_Ops={'fmac 0 32'};
WriteLinopToFile(TGV_Denoise_MainLinop_FN,TGV_Main_Ops);
ChooseFirstRes=bart(['linopScript ' TGV_Denoise_MainLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP);
%% First regulatizer linop is choosing u,v
%  and applying symmetric differetiation,
% i.e. du/dx, dv/dy, (du/dy+dv/dx)/2
ChooseLastTwo=[0 0; 1 0; 0 1];
ChooseLastTwoP=permute(ChooseLastTwo,[7 6 5 4 3 1 2]);
Img(:,:,:,:,:,2)=checkerboard(N/8);
Img(:,:,:,:,:,3)=checkerboard(N/8)+1;
TGV_Denoise_uvLinop_FN='TGV_Denoise_uv.txt';
% after choosing 2,3 is X   Y     1     1     1     1     2
% After grad we're X Y 1 1 1 1 2 ... 2
% Multiply that by [1 0.5; 0.5 1]
SymDiffWeight=[1 0.5; 0.5 1];
SymDiffWeightP=permute(SymDiffWeight,[8 7 6 5 4 3 1 16 15 14 13 12 11 10 9 2]);
TGV_uv_Ops={'fmac 1 32','grad 3','fmac 2 0'};
WriteLinopToFile(TGV_Denoise_uvLinop_FN,TGV_uv_Ops);
uvLinopRes=bart(['linopScript ' TGV_Denoise_uvLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP);
fgmontage(uvLinopRes);
%% 2nd regularizer is on dI/dx-u, dI/dy-v
% we work as follows:
% 
% Increase [I,u,v] to [I 2I; u 2u; v 2v], Over a new dim AUX (sized 2) 
% 
% diff along X,Y, AUX. So we're left with
% 
% [X Y 1 1 1 TGV AUX 1 1 1 1 1 1 1 1 GRAD], GRAD of size 3
% 
% looking only at TGV,GRAD dims, we have: (? = not important)
% 
% M= [dI/dx dI/dy    ?                    -> GRAD dim
% 
%       ?     ?      du/dAUX=-u            rows:TGV dim
% 
%       ?     ?      dv/dAUX=-v    ]
% 
% On the AUX dim we'll have [-1, 1] multiplications of this matrix. 
% 
% From this matrix we'll choose dI/dx-u, dI/dx-v, into AUX2 (dim 8) to
% 
% [X Y 1 1 1 1 1 AUX2], where AUX2 contains dI/dx-u, dI/dy-v as needed
AUX=[1;2];
AUXP=permute(AUX,[7 6 5 4 3 2 1]);
clear ChooseFromM
ChooseFromM(:,:,1,2)=[1 0 0; 0 0 1; 0 0 0];
ChooseFromM(:,:,2,2)=[0 1 0; 0 0 0; 0 0 1];
ChooseFromMP=permute(ChooseFromM,[9 8 7 6 5 1 4 3 10 11 12 13 14 15 16 2]); % 2^15+2^5+2^6 =32864
TGV_Denoise_dIuvLinop_FN='TGV_Denoise_dIuv.txt';
TGV_dIuv_Ops={'fmac 3 0','grad 67','fmac 4 32864'};
WriteLinopToFile(TGV_Denoise_dIuvLinop_FN,TGV_dIuv_Ops);
dIuvLinopRes=bart(['linopScript ' TGV_Denoise_dIuvLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
fgmontage(dIuvLinopRes);
%% In one multi-linop file
% Here we create a multi-linop file, and apply each of the linops on the
% input Img, by calling linopScript -L
TGV_Denoise_Linops_FN='TGV_Denoise_linops.txt';
WriteLinopToFile(TGV_Denoise_Linops_FN,{TGV_Main_Ops,TGV_uv_Ops,TGV_dIuv_Ops });

V1=bart(['linopScript -L 0 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
V2=bart(['linopScript -L 1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(V1)),V1,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
V3=bart(['linopScript -L 2 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(V2)),V2,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);

size(V1)
size(V2)
size(V3)
%% TGV denoising example
% Now we call the parallel-imaging compressed-sensing recon with
% linopScript by calling picsS, with the multi-linop Script file.
% The first linop will be used as image-to-signal operator,
% and the rest will used for the regularizers, by their order. So, the 1st
% additional linop will be used for the first regularized defined in the
% picsS -R part.
%
% We follow the authors' suggestion and set the regularization for the
% first term to be double that of the second.
[X1 Y1]=meshgrid(1:N,1:N);
[X2 Y2]=meshgrid(1:2:N,1:2:N);
TGVexample=X1;
TGVexample((N/4)+(1:HalfN),(N/4)+(1:HalfN))=fliplr(X2);

Noised=TGVexample+randn(N)*10;

IdentitiyFN='Identity.txt';
WriteLinopToFile(IdentitiyFN,'ident');

DenoisedTV=bart(['picsS -R T:3:3:4 ' IdentitiyFN],FillOnesTo16(size(Noised)),Noised);
DenoisedTVx=bart(['picsS -R T:3:3:8 ' IdentitiyFN],FillOnesTo16(size(Noised)),Noised);

DenoisedTGV=bart(['picsS -m -R 1:0:2 -R 1:131:1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Noised,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
DenoisedTGVx=bart(['picsS -m -R 1:0:4 -R 1:131:2 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Noised,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);

%% Resulting denoised images
fgmontage(cat(3,DenoisedTV,DenoisedTGV(:,:,1)));title('TV                                                   TGV','FontSize',16);

%% Comparison
% while both TV and TGV show similar smoothness, the TGV preserves the
% edges better. So maybe some of the math used here was right.
figure;
plot(squeeze(DenoisedTGV(HalfN,:,1)),'k');hold on;
plot(squeeze(DenoisedTGVx(HalfN,:,1)),'r');hold on;
plot(squeeze(DenoisedTV(HalfN,:,1)),'k:','LineWidth',2);hold on;
plot(squeeze(DenoisedTVx(HalfN,:,1)),'r:','LineWidth',2);hold on;
legend({'TGV weak','TGV strong','TV weak','TV strong'},'Location','Best')