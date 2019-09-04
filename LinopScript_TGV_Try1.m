% Image
% Channels
% 
% ChooseFirst=[1;0;0];
% ChooseTwoThree=[0 0; 1 0; 0 1];
%%
setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/LSBart/bart-toUpload')
%% TGV denoising
% image is [X Y 1 1 1 TGV] where TGV is I,u,v
% "Signal" is simply the image, [X Y]
Signal=phantom(128);
Img=Signal;
Img(:,:,1,1,1,3)=0;
% TGV image-to-signal is just choosing the first of TGV, we do that by
% multiplying with [1 0 0]
ChooseFirst=[1;0;0];
ChooseFirstP=permute(ChooseFirst,[6 5 4 3 2 1]);
TGV_Denoise_MainLinop_FN='TGV_Denoise_Main.txt';
TGV_Main_Ops={'fmac 0 32'};
WriteLinopToFile(TGV_Denoise_MainLinop_FN,TGV_Main_Ops);
ChooseFirstRes=bart(['linopScript ' TGV_Denoise_MainLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP);
%% First regulatizer linop is choosing u,v
ChooseLastTwo=[0 0; 1 0; 0 1];
ChooseLastTwoP=permute(ChooseLastTwo,[7 6 5 4 3 1 2]);
Img(:,:,:,:,:,2)=checkerboard(16);
Img(:,:,:,:,:,3)=checkerboard(16)+1;
TGV_Denoise_uvLinop_FN='TGV_Denoise_uv.txt';
TGV_uv_Ops={'fmac 1 32'};
WriteLinopToFile(TGV_Denoise_uvLinop_FN,TGV_uv_Ops);
uvLinopRes=bart(['linopScript ' TGV_Denoise_uvLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP);
fgmontage(uvLinopRes);
%% Better du/dx, dv/dy, (du/dy+dv/dx)/2
ChooseLastTwo=[0 0; 1 0; 0 1];
ChooseLastTwoP=permute(ChooseLastTwo,[7 6 5 4 3 1 2]);
Img(:,:,:,:,:,2)=checkerboard(16);
Img(:,:,:,:,:,3)=checkerboard(16)+1;
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
% Increase [I,u,v] to [I 2I; u 2u; v 2v], Over a new dim AUX (sized 2) 
% diff along X,Y, AUX. So we're left with
% [X Y 1 1 1 TGV AUX 1 1 1 1 1 1 1 1 GRAD], GRAD of size 3
% looking only at TGV,GRAD dims, we have: (? = not important)
% M= [dI/dx dI/dy    ?                    -> GRAD dim
%       ?     ?      du/dAUX=-u            rows:TGV dim
%       ?     ?      dv/dAUX=-v    ]
% On the AUX dim we'll have [-1, 1] multiplications of this mat. 
% From this mat we'll choose dI/dx-u, dI/dx-v, into AUX2 (dim 8) to
% [X Y 1 1 1 1 1 AUX2], where AUX2 contains dI/dx-u, dI/dy-v as needed
AUX=[1;2];
AUXP=permute(AUX,[7 6 5 4 3 2 1]);
clear ChooseFromM
ChooseFromM(:,:,1,2)=[1 0 0; 0 0 1; 0 0 0];
ChooseFromM(:,:,2,2)=[0 1 0; 0 0 0; 0 0 1];
ChooseFromMP=permute(ChooseFromM,[9 8 7 6 5 1 4 3 10 11 12 13 14 15 16 2]); % 2^15+2^5+2^6 =32864
TGV_Denoise_dIuvLinop_FN='TGV_Denoise_dIuv.txt';
% TGV_dIuv_Ops={'fmac 2 0','grad 67','fmac 3 32864'};
TGV_dIuv_Ops={'fmac 3 0','grad 67','fmac 4 32864'};
WriteLinopToFile(TGV_Denoise_dIuvLinop_FN,TGV_dIuv_Ops);
% fid=fopen(TGV_Denoise_dIuvLinop_FN,'wt');fprintf(fid,"fmac 2 0\ngrad 67\nfmac 3 32864\n");fclose(fid);
% fid=fopen(TGV_Denoise_dIuvLinop_FN,'wt');fprintf(fid,"fmac 2 0\ngrad 67\n");fclose(fid);
% fid=fopen(TGV_Denoise_dIuvLinop_FN,'wt');fprintf(fid,"fmac 2 0\n");fclose(fid);
dIuvLinopRes=bart(['linopScript ' TGV_Denoise_dIuvLinop_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
fgmontage(dIuvLinopRes);
%% In one multi-linop file
TGV_Denoise_Linops_FN='TGV_Denoise_linops.txt';
WriteLinopToFile(TGV_Denoise_Linops_FN,{TGV_Main_Ops,TGV_uv_Ops,TGV_dIuv_Ops });

V1=bart(['linopScript -L 0 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Img,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
V2=bart(['linopScript -L 1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(V1)),V1,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);
V3=bart(['linopScript -L 2 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(V2)),V2,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);

size(V1)
size(V2)
size(V3)
%% Now picsS
Img=Signal;
Img(:,:,1,1,1,3)=0;

Rec=bart(['picsS -m -R 1:3:67:0.2 -R 1:3:131:0.1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Signal,ChooseFirstP,ChooseLastTwoP,AUXP,ChooseFromMP);
%% TGV denoising example
[X1 Y1]=meshgrid(1:128,1:128);
[X2 Y2]=meshgrid(1:2:128,1:2:128);
TGVexample=X1;
TGVexample(32+(1:64),32+(1:64))=fliplr(X2);

Noised=TGVexample+randn(128)*10;

IdentitiyFN='Identity.txt';
WriteLinopToFile(IdentitiyFN,'ident');

DenoisedTV=bart(['picsS -R T:3:3:4 ' IdentitiyFN],FillOnesTo16(size(Noised)),Noised);

DenoisedTGV=bart(['picsS -m -R 1:0:2 -R 1:131:1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Noised,ChooseFirstP,ChooseLastTwoP,SymDiffWeightP,AUXP,ChooseFromMP);

% DenoisedTGV=bart(['picsS -m -R T:3:67:1 -R 1:3:131:1 ' TGV_Denoise_Linops_FN],FillOnesTo16(size(Img)),Noised,ChooseFirstP,ChooseLastTwoP,AUXP,ChooseFromMP);
%% GRAD test
I=phantom(200);
GradSciptFN='grad.txt';
fid=fopen(GradSciptFN,'wt');fprintf(fid,"grad 3\n");fclose(fid);
DelI=bart(['linopScript ' GradSciptFN],FillOnesTo16(size(I)),I);
%%
I=[1;3;7];
GradSciptFN='grad.txt';
fid=fopen(GradSciptFN,'wt');fprintf(fid,"grad 1\n");fclose(fid);
DelI=bart(['linopScript ' GradSciptFN],FillOnesTo16(size(I)),I);
% DelI = [-6;2;4]