A=rgb2gray(imread('/home/g/Downloads/SBbatch017314_out_0.0208.png'));
CurI=7;

H=896;
W=1792;
B=A(H*(CurI-1)+(1:H),:);
C=B(1:128,1:(W/2));
CP=B(1:128,W+(1:(W/2)));
CC=[C;CP];
% fgmontage(CC)
% axis equal
% setYaxis([0 128*2])
%
D=B([1:256 128*5+(1:128)],W/2+(1:(W/2)));
% fgmontage(D)
% axis equal
% setYaxis([0 128*3])
E=[CC;D];
fgmontage(E)
axis equal
setYaxis([0 128*5])
%%
A=rgb2gray(imread('/home/g/Downloads/MBbatch022715_out_0.0383.png'));
%%
CurI=5;

H=480;
W=2688;
HH=96;
B=A(H*(CurI-1)+(1:H),:);
C=B(1:HH,1:(W/2));
CP=B(1:HH,W+(1:(W/2)));
CC=[C;CP];
% fgmontage(CC)
% axis equal
% setYaxis([0 128*2])
%
D=B([1:HH*2 HH*4+(1:HH)],W/2+(1:(W/2)));
% fgmontage(D)
% axis equal
% setYaxis([0 128*3])
E=[CC;D];
fgmontage(E)
axis equal
setYaxis([0 HH*5])