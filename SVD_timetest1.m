Sz=[128 128];
N=10000;
NN=prod(Sz);
A=randn([N,NN])+1i*randn([N,NN]);
NSV=200;
%%
disp(datestr(now));
tic
[U,S,V]=svds(A,NSV);
t=toc;
disp([num2str(Sz) ' ' num2str(N) ' ' num2str(NSV) ' ' num2str(t)]);
% tiger:
% 128  128 4000 50 55.6665
% 128  128 4000 200 94.4159
% 128  128 10000 200 255.1676
%% Prepare database
% HCPData_256x256_int16
% random flip rot crop 
RandomPhaseLinearFac=2;
RandomPhaseQuadraticFac=0.05;
RandomPhaseScaleFac=1.4;
P=zeros([N Sz]);
for i=1:N
    disp(i);
    P(i,:,:)=GenerateRandomSinPhase(Sz,RandomPhaseLinearFac,RandomPhaseQuadraticFac,RandomPhaseScaleFac);
end
%%
HCP=load('/autofs/cluster/kawin/Gilad/HCPData_256x256_int16');
HCP=HCP.HCPData(:,:,1:N);
%% magnitude
SzHCP=gsize(HCP,1:2);
CropN=SzHCP(1)-Sz(1);
M=single(zeros([N Sz]));
for i=1:N
    disp(i);
    tmp=squeeze(single(HCP(:,:,i)));
    if(rand<0.5)
        tmp=gflip(tmp,1);
    end
    if(rand<0.5)
        tmp=gflip(tmp,2);
    end
    if(rand<0.5)
        tmp=rot90(tmp);
    end
    SIdx=randi(CropN,[1 2]);
    M(i,:,:)=tmp(SIdx(1)+(1:Sz(1)),SIdx(2)+(1:Sz(2)));
end
%%
C=double(M.*P);
C=reshape(C,[N prod(Sz)]);
NSVc=10000;
disp(datestr(now));
tic
[~,S,V]=svds(C,NSVc);
t=toc;
disp([num2str(Sz) ' ' num2str(N) ' ' num2str(NSVc) ' ' num2str(t)]);
save('SV128x128x10k_10k','S','V');
%%
VM=reshape(conj(V),[Sz NSVc]);
S=diag(S);
save('SVx128x128x10k','S','VM');
%%
VMS=VM.*perm31(S);
%% Recon with those
NSVsToUse=500;
SVScriptFN=[pwd filesep 'SCScript.txt'];
ResSz16=FillOnesTo16([ones(1,6) NSVsToUse]);
Acc=2;
CenterSize=4;
Msk=squeeze(bart(['poisson -Y ' num2str(Sz(1)) ' -Z ' num2str(Sz(2)) ' -y ' num2str(Acc) ' -z ' num2str(Acc) ' -C ' num2str(CenterSize)]));
Msk=Msk*0+1;
Sens=1;
% sensfile is 0
% Msk is 1
% SV is 2
% SV dim is 7, .e. 64
WriteLinopToFile(SVScriptFN,{'fmac 2 64','fmac 0 16','fftc 3','fmac 1 0'});
SVP=perm73(VMS(:,:,1:NSVsToUse));

I=phantom(128);
FI=fft2cg(I);
FIM=FI.*Msk;

Rec=bart(['picsS -w 1 -d 4  ' SVScriptFN],ResSz16,FIM,Sens,Msk,SVP);
RecX=sum(Rec.*SVP,7);
%%
L2Lambda=0.01;
Rec=bart(['picsS -m -w 1 -d 4 -R 2:' num2str(L2Lambda) ' ' SVScriptFN],ResSz16,FIM,Sens,Msk,SVP);
%%
ICoeffs=gsum(I.*conj(VM),1:2);
IRec=sum(ICoeffs.*VM,3);