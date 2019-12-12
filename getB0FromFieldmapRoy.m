BaseP='/media/a/DATA/FromRoy/S034/niftis/fieldmap/';
M0='RoyHaa_200616_S034_20160620_001_007_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E00_M.nii.gz';
P0='RoyHaa_200616_S034_20160620_001_009_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E00_P.nii.gz';
M1='RoyHaa_200616_S034_20160620_001_011_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E01_M.nii.gz';
P1='RoyHaa_200616_S034_20160620_001_013_BP_fieldmap32Ch_MGE3_G4_BP_fieldmap32Ch_MGE3_G4_E01_P.nii.gz';

M0I=loadniidata([BaseP M0]);
P0I=loadniidata([BaseP P0]);
C0=M0I.*exp(1i*2*pi*P0I/4095);

M1I=loadniidata([BaseP M1]);
P1I=loadniidata([BaseP P1]);
C1=M1I.*exp(1i*2*pi*P1I/4095);

I=cat(4,C0,C1);
%%
nSlices=size(I,3);
%%
TEs_us=[4600 9340];
%%
save([BaseP FN '.mat'],'I','sTwix');

%%
WhichTwo=[1 2];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
Mag=squeeze(sqrt(prod(abs(I(:,:,:,WhichTwo)),4)));

Combined=Mag.*exp(1i*angle(M));
gammaMHz=42.5774806;
gammaHz=gammaMHz*1e6;

deltaTE_us=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle=double(angle(Combined));
B0_Hz=dAngle/(2*pi*deltaTE_us/1e6);

Mg=Mag;
%%
WhichSli=1:nSlices;
WhichSli=1:38;
% WhichSli=15:20;
[unwrapped] = cusackUnwrap(dAngle(:,:,WhichSli), Mg(:,:,WhichSli)/3000);
% fgmontage(unwrapped)

B0_HzU=unwrapped/(2*pi*deltaTE_us/1e6);
fgmontage(B0_HzU)
%%
SliI=16;

ToFit=B0_HzU(:,:,SliI);
W=abs(Mg(:,:,SliI));

ToFit(1,:)=ToFit(2,:);

Dx=diff(ToFit,[],1)>30;
Dx=[Dx(1,:); Dx];

Dy=diff(ToFit,[],2)>30;
Dy=[Dy(:,1) Dy];
D=Dx | Dy;
D=imdilate(D,ones(2));
W(D)=0;
%%
fgmontage(cat(3,ToFit,W));
%%
[x,y]=meshgrid(1:128,1:128);
xF=x(:);
yF=y(:);
wF=W(:);
zF=ToFit(:);

DForFit=3;
x=xF(1:DForFit:end);
y=yF(1:DForFit:end);
z=zF(1:DForFit:end);
w=wF(1:DForFit:end);
fo = fitoptions('Weights',w,'Method','LowessFit','Span',0.01);
tic
sf = fit([x, y],z,'lowess',fo)
% sf = fit([x, y],z,fo)
toc
% figure;plot(sf,[x,y],z)
tic
X=sf([xF(1:1:end),yF(1:1:end)]);
toc
X2=reshape(X,[128 128]);

fgmontage(ToFit,[-800 400])
fgmontage(X2,[-800 400])




%%
Mg(Mg==0)=eps;
WhichSli=15;
for WhichSli=1:nSlices
    [PhiCostantini] = cunwrap(dAngle(:,:,WhichSli), struct('weight',Mg(:,:,WhichSli)/3000,'RoundK',false,'maxblocksize',10));
    DA=angle(gsum(Mg(:,:,WhichSli).*exp(1i*(PhiCostantini-dAngle(:,:,WhichSli)))));
    
    PhiCostantiniB=PhiCostantini-DA;
    L=round((PhiCostantiniB-dAngle(:,:,WhichSli))/(2*pi));
    UL=unique(L);
    clear Scr
    for i=1:numel(UL)
        B=L==UL(i);
        Curmg=Mg(:,:,WhichSli);
        Scr(i)=sum(Curmg(B));
    end
    [~,MI]=max(Scr);
    PhiCostantiniB=PhiCostantiniB-UL(MI)*2*pi;
    
    
%     PhiCostantiniM(:,:,WhichSli)=PhiCostantini;
    PhiCostantiniBM(:,:,WhichSli)=PhiCostantiniB;
    % fgmontage(angle(exp(1i*dAngle(:,:,WhichSli))),[-pi pi])
    % fgmontage(angle(exp(1i*PhiCostantiniB)),[-pi pi])
    %
    % fgmontage(PhiCostantiniB,[-2*pi 2*pi])
    % fgmontage(dAngle(:,:,WhichSli),[-2*pi 2*pi])
    B0_HzU(:,:,WhichSli)=PhiCostantiniB/(2*pi*deltaTE_us/1e6);
end
%%
fgmontage(PhiCostantini)
%%
Cx=padRight(padLeft(Combined(:,:,WhichSli),5,2),6,2);
B0_HzUx=padRight(padLeft(B0_HzU,5,2),6,2);
B0Q=imresizeBySlices(B0_HzUx,Sz2);
B0Q=imresizeBySlices(B0_HzU,SnsSzB);
%%
for i=1:numel(sTwix.hdr.Phoenix.sSliceArray.asSlice)
    LocFieldMap(i)=sTwix.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition.dTra;
end
%%
%%
fgmontage(permute(Mg,[3 1 2]))
%%
SnsSz=[96 96];
for SliI=1:nSlices
    disp(SliI);
    Sens(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSz);
end
save([BaseP FN '.mat'],'I','sTwix','B0_Hz','Sens');
%%
SnsSzB=[128 128];
for SliI=1:nSlices
    disp(SliI);
    SensB(:,:,:,SliI)=RunESPIRiTForSensMaps(FirstEcho(:,:,:,SliI),0,SnsSzB);
end
%%
save([BaseP FN '.mat']);
%%
WhichTwo=[1 2];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
Mag=squeeze(mean(abs(I(:,:,:,WhichTwo,:)),4));

Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));

deltaTE_us=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us/1e6)));
B0_Hz=dAngle/(2*pi*deltaTE_us/1e6);
TwoPi_Hz=1e6/(deltaTE_us);
%%
WhichTwo=[2 3];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
Mag=squeeze(mean(abs(I(:,:,:,WhichTwo,:)),4));

Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));

deltaTE_us2=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle2=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us2/1e6)));
B0_Hz2=dAngle2/(2*pi*deltaTE_us2/1e6);
TwoPi_Hz2=1e6/(deltaTE_us2);
%%
WhichTwo=[3 4];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
Mag=squeeze(mean(abs(I(:,:,:,WhichTwo,:)),4));

Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));

deltaTE_us3=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle3=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us3/1e6)));
B0_Hz3=dAngle3/(2*pi*deltaTE_us3/1e6);
TwoPi_Hz3=1e6/(deltaTE_us3);
%%
WhichTwo=[1 3];

M=squeeze(I(:,:,:,WhichTwo(1),:)./I(:,:,:,WhichTwo(2),:));
Mag=squeeze(mean(abs(I(:,:,:,WhichTwo,:)),4));

Combined=squeeze(sum(Mag.*exp(1i*angle(M)),3));

deltaTE_us4=TEs_us(WhichTwo(2))-TEs_us(WhichTwo(1));

dAngle4=double(angle(Combined.*exp(1i*0*2*pi*scanFreq_Hz*deltaTE_us4/1e6)));
B0_Hz4=dAngle4/(2*pi*deltaTE_us4/1e6);
TwoPi_Hz4=1e6/(deltaTE_us4);
%%
ddeltaTE_us=deltaTE_us2-deltaTE_us;
ddAngle=angle(exp(1i*(dAngle2-dAngle+0*2*pi*scanFreq_Hz*ddeltaTE_us/1e6)));
B0_Hzx=ddAngle/(2*pi*ddeltaTE_us/1e6);