% ScriptFN=[BaseSP 'MPBDScriptP.txt'];
ScriptFN=[BaseSP 'MPBDScriptCartMP.txt'];
% Open over PD and then like ITS
MP_Cmnds=[{'fmac 2 0'}, ITS_Cmnds];
WriteLinopToFile(ScriptFN,MP_Cmnds);
% SensCSMap is 0
% mask_sampleP is 1
% PD Decay is 2

MPImSize16=ImSz16;
MPImSize16(TS_Dim)=1;
BARTS=struct();
BARTS.cmd=['picsS -w 1 -R W:3:0:0.001 -m ' ScriptFN];
BARTS.ImSz16=MPImSize16;
BARTS.sig=sig;
BARTS.Others={SensCSMap mask_sample};

BaseFP='/tmp/';

Flds=setdiff(fieldnames(BARTS),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTS.(Flds{i})))
        for j=1:numel(BARTS.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTS.(Flds{i}){j});
            BARTS.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTS.(Flds{i}));
        BARTS.(Flds{i})=tmpFN;
    end
end

BARTS_Aop=struct();
BARTS_Aop.cmd=['linopScript -N ' LS_ScriptFN];
BARTS_Aop.ImSz16=ImSz16;
BARTS_Aop.Others={SensCSMap mask_sample};

BARTStruct=BARTS_Aop;
Flds=setdiff(fieldnames(BARTStruct),'cmd');
for i=1:numel(Flds)
    if(iscell(BARTStruct.(Flds{i})))
        for j=1:numel(BARTStruct.(Flds{i}))
            tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
            writecfl(tmpFN,BARTStruct.(Flds{i}){j});
            BARTStruct.(Flds{i}){j}=tmpFN;
        end
    else
        tmpFN=[BaseFP 'tmp' num2str(randi(1000000000))];
        writecfl(tmpFN,BARTStruct.(Flds{i}));
        BARTStruct.(Flds{i})=tmpFN;
    end
end
BARTS_Aop=BARTStruct;

ksp_adj=bart(['linopScript -A ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

disp('Prepared for BART_Aop');
%%
CurTestElemIdx=[1:4];
% CurTestElemIdx=[1:2];
% CurTestElemIdx=[1 3 4];
% CurTestElemIdx=[3];
% CurTestElemIdx=2;

ElemNames={'m' 'p' 'b' 't'};
nElements=numel(ElemNames);
ElemGT={mGT,pGT,bGT,tGT};
ElemL={M,P,B,TT};

fIdent=@(x) x;
fOne=@(x) 1;
fExp=@(x) exp(x);
fInv=@(x) 1./x;
fMinusInvSqr=@(x) -1./(x.^2);

Elem_u={fIdent,fIdent,fIdent,fInv};
% Elem_du={fOne,fOne,fOne,fMinusInvSqr};
Elem_s={fIdent,fExp,fExp,fExp};
% Elem_ds={fOne,fExp,fExp,fExp};

Elem0=ElemGT;

Distort=true;
if(Distort)
    [X, Y]=ndgrid(linspace(-1,1,N),linspace(-1,1,N));
    cGT=mGT.*exp(1i*pGT);
    SmoothR=21;
    SmoothS=3;
    SmoothR=11;
    SmoothS=1;
    scGT=SmoothBySlices(cGT,SmoothR,SmoothS);
    sbGT=SmoothBySlices(bGT,SmoothR,SmoothS);
    stGT=SmoothBySlices(tGT,SmoothR,SmoothS);

    Elem0{1}=abs(scGT);
    Elem0{2}=angle(scGT);
    Elem0{3}=sbGT;
    Elem0{4}=stGT;

    % for i=1:numel(CurTestElemIdx)
    % %     Elem0{CurTestElemIdx(i)}=SmoothBySlices(ElemGT{CurTestElemIdx(i)},[20 20],3);
    %     Elem0{CurTestElemIdx(i)}=SmoothBySlices(ElemGT{CurTestElemIdx(i)},[11 11],1);
    % end

    % Elem0{1}=SmoothBySlices(ElemGT{1},[20 20],3);
    Elem0{2}=Elem0{2}+X.*1+Y.*0.6-0.2;
    % Elem0{2}=Elem0{2}+X.*3+Y.*2.6-1.5;
    % Elem0{CurTestElemIdx(i)}=SmoothBySlices(ElemGT{CurTestElemIdx(i)},[20 20],1);

    Elem0{2}=angle(exp(1i*Elem0{2}));

    Elem0{3}=Elem0{3}+X.*4+Y.*5.5-1;
end

% if(nCS>1)
%     m0=cat(CS_Dim,m0,m0*0);
%     p0=cat(CS_Dim,p0,p0*0);
%     b0=cat(CS_Dim,b0,b0);
%     t0=cat(CS_Dim,t0,t0*0);
% end

for i=1:nElements
    Elem_sLu{i}=Elem_s{i}(ElemL{i}*(Elem_u{i}(Elem0{i})));
end

Est0=1;
for i=1:nElements
    Est0=Est0.*Elem_sLu{i};
end
Est0X=squeeze(sum(Est0,CS_Dim));

MPBT=RefValueMPBT;
for i=1:nElements
    MPBT.(ElemNames{i})=Elem0{i};
end

disp('ok initials');
%%
% ElemsAlphas=[0.0323, 1.5366e-07, 1.3036e-05, 0.0015];
% AlphaFacMPBT=struct('m',1,'p',100,'b',1e-2,'t',1e7);
AFacMPBT=[1e-0,1000,10,1e7];
% AFacMPBT=[3,1000,10,1e7];
% AFacMPBT=AFacMPBT/2;
% AlphaFacMPBT=struct('m',1e-0,'p',1000,'b',10,'t',1e7); % Very good for each one separately
AlphaFacMPBT=struct('m',AFacMPBT(1),'p',AFacMPBT(2),'b',AFacMPBT(3),'t',AFacMPBT(4));
% AlphaFacMPBT=struct('m',1e-1,'p',100,'b',1,'t',1e6);
% AlphaFacMPBT=struct('m',1e-0,'p',1e+2,'b',10,'t',1e7);
% AlphaFacMPBT=struct('m',1e-1,'p',1e+2,'b',1,'t',1e6);
% AlphaFacMPBT=struct('m',1,'p',-1e6,'b',-1e+0,'t',1e7);
ElemsLambda=[0.015,0.05,0.05,0.0003]; % Very good for each one separately
ElemsLambda=[0.02,0.05,0.05,0.0003];
ElemsLambda=ElemsLambda/10;
% ElemsLambda=-[25e-3,5e-2,0.05,0.0003];
% ElemsLambda=[0.015,0,0.05,0.003];
% ElemsLambda=[0.015,0.05,0.05,300];

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
end
ElemsL={T,1i*T,-1i*2*pi*TDx/1000,-TDx};
for i=1:numel(ElemsL)
%     writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],repmat(ElemsL{i},[N N 1 1 1 1 1]));
end
Mm=M*m0;
expPp0 = exp(P * p0);
expBb0 = exp(B * b0);
expTt0 = exp(TT * (1./t0));

Othersp=Mm.*expBb0.*expTt0;
Othersb=Mm.*expPp0.*expTt0;
Otherst=Mm.*expPp0.*expBb0;
ElemsAlphas=[AlphaFacMPBT.m/ lipschitz(M),AlphaFacMPBT.p/ lipschitz(P) / (max(abs(Othersp(:)))^2 + eps),...
    AlphaFacMPBT.b / lipschitz(B) / (max(abs(Othersb(:)))^2 + eps),...
    AlphaFacMPBT.t / lipschitz(TT) / (max(abs(Otherst(:)))^2 + eps)];

% alphatbase=1e11;
% alphatbase=1e7;
% ElemsAlphas(4)=alphatbase / lipschitz(TT) / (max(abs(Otherst(:)))^2 + eps);
writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
% writecfl([ToBARTP 'ninneriter'],ninneriter.');
% Types: M=1 P=2 D/B=3 T=4 C=5
ElementTypes=[1 2 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');


writecfl([ToBARTP 'sig_adj'],ksp_adj);

ninneriterBART=[0 0 0 0];
for i=1:numel(CurTestElemIdx)
    ninneriterBART(CurTestElemIdx(i))=2;
end
% ninneriterBART(1)=7;
% ninneriterBART(4)=0;
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:4
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%
writecfl([ToBARTP 'GT'],GT);
%%

% Continue?
% for i=1:numel(Elem0)
%     writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
% end
QQ=bart(['splitProx -i 5 -s 60 -t -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
% bart(['splitProx -i 5 -s 60 -t -d 2 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
%%
% save('ForCallx.mat','ToBARTP','BARTS_Aop');
% load('ForCall.mat','ToBARTP','BARTS_Aop');

ContinueRun=true;
if(ContinueRun)
    % Continue?
    for i=1:numel(Elem0)
        writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
    end
end

QQ=bart(['splitProx -i 150 -s 60 -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});
%%
for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
%%
ElemRanges={[0 10],[-pi pi],[-400 400],[0 100]};

ElemRanges={[0 1000],[-pi pi],[-400 400],[0 100]};
mThresh=3e-8;
figure;
for i=1:4
    subplot(2,2,i);
    gmontage(Maps{i},ElemRanges{i});removeTicks;
end
%%
%     alphab = 1.0*AlphaFacMPBT.b / lipschitz(B) / (max(abs(Mm0(:)))^2 + eps) ;
%         expBb0 = exp(B * b0);
%         CurEst=Mm0.*expPp0.*expTt0 .* expBb0;
%         AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
% r=ksp_adj- AHACurEst;
% QQQ=conj(CurEst) .* r;
% AfterAdjL=B' * QQQ;
% Newb=b0+alphab*real(AfterAdjL);
% % b = Pb(b + alphab * real(B' * (conj(Mm) .* conj(expBb) .* r)), alphab);
%         
% XDimsOut=readcfl([ToBARTP 'XDimsOut']);
% mDimsOut=readcfl([ToBARTP 'mDimsOut']);
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<=0,1))-1);
% figure;plot(ErrVec)
%
Maps=cell(1,4);
for i=1:4
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
ElemRanges={[0 10],[-pi pi],[-400 400],[0 100]};
mThresh=50;
mThresh=3e-8;

CurElem=CurTestElemIdx;
ElemRanges={[0 700],[-pi pi],[-400 400],[0 100]};
%
figure;
for i=1:4
    subplot(2,2,i);
%     gmontage(Maps{i},ElemRanges{i});removeTicks;
%     gmontage(Maps{i}.*(m0>mThresh),ElemRanges{i});removeTicks;
%     gmontage(Elem0{i}.*(m0>mThresh),ElemRanges{i});removeTicks;
    gmontage(ElemGT{i}.*(m0>mThresh),ElemRanges{i});removeTicks;
end
%%
CurElem=3;
ElemRanges={[0 700],[-pi pi],[-10 10],[0 100]};
ElemRecs=cat(3,ElemGT{CurElem},Maps{CurElem},Elem0{CurElem});
% ElemRecs=cat(3,ElemGT{CurElem},Maps{CurElem},Elem0{2});
% ElemRecs=angle(exp(1i*ElemRecs));
ElemRecs=ElemRecs-ElemGT{CurElem};
% ElemRecs=ElemRecs.*(m0>50);
fgmontagex(ElemRecs,ElemRanges{CurElem},'Size',[1 3]);title([PadStringWithBlanks('GT',70) PadStringWithBlanks('splitProx',70) 'Start point'],'FontSize',16);


%% Try CBT
% ElemsAlphas=[0.0323 1.3036e-5 0.0015];
ElemsAlphas=[0.000323 0.3036e-6 0.00003];
ElemsAlphas=ElemsAlphas*50;
ElemsLambda=[25e-2,0.05,0.0003];

ToBARTP='/autofs/space/daisy_002/users/Gilad/gUM/';

ElemsL={T,-1i*2*pi*TDx/1000,-TDx};

Elem0{1}=Elem0{1}.*exp(1i*Elem0{2});
Elem0{2}=Elem0{3};
Elem0{3}=Elem0{4};
Elem0=Elem0(1:3);
for i=1:numel(Elem0)
    writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Elem0{i});
    writecfl([ToBARTP 'ElemsL_' num2str(i-1)],ElemsL{i});
end

writecfl([ToBARTP 'ElemsAlpha'],ElemsAlphas.');
writecfl([ToBARTP 'ElemsLambda'],ElemsLambda.');
% Types: M=1 P=2 D/B=3 T=4 C=5
ElementTypes=[5 3 4];
writecfl([ToBARTP 'ElementTypes'],ElementTypes.');

writecfl([ToBARTP 'sig_adj'],ksp_adj);

CurTestElemIdx=[1:3];
% CurTestElemIdx=2;
ninneriterBART=[0 0 0 0];
for i=1:numel(CurTestElemIdx)
    ninneriterBART(CurTestElemIdx(i))=2;
end
writecfl([ToBARTP 'ninneriter'],ninneriterBART);
disp('saved all');

for i=1:3
    delete([ToBARTP 'Elem' num2str(i-1) '.hdr']);
    delete([ToBARTP 'Elem' num2str(i-1) '.cfl']);
end
%
writecfl([ToBARTP 'GT'],GT);
%

% Continue?
% for i=1:numel(Maps)
%     writecfl([ToBARTP 'ElemsWS_' num2str(i-1)],Maps{i});
% end
QQ=bart(['splitProx -i 2000 -s 60 -g -F ' ToBARTP ' ' LS_ScriptFN],BARTS_Aop.ImSz16,sig,BARTS_Aop.Others{:});

% XDimsOut=readcfl([ToBARTP 'XDimsOut']);
% mDimsOut=readcfl([ToBARTP 'mDimsOut']);
ErrVec=readcfl([ToBARTP 'ErrVec']);
ErrVec=ErrVec(1:(find(ErrVec<1,1))-1);
% figure;plot(ErrVec)
%
Maps=cell(1,3);
for i=1:3
    Maps{i}=readcfl([ToBARTP 'Elem' num2str(i-1)]);
end
Maps{4}=angle(Maps{1});
Maps{1}=abs(Maps{1});
Maps=Maps([1 4 2 3]);

CurElem=CurTestElemIdx;
ElemRanges={[0 700],[-pi pi],[-400 400],[0 100]};
%
figure;
for i=1:4
    subplot(2,2,i);
%     gmontage(Maps{i}.*(m0>50),ElemRanges{i});removeTicks;
%     gmontage(Elem0{i}.*(m0>50),ElemRanges{i});title('WS');removeTicks;
    gmontage(ElemGT{i},ElemRanges{i});title('GT');removeTicks;
end