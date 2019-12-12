KernsP_TSTHLR=getKernsFromTrajM(TrajM(RepsForKerns,TrajPartMed),Sz,TSB_THLR);
%%
tic
tmp=getKernsFromTrajM(TrajM(1:2,:),Sz,TSB_THLR(:,1:2));
t=toc
tic
tmp=getKernsFromTrajM(repmat(TrajM(1:1,:),[1 2]),Sz,repmat(TSB_THLR(:,1:1),[2 1]));
t=toc
%%
QQTraj=20+35i;
% QQ=getKernsFromTrajM(QQTraj,Sz,1);
% fgmontagex(log(abs(QQ(round(real(QQTraj)*2+1)+(-11:11),round(imag(QQTraj)*2+1)+(-11:11)))))

[X Y]=ndgrid(0:0.1:(1-0.0001),0:0.1:(1-0.0001));
QQTrajX=[QQTraj QQTraj+0.01 QQTraj+0.005 QQTraj+0.9 QQTraj+0.1i QQTraj+0.9i];

% QQTrajX=Row((QQTraj+X+1i*Y));
CTo2Rows=@(X) [real(X);imag(X)];

% SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(QQTrajX),Sz), Sz, [6 6], Sz*2); % st.om
SnufftStruct_CurRep = nufft_init(BART2Fes_NUFT_Idxs(CTo2Rows(QQTrajX),Sz), Sz, [7 7], Sz*2); % st.om
P=perm31(reshape(full(SnufftStruct_CurRep.p),[numel(QQTrajX) Sz*2]));
grmss(mean(P(:,:,1:2),3)-P(:,:,3))./grmss(P(:,:,3))
Range=4;
tmp=(P(round(imag(QQTraj)*2+1)+(-Range:Range),round(real(QQTraj)*2+1)+(-Range:Range),:));
fgmontage(log(abs(tmp)))