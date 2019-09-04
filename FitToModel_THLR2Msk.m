X0NoBN_TH=H_AllTS2*squeeze(In);

[ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);

U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
VH1=permute43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));

U2=U_tmp(:,:,:,2).*(s_tmpP(:,:,1,2));
VH2=permute43(conj(V_tmp(:,:,:,2))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));

LR1Est=permute(sum(U1.*permute(VH1,[1 2 5 3 4]),4),[1 2 3 5 4]);
LR2Est=permute(sum(U2.*permute(VH2,[1 2 5 3 4]),4),[1 2 3 5 4]);
LREst=LR1Est+LR2Est.*Use2ndSVMap;

% Db=X0NoBN_TH-LREst;
% Db_rms(iter)=grmss(Db(~isnan(Db)));
    
Out=(H_AllTS2'*LREst);