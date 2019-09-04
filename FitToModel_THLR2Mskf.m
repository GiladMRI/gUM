function [Out, Db_rms]=FitToModel_THLR2Mskf(In,H,Use2ndSVMap)

% syms a b c d w
% eqn=(a.*w+c.*(1-w))^2==b.*w+d.*(1-w);
% sol=solve(eqn,w);
% qs=a*sol+c*(1-sol);
% sqs=simplify(qs)
% %   (b - d + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
% %  -(d - b + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
% (b - d + (4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2)^(1/2))/(2*(a - c))
% Denom=(2*(a - c))
% BasePart=b-d
% Discriminant=4*a^2*d - 4*a*b*c - 4*a*c*d + b^2 + 4*b*c^2 - 2*b*d + d^2
% 
% Sol1=(BasePart+sqrt(Discriminant))/Denom;
% Sol2=(BasePart-sqrt(Discriminant))/Denom;

X0NoBN_TH=H*squeeze(In);

[ U_tmp, s_tmp, V_tmp] = batch_svd(X0NoBN_TH);
s_tmpP=permute43(s_tmp);%     XX=permute(X0N,[1 2 7 6 5 4 3]);

U1=U_tmp(:,:,:,1).*(s_tmpP(:,:,1,1));
VH1=permute43(conj(V_tmp(:,:,:,1))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));

U2=U_tmp(:,:,:,2).*(s_tmpP(:,:,1,2));
VH2=permute43(conj(V_tmp(:,:,:,2))); % .*sqrt(s_tmpP(:,:,1,1:nHComps))));

LR1Est=permute(sum(U1.*permute(VH1,[1 2 5 3 4]),4),[1 2 3 5 4]);
LR2Est=permute(sum(U2.*permute(VH2,[1 2 5 3 4]),4),[1 2 3 5 4]);
LREst=LR1Est+LR2Est.*Use2ndSVMap;

Db=X0NoBN_TH-LREst;
Db_rms=grmss(Db(~isnan(Db)));
    
Out=(H'*LREst);