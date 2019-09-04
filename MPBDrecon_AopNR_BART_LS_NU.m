function [MPBD, Relr]= MPBDrecon_AopNR_BART_LS_NU(ksp_adj, BARTS_Aop, M, P, B, D, Pm, Pp, Pb, Pd, MPBD, W, niter, ninneriter, dohogwild, doplot,GT, BARTS, BARTk,AlphaFacMPBD)
global Relr

CS_Dim=5;

ksp_adj_rms=grmss(ksp_adj);

if(~exist('GT','var'))
    GT=[];
end

DoW=(numel(W)==1);
if(DoW)
    nW=W;
end

m=MPBD.m;
p=MPBD.p;
b=MPBD.b;
d=MPBD.d;

MmF=M*m;
expPpF = exp(P * p);
expDdF = exp(D * d);
expBbF = exp(B * b);

if(~isempty(GT))
    EstF=MmF.*expPpF.*expDdF.*expBbF;
    EstF2=sum(EstF,CS_Dim);
    disp(['Start ErrGT: ' num2str(grmss(GT-EstF2)/grmss(GT),'%.3g')]);
end

k = 0;
K = 10;
h = 1;

giter=0;
% Start iteration
for it = 1:niter
    disp([num2str(it) ' ' datestr(now)]);
    
    if(DoW)
        W=cell(1,nW);
        for i=1:nW
            tmp=angle(exp(1i*(p + (rand*2-1)*pi)))-p;
            W{i}=tmp;
        end
    end
    
    if(mod(it-1,BARTk)==0)
        expDdF = exp(D * d);
        expBbF = exp(B * b);
        
        BD=expDdF.*expBbF;
%         PDP=permute(BD,[1 2 3 6 7 4 5]);

%         RecM=bart(BARTS.cmd,BARTS.ImSize16,BARTS.yP,BARTS.TP,PDP,BARTS.Others{:});
%         RecM=bart(BARTS.cmd,BARTS.ImSz16,BARTS.sig,BARTS.Others{:},BD);
        RecM=bart(BARTS.cmd,BARTS.ImSz16,BARTS.sig,BARTS.Others{1:3},BD,BARTS.Others{5});
        
        m=abs(RecM);
        p=angle(RecM);
        
        MPBD.m=m;
        MPBD.p=p;
        
        if(~isempty(GT))
            MmF=M*m;
            expPpF = exp(P * p);
            
            EstF=MmF.*expPpF.*expDdF.*expBbF;
            EstF2=sum(EstF,CS_Dim);
            disp(['ErrGT after BART: ' num2str(grmss(GT-EstF2)/grmss(GT),'%.3g')]);
            aa=5;
        end
    else
        expPp = exp(P * p);
        expBb = exp(B * b);
        expDd = exp(D * d);
        expPp = expPp  .* expDd .* expBb;
        alpham = 1.0 *AlphaFacMPBD.m/ lipschitz(M) * h;
        disp('m iters');
        for itinner = 1:ninneriter(1)
            Mm = M * m;
            CurEst=Mm .* expPp;
            AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
            r=ksp_adj- AHACurEst;
%             disp(grmss(r)/ksp_adj_rms);
%             r = A' * (y - A * (Mm .* expPp));
            tmp=real(M' * (conj(expPp) .* r ));
            giter=giter+1;
            Relr(giter)=grmss(r)/ksp_adj_rms;
            disp([grmss(r)/ksp_adj_rms grmss(alpham*tmp)/grmss(m)]);
            m = Pm(m + alpham * tmp, alpham);
            
            % Plot
            if (doplot)
                figure(32);
                subplot(2,2,1),
                imshow3(abs(m))
                titlef(it);
                subplot(2,2,2),
                imshow3(real(p))
                titlef(it);
                subplot(2,2,3),
                imshow3(real(b))
                titlef(it);
                subplot(2,2,4),
                imshow3(real(d))
                titlef(it);
                subplot(2,2,1),
                titlef(itinner);
                drawnow
            end
            
            if(size(m,3)>1)
                BB=abs(m(:,:,1))>abs(m(:,:,2));
                tmp=m;
                m(:,:,1)=tmp(:,:,1).*BB+tmp(:,:,2).*(1-BB);
                m(:,:,2)=tmp(:,:,2).*BB+tmp(:,:,1).*(1-BB);
                tmp=d;
                d(:,:,1)=tmp(:,:,1).*BB+tmp(:,:,2).*(1-BB);
                d(:,:,2)=tmp(:,:,2).*BB+tmp(:,:,1).*(1-BB);
                tmp=p;
                p(:,:,1)=tmp(:,:,1).*BB+tmp(:,:,2).*(1-BB);
                p(:,:,2)=tmp(:,:,2).*BB+tmp(:,:,1).*(1-BB);
                tmp=b;
                b(:,:,1)=tmp(:,:,1).*BB+tmp(:,:,2).*(1-BB);
                b(:,:,2)=tmp(:,:,2).*BB+tmp(:,:,1).*(1-BB);
                expPp = exp(P * p);
                expDd = exp(D * d);
                expBb = exp(B * b);
                expPp = expPp  .* expDd .* expBb;
            end
        end
        MPBD.m=m;
        
        Mm = M * m;
        Mm=Mm.*expDd.*expBb;
        alphap = 1.0 *AlphaFacMPBD.p/ lipschitz(P) / (max(abs(Mm(:)))^2 + eps) * h;
        disp('p iters');
        for itinner = 1:ninneriter(2)
            % Get random phase wraps
            if isempty(W)
                w = 0;
            else
                t = randi(length(W));
                w = W{t};
            end
            
            expPp = exp(P * p);
            CurEst=Mm .* expPp;
            AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
            r=ksp_adj- AHACurEst;
%             disp(grmss(r)/ksp_adj_rms);
%             r = A' * (y - A * (Mm .* expPp));
            tmp=real (P' * (conj(Mm) .* conj(expPp) .* r));
%             disp(grmss(tmp)/grmss(p));
            giter=giter+1;
            Relr(giter)=grmss(r)/ksp_adj_rms;
            disp([grmss(r)/ksp_adj_rms grmss(alphap*tmp)/grmss(p)]);
%             p = Pp(p + alphap * tmp + w, alphap) - w;
            CurW=rand*2*pi-pi;
            ppw=angle(exp(1i*(p+CurW)));
            ppwp = Pp(ppw + alphap * tmp, alphap);
            p=angle(exp(1i*(ppwp-CurW)));
            
            % Plot
            if (doplot)
                figure(32);
                subplot(2,2,1),
                imshow3(abs(m))
                titlef(it);
                subplot(2,2,2),
                imshow3(real(p))
                titlef(it);
                subplot(2,2,3),
                imshow3(real(b))
                titlef(it);
                subplot(2,2,4),
                imshow3(real(d))
                titlef(it);
                subplot(2,2,2),
                titlef(itinner);
                drawnow
            end
        end
        MPBD.p=p;
    end
    
    % Now B0
    Mm = M * m;
    expPp = exp(P * p);
    expBb = exp(B * b);
    expDd = exp(D * d);
    Mm=Mm.*expPp.*expDd;
    alphab = 1.0*AlphaFacMPBD.b / lipschitz(B) / (max(abs(Mm(:)))^2 + eps) * h;
    disp('b iters');
    for itinner = 1:ninneriter(3)
        expBb = exp(B * b);
        CurEst=Mm .* expBb;
        AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
        r=ksp_adj- AHACurEst;
%         disp(grmss(r)/ksp_adj_rms);
%         r = A' * (y - A * (Mm .* expBb));
        tmp=real(B' * (conj(Mm) .* conj(expBb) .* r));
%         disp(grmss(tmp)/grmss(b));
        giter=giter+1;
        Relr(giter)=grmss(r)/ksp_adj_rms;
        disp([grmss(r)/ksp_adj_rms grmss(alphab*tmp)/grmss(b)]);
        b = Pb(b + alphab * tmp, alphab);
        
        % Plot
        if (doplot)
            figure(32);
            subplot(2,2,1),
            imshow3(abs(m))
            titlef(it);
            subplot(2,2,2),
            imshow3(real(p))
            titlef(it);
            subplot(2,2,3),
            imshow3(real(b))
            titlef(it);
            subplot(2,2,4),
            imshow3(real(d))
            titlef(it);
            subplot(2,2,3),
            titlef(itinner);
            drawnow
        end
    end
    MPBD.b=b;
    
    % Now decay
    Mm = M * m;
    Mm=Mm.*expPp.*expBb;
    alphad = 1.0*AlphaFacMPBD.d / lipschitz(D) / (max(abs(Mm(:)))^2 + eps) * h;
    disp('d iters');
    for itinner = 1:ninneriter(4)
        expDd = exp(D * d);
        CurEst=Mm .* expDd;
        AHACurEst= bart(BARTS_Aop.cmd,BARTS_Aop.ImSz16,CurEst,BARTS_Aop.Others{:});
        r=ksp_adj- AHACurEst;
%         disp(grmss(r)/ksp_adj_rms);
%         r = A' * (y - A * (Mm .* expDd));
        tmp=real(D' * (conj(Mm) .* conj(expDd) .* r));
%         disp(grmss(tmp)/grmss(d));
        giter=giter+1;
        Relr(giter)=grmss(r)/ksp_adj_rms;
        disp([grmss(r)/ksp_adj_rms grmss(alphad*tmp)/grmss(d)]);
        d = Pd(d + alphad * tmp, alphad);
        
        % Plot
        if (doplot)
            figure(32);
            subplot(2,2,1),
            imshow3(abs(m))
            titlef(it);
            subplot(2,2,2),
            imshow3(real(p))
            titlef(it);
            subplot(2,2,3),
            imshow3(real(b))
            titlef(it);
            subplot(2,2,4),
            imshow3(real(d))
            titlef(it);
            subplot(2,2,4),
            titlef(itinner);
            drawnow
        end
    end
%     d=min(d,0);
    d=max(d,1e-4);
    
    MPBD.d=d;
    
    if dohogwild
        k = k + 1;
        if k == K
            k = 0;
            K = K * 2;
            h = h / 2;
        end
    end
    
    if(~isempty(GT))
        MmF=M*m;
        expPpF = exp(P * p);
        expDdF = exp(D * d);
        expBbF = exp(B * b);
        
        EstF=MmF.*expPpF.*expDdF.*expBbF;
        EstF2=sum(EstF,CS_Dim);
        disp(['ErrGT: ' num2str(grmss(GT-EstF2)/grmss(GT),'%.3g')]);
        aa=5;
    end
end