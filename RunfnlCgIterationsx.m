param=ExtendStruct(struct('pNorm',1,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',XFM,'TV',TVOP,'xfmWeight',xfmWeight,'ShowFig',ShowFig),init);
param.data =     DataP;
nfnlCgIters=40;
RunFnlViewAmp=1;
res=zeros(Sz2);
FigH=4000;
if(param.ShowFig)
    figure(FigH);close(FigH);
end
StartTime_fnl=now;
param.Verbose=false;
clear ObjConv Score
for n=1:nfnlCgIters
    [res, CurObj] = fnlCg(res,param);
    ObjConv(n)=CurObj;
    im_res = param.XFM'*res;
    if(param.ShowFig)
        figure(FigH); subplot(1,3,1);
        gmontage(abs(gflip(im_res,[]))); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,3,2);
        gmontage(angle(gflip(im_res,[]))); drawnow;% title(qq)
        subplot(1,3,3);
        plot(ObjConv);setYaxis([0 CurObj*3]);if(n>1), setXaxis([1 n]);end
    end
%     t=toc;
    if(n>1)
        dObjNom=ObjConv(n-1)-CurObj;
        dObjP=dObjNom/ObjConv(n-1);
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g') ' dObjP ' num2str(dObjP,'%g')]);
        if(dObjP<2e-3 || dObjNom<1e-16)
            disp('Not advancing. Stopping.');
            break;
        end
    else
        disp(['Iter #' num2str(n,'%02d') ' ' datestr(now) ' ' num2str(CurObj,'%5.3g')]);
    end
end
if(param.ShowFig)
    close(FigH)
end
disp('ok im_res');