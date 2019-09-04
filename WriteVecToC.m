FNBase='SpiGradVec';

% M=rand(3,12);
% M=[GAll(:,:,1) GAll(:,:,2) GAll(:,:,3)].';
M=CombineDims(GAll,[3 2]).';
MaxVals=max(abs(M),[],2);
M=M./MaxVals;
M=[real(M) imag(M)];
M(:,end)=MaxVals;

[nRows,nCols]=size(M);
fid=fopen([FNBase '.cpp'],'wt');
fprintf(fid,'#include "%s.h"\r\nstatic float %s[%d][%d] = {\n\n{',FNBase,FNBase,nRows,nCols);
for i=1:nRows
    for j=1:nCols-1
        fprintf(fid,'%.5ff,',M(i,j));
    end
    if(i==nRows)
        fprintf(fid,'%.5ff}',M(i,end));
    else
        fprintf(fid,'%.5ff},\n{',M(i,end));
    end
end
fprintf(fid,'};\nvoid get%s(float *&Out, long l) {\nOut=%s[l]; };\n',FNBase,FNBase);
fclose(fid);

fid=fopen([FNBase '.h'],'wt');
fprintf(fid,'#ifndef ___%s__\r\n#define ___%s__ 1\nvoid get%s(float *&Out, long l);\n#endif',FNBase,FNBase,FNBase);
fclose(fid);