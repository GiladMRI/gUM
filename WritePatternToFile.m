function WritePatternToFile(FN,M,N)
% M=randi(57,3,4);
% FN='Pattern2.txt';
% Ttl='Pattern 2, 11 patterns for 120';
fid=fopen(FN,'wt');
% fprintf(fid,'%s\n',Ttl);
fprintf(fid,'%d %d %d %d\n',size(M,1),size(M,2),size(M,3),N);
for i=1:size(M,1)
    for j=1:size(M,2)
        fprintf(fid,' %d',M(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%
% FN='Pattern2.txt';
% fid=fopen(FN,'rt');
% % TtlX=fscanf(fid,'%s\n');
% Sz=fscanf(fid,'%d %d %d\n',3);
% for i=1:size(M,1)
%     for j=1:size(M,2)
%         MX(i,j)=fscanf(fid,' %d',1);
%     end
%     fscanf(fid,'\n');
% end
% fclose(fid);