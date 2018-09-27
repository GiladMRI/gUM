BaseP='C:\STRDCE\John\Database\DCEOut\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
D={D.name}';
%%
for i=1:numel(D)
    disp([i numel(D)]);
    WorkingP=[BaseP D{i} filesep];
    MakeReport;
end
disp('Finished');