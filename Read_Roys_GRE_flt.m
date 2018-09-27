BaseP='/media/a/DATA/Roys_GRE_ME_07_14Subjects/Recond/';
D=dir([BaseP '*.dat']);
D=D([D.bytes]>1e10);
DNames={D.name}';
%%
for i=13:numel(DNames)
    for e=1:4
        disp([i e]);
        disp(DNames{i});
        [ data_matrix, filename ]  = read_flt_Roy([BaseP DNames{i}], e);
        data_matrix=permute(data_matrix,[1 3 4 2]);
        save([BaseP DNames{i}(1:end-4) '_E' num2str(e) '.mat'],'data_matrix');
    end
end
%%
fid=fopen(FN,'rb');

XX = fread(fid, 1,'uint32');

[prot,rstraj] = read_twix_hdr(fid);