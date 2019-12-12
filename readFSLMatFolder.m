function Out=readFSLMatFolder(FNP);
% FNP='Movingx_mcf.mat';
D=dir([FNP filesep 'MAT_*']);
Out=zeros(4,4,numel(D));
for i=1:numel(D)
    Out(:,:,i)=readFSLMatFile([FNP filesep D(i).name]);
end