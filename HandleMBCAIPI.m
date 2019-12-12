nukData=ADataIsPy(:,:,SliI,1).';
nukData=nukData(:,3:end);
%     nukData=nukData.*T2SEstComp;
nukDataP=permute(nukData,[3 2 4 1]);
SensP=permute(SensX(:,:,:,:,SliI),[1 2 5 3 4]);

Out=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],BARTTrajAct, nukDataP, SensP(:,:,:,:,1));
%%
SensP=permute(SensX(:,:,:,1,[SliI SliI+12]),[1 2 5 3 4]); % 128   128     2    32
SensP=permute(SensP,[1 2 5 4 3]);
Out2=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],BARTTrajAct, nukDataP, SensP);
%%
for SliI=1:12
    nukData=ADataIsPy(:,:,SliI,1).';
    nukData=nukData(:,3:end);
    %     nukData=nukData.*T2SEstComp;
    nukDataP=permute(nukData,[3 2 4 1]);

    SensP=permute(SensX(:,:,:,1,[SliI SliI+12]),[1 2 5 3 4]); % 128   128     2    32
    SensP=permute(SensP,[1 2 5 4 3]);
    OutMS(:,:,SliI,:,:)=bart(['pics -S -m -R W:7:0:' num2str(1e-5) ' -t'],BARTTrajAct, nukDataP, SensP);
end
OutMSX=CombineDims(OutMS,[4 3]);
%%
fgmontage(OutMSX)

fgmontage(permute(OutMSX,[1 3 2]))

fgmontage(permute(OutMSX,[3 1 2]))