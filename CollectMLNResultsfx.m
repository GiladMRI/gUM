function [MLNRXS, SlicesFound]=CollectMLNResultsfx(BaseTFP,MLNPrefix,nSlices)
% get MLN result
clear MLNRXS
SlicesFound=[];
for SliI=1:nSlices
    try
%     BaseTFP='/autofs/cluster/kawin/Gilad/TF/srezN/';
    DD=dir(BaseTFP);
    DD=DD(strhas({DD.name},'train'));
    DD=DD(strhas({DD.name},[MLNPrefix num2str(SliI) '_']));
    CurSRP=[BaseTFP DD( ([DD.datenum]==max([DD.datenum])) ).name filesep];
    DD2=dir(CurSRP);
    DD2=DD2(strhas({DD2.name},'png'));
    [DDS DDOrd]=sort([DD2.datenum]);
    
    MLNR=imread([CurSRP DD2(DDOrd(end)).name]);
    MLNSz=gsize(MLNR,1:2);
    MLNSz(1)=MLNSz(1)/8;
    MLNSz(2)=MLNSz(2)/(8*8);
    MLNRX=MLNR(1:MLNSz(1),MLNSz(2)*8*5+(1:(MLNSz(2)*8)),1);
    MLNRX=PartitionDim(MLNRX,2,8);
    MLNRXS(:,:,:,SliI)=MLNRX;
    disp(SliI);
    SlicesFound=[SlicesFound SliI];
    catch
    end
end
MLNRXS(1,1,1,nSlices)=0;
MLNRXS=double(MLNRXS);
% MLNRXS=MLNRXS(:,:,:,ROrd);
% SlicesFound=sort(Ord(SlicesFound));