% MMM=repmat(Rec_CompgB0_RSS_MXCMY(:,:,SliIs,5:8),[1 1 5 1]);
% RRR=repmat(RefSV1x(:,:,Ord(SliIs)),[1 1 5 1]);

MMM=repmat(ToMove,[1 1 5 1]);
RRR=repmat(ForRef,[1 1 5 1]);

MovFN='Movingx';
RefFN='Ref';

delete([MovFN '.nii']);
delete([RefFN '.nii']);

Raw2Nii(MMM,[MovFN '.nii'],'float32');
Raw2Nii(RRR,[RefFN '.nii'],'float32');

CostTypeMCF='normcorr';
% CostTypeMCF='normmi';
CostType='normcorr';

system(['rm -r ' MovFN '_mcf.mat*']);
%% First mcflirt into median volume
MCF_FN=[MovFN '_mcf.nii.gz'];
delete(MCF_FN);
Extra='';

Cmd=['mcflirt -in ' MovFN '.nii -cost ' CostTypeMCF ' -dof 6 -mats -refvol 1 '];
system(Cmd);
% MCF=loadniidata(MCF_FN);
% % fgmontagex(MCF(:,:,3,:));title('MCF')
% fgmontagex(MCF(:,:,3,1));title('MCF1')
% fgmontagex(MCF(:,:,3,2));title('MCF2')
% fgmontagex(MCF(:,:,3,3));title('MCF3')
% fgmontagex(MCF(:,:,3,4));title('MCF4')
% 
% fgmontagex(ToMove(:,:,1));title('ToMove1')
% fgmontagex(ToMove(:,:,2));title('ToMove2')
%% Average across time
TmeanFN=[MovFN '_mcf_mean.nii.gz'];
delete(TmeanFN);

Cmd=['fslmaths ' MovFN '_mcf.nii.gz -Tmean ' MovFN '_mcf_mean.nii.gz'];
system(Cmd);
% Tmean=loadniidata(TmeanFN);
%% flirt on Tmean
TmeanMoved_FN=[MovFN '_mcf_mean_to_' RefFN '.nii.gz'];
delete(TmeanMoved_FN);

Cmd=['flirt -in ' MovFN '_mcf_mean.nii.gz -ref ' RefFN,...
    '.nii -cost ' CostType ' -2D -out ' MovFN '_mcf_mean_to_' RefFN '.nii -omat ' MovFN '_mcf_mean_to_' RefFN '.mat'];% ...
%     ' -omat ' MovFN '.mat'];
system(Cmd);
% TmeanMoved=loadniidata(TmeanMoved_FN);
%% Apply second transormations
Moved_FN=[MovFN '_mcf_to_' RefFN '.nii.gz'];
delete(Moved_FN);
Cmd=['flirt -in ' MovFN '_mcf.nii.gz -ref ' RefFN '.nii -applyxfm -init ' MovFN '_mcf_mean_to_' RefFN '.mat -out ',...
    MovFN '_mcf_to_' RefFN '.nii'];
system(Cmd);
%%
MMMC=loadniidata(Moved_FN);
Moved=MMMC(:,:,3,:);
%%
% OutFN=['Correctedx_' CostType Extra(2:end)];
% delete([OutFN '.nii']);
% delete([OutFN '.nii.gz']);
% Cmd=['/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/6.0.1/bin/mcflirt -in ', MovFN,...
%     ' -reffile ', RefFN,...
%     ' -out ' OutFN ' ',...
%     ' -bins 256 ',...
%     ' -cost ' CostType,...
%     ' -dof 6 ',Extra,...
%     ' -mats -plots '];
%       
% system(Cmd);
% MMMC=loadniidata([OutFN '.nii.gz']);
% disp('Ran mclirt');
% %%
% 
% % /autofs/cluster/kawin/Gilad/Correctedx.nii
% % /autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/6.0.1/bin/flirt -in /autofs/cluster/kawin/Gilad/Moving.nii -ref /autofs/cluster/kawin/Gilad/First.nii -out /autofs/cluster/kawin/Gilad/Corrected.nii -omat /autofs/cluster/kawin/Gilad/Corrected.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -2D -dof 12  -interp trilinear
% % -cost {mutualinfo,woods,corratio,normcorr,normmi,leastsquares} : default is normcorr
% % http://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/MCFLIRT.html
% % -smooth <num>   
% CostType='normcorr';
% % CostType='normcorr';
% Extra='';
% % Extra='-gdt';
% % /autofs/cluster/kawin/Gilad/
% OutFN=['Correctedx_' CostType Extra(2:end)];
% delete([OutFN '.nii']);
% delete([OutFN '.nii.gz']);
% Cmd=['/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/6.0.1/bin/mcflirt -in ', MovFN,...
%     ' -reffile ', RefFN,...
%     ' -out ' OutFN ' ',...
%     ' -bins 256 ',...
%     ' -cost ' CostType,...
%     ' -dof 6 ',Extra,...
%     ' -mats -plots '];
%       
% system(Cmd);
% MMMC=loadniidata([OutFN '.nii.gz']);
% disp('Ran mclirt');
% 
% RRR4=CombineDims(CombineDims(PartitionDim(RRR(:,:,1:4),3,2),[3 1]),[3 2]);
% MMM4=CombineDims(CombineDims(PartitionDim(squeeze(MMM(:,:,3,1:4)),3,2),[3 1]),[3 2]);
% MMMC4=CombineDims(CombineDims(PartitionDim(squeeze(MMMC(:,:,3,1:4)),3,2),[3 1]),[3 2]);
% RRR4=RRR4./grmss(RRR4);
% MMM4=MMM4./grmss(MMM4);
% MMMC4=MMMC4./grmss(MMMC4);
% 
% % figure;imshowpair(MMM4,RRR4,'Checkerboard');
% % figure;imshowpair(MMMC4,RRR4,'Checkerboard');title([CostType ' ' Extra]);
% % fgmontagex(RRR4);title('Ref');
% fgmontagex(MMMC4);title([CostType ' ' Extra]);