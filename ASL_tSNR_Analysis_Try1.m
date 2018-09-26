addpath('/home/deni/gUM/NIFTI_20140122/');
M=rot90(loadniidata('/media/deni/c78a9273-3214-4387-9f72-4cdc3adef255/OnBP_9Aug18/BP_EPI_ASL_DCM/20180809_154152BPFAIR2x2z3RSs012a001.nii'));
%% Very simple Perfusion tSNR analysis
MEven=M(:,:,:,2:2:end);
MOdd=M(:,:,:,1:2:end);
MPerf=MEven-MOdd;
MeanMPerf=mean(MPerf,4);
StdMPerf=std(MPerf,0,4);
Perf_tSNR=MeanMPerf./StdMPerf;

fgmontage(Perf_tSNR,[0 1])

% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide
FN='asd';
Cmd=['fsl fnirt ' FN ' asd'];
system(Cmd);


