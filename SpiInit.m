% opengl save software

clear

BasePA='/media/a/ec52f4a8-12c2-4a24-9e6f-d65aaa410529/Gilad/';
BasePA='/media/g/301A1F881A1F49EC/';
BasePA='/autofs/space/daisy_002/users/Gilad/';
% cd('~/HomeA/gUM/Wavelab850');
% WavePath;

run ([BasePA 'irt/setup.m']);

cd([BasePA 'gUM']);
BasePX=[pwd filesep];
BaseP=[BasePX '4GL/'];

addpath([BaseP 'Spiral_recon_T1/'])
addpath(genpath([BaseP 'Spiral_recon_T1/raw_header']))
addpath([BaseP 'Spiral_recon_T1/misc'])
addpath([BaseP 'Spiral_recon_T1/io'])

addpath(genpath([BasePA 'SPENOnline']))
addpath(genpath(['/media/g/aa384808-2266-4edb-87a7-637bc772dc2f/SPENOnline']));

addpath(genpath([BasePA 'gUM/EPFLSpiral']))

addpath(genpath([BasePA 'gpuNUFFT-master']))

addpath(genpath([BasePA 'Tools/NIfTI_20140122']));

addpath([BasePA 'gUM']);

addpath(genpath([BasePA 'gUM/ASPIRE-master']))

setenv('TOOLBOX_PATH',[BasePA 'bart-0.4.03'])

!synclient HorizTwoFingerScroll=0

addpath('/autofs/space/daisy_002/users/Gilad/t2shuffling-support-master/src/utils');

addpath(BasePA);

% setenv('TOOLBOX_PATH','/autofs/cluster/kawin/sid/bart')
setenv('TOOLBOX_PATH','/autofs/space/daisy_002/users/Gilad/bart-0.4.04b')
setenv('SHELL','/bin/bash')
setenv('MATLAB_SHELL','/bin/bash')
system('/bin/bash source /usr/pubsw/packages/mkl/2019/bin/compilervars.sh intel64')
setenv('LD_LIBRARY_PATH','/usr/pubsw/packages/mkl/2019/compilers_and_libraries_2019.0.117/linux/compiler/lib/intel64_lin:/autofs/cluster/pubsw/arch/CentOS7-x86_64/packages/mkl/2019//compilers_and_libraries_2019.0.117/linux/mpi/intel64/lib/release:/autofs/cluster/pubsw/arch/CentOS7-x86_64/packages/mkl/2019//compilers_and_libraries_2019.0.117/linux/mpi/intel64/lib:/usr/pubsw/packages/mkl/2019/compilers_and_libraries_2019.0.117/linux/compiler/lib/intel64_lin:/usr/pubsw/packages/mkl/2019/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin:/autofs/cluster/pubsw/arch/CentOS7-x86_64/packages/mkl/2019/compilers_and_libraries_2019.0.117/linux/tbb/lib/intel64/gcc4.7:/autofs/cluster/pubsw/arch/CentOS7-x86_64/packages/mkl/2019/compilers_and_libraries_2019.0.117/linux/tbb/lib/intel64/gcc4.7')

addpath([BasePA 'safe_pns_prediction-master']);

PhaseCyclingBaseP=[BasePA 'phase_cyclingD/'];
addpath(PhaseCyclingBaseP);
addpath([PhaseCyclingBaseP 'utils'])
addpath([PhaseCyclingBaseP 'data']);
addpath([PhaseCyclingBaseP 'linops']);
addpath([PhaseCyclingBaseP 'prox']);
addpath([PhaseCyclingBaseP 'phase_cycling']);

addpath('/autofs/space/daisy_002/users/Gilad/gUM/tOptGrad_V0.2/minTimeGradient/Matlab');
addpath(genpath('/autofs/space/daisy_002/users/Gilad/mintgrad/'));

PWD=pwd;
cd('/autofs/space/daisy_002/users/Gilad/cvx')
cvx_startup
cd(PWD)

addpath([BasePA 'gUM/gif/']);

% getenv('HOME')
% char(java.net.InetAddress.getLocalHost.getHostName)

% addpath(genpath('/media/a/DATA/StreamComplexImages - binaries, Matlab, everthing/MatlabCode/'))