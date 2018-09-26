% opengl save software

clear

UserName=char(java.net.InetAddress.getLocalHost.getHostName);
if(strhas(UserName,'deni'))
    BaseTFFolder='~/TF/';
    TFFolder='/home/deni/TF/srez/';
    cd('/home/deni/gUM/Wavelab850/');
else % gilad
    TFFolder='/home/a/HomeA/TF/srez/';
    cd([BasePA 'gUM/Wavelab850/']);
end

% BasePA='~/HomeA/';
BasePA='~/';

WavePath;

run ([BasePA 'irt/setup.m']);

cd([BasePA 'gUM']);
BasePX=[pwd filesep];
BaseP=[BasePX '4GL/'];

addpath([BaseP 'Spiral_recon_T1/'])
addpath(genpath([BaseP 'Spiral_recon_T1/raw_header']))
addpath([BaseP 'Spiral_recon_T1/misc'])
addpath([BaseP 'Spiral_recon_T1/io'])

addpath(genpath([BasePA 'SPENOnline']))

addpath(genpath([BasePA 'gUM/EPFLSpiral']))

addpath(genpath([BasePA 'gpuNUFFT-master']))

addpath(genpath([BasePA 'Tools/NIfTI_20140122']));

addpath([BasePA 'gUM']);

addpath(genpath([BasePA 'gUM/ASPIRE-master']))

setenv('TOOLBOX_PATH',[BasePA 'bart-0.4.03'])

addpath([BasePA 'gUM/NIFTI_20140122/']);

