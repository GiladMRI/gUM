setenv('TOOLBOX_PATH','~/HomeA/bart-0.4.03MultiKernel')

% a=magic(3)+1i*randn(3);
% a=repmat((1:9).',[1 1 7]);
% clear Kernels
% Kernels(:,:,1)=[-1 1 0; -1 1 0; -1 1 0];
% Kernels(:,:,2)=Kernels(:,:,1).';

% Kernels=DD(:,:,1:18);
Kernels=DD;
% Kernels=randn([3 3 2]);

Kernels=Kernels*0.;

writecfl('/tmp/Kernels',Kernels);
Rec=bart(['pics -S -m -R T:3:0:' num2str(1e-4) ' -p'],Q, FI.*Q, FI*0+1);

fgmontage(Rec,[0 130]);
%%
% b=readcfl('/tmp/qqq');