A=load('/media/a/H1/HCPData_256x256_int16.mat');
I=permute(A.HCPData(:,:,1:10000),[3 1 2]);
save('/media/a/H1/First10kIm256x256Magint16.mat','I');
disp('10')
I=permute(A.HCPData(:,:,1:20000),[3 1 2]);
save('/media/a/H1/First20kIm256x256Magint16.mat','I','-v7.3');
disp('20')
I=permute(A.HCPData(:,:,:),[3 1 2]);
save('/media/a/H1/AllIm256x256Magint16.mat','I','-v7.3');
disp('All')