load('Maastricht_Spiral.mat');

SensP=ones(Sz1);

for c=1:NCha
    disp(c);
    DataP=conj(projC(:,:,:,c));
    recoC(:,:,c) = bart('pics -r:0.0001 -R T:7:0:0.1 -t ',TrajectoryP, DataP, SensP);
end

ESPIRIT_val=15;
FF=bart('fft 6',permute(recoC,[4 1 2 3]));
calib = bart(['ecalib -r ' num2str(ESPIRIT_val)], FF(1,:,:,:));
SensB(:,:,:,1) = permute(bart('slice 4 0', calib),[2 3 4 1]);
SensP=permute(SensB,[1 2 4 3]);

DataP=conj(projC);
PICSprm='-m -r:0.0001 -R W:7:0:0.000001 -t';
recoX = bart(['pics ' PICSprm],TrajectoryP, DataP, SensP);
figure;imagesc(abs(recoX));colormap gray;
title(PICSprm);axis equal