function gmkdir(P)
if(~exist(P,'dir'))
    mkdir(P);
end
system(['chmod +777 -R ' P]);
disp(['Created ' P]);